#!/usr/bin/env python3
"""
pep2rna.py — Detect peptides encoded by annotated ORFs in RNA-seq reads.

PURPOSE
-------
Given a list of peptides (typically derived from a target protein or a
peptide-array experiment), search RNA-seq BAM files for reads whose
translated sequence (in any of the three reading frames) contains one of
those peptides. This is used to look for evidence that a peptide is
actually expressed at the RNA level, based on annotated nuORF and
MetamORF loci.

METHOD
------
1. Read a TSV mapping peptides -> target identifiers (e.g. AltP43s).
2. Build an Aho-Corasick automaton over the peptide patterns.
3. Scan two reference FASTA files of annotated ORFs (nuORFdb, MetamORF)
   to find, for each peptide, the genomic loci (chrom, start, end, strand)
   where it appears.
4. Merge loci coming from both sources and collapse overlapping intervals
   per (chrom, strand).
5. For each BAM file, for each peptide, for each locus, fetch the reads
   overlapping the locus, filter by strand and MAPQ, reconstruct the
   read sequence in mRNA orientation, translate it in all three frames,
   and count how many reads contain the peptide.
6. Report per-BAM, per-locus hit counts.

INPUTS
------
Required:
  --peptides      TSV file, tab-separated, with header. Columns:
                    col 1: peptide sequence
                    col 2: target identifier (e.g. AltP43s)
  --nuorf         FASTA of nuORFdb protein sequences. Headers are expected
                  to look like:  <core>:<start>-<end>:<strand>|<...>
  --metamorf      FASTA of MetamORF protein sequences. Headers are expected
                  to contain:  chr=<chrom>  strand=<+|->
                                start_pos=<int>  stop_pos=<int>
  --bam-dir       Directory containing the .bam files to scan.

Optional:
  --pattern-ids   Comma-separated subset of identifiers to keep.
                  Default: the built-in "short" set (AltP43s, AltP11s, ...).
  --output        Output TSV path. Default: hits.tsv
  --workers       Number of parallel worker processes. Default: 28
  --min-mapq      Minimum mapping quality for a read to be considered.
                  Default: 0
  --quiet         Do not echo the report to stdout.

OUTPUTS
-------
A TSV file (path given by --output) with one row per (bam, peptide, locus)
that had at least one hit. Columns:

  bam_file    peptide    chrom    start    end    strand    reads    hits

Plus a final line on stdout summarizing the total number of hits.

USAGE
-----
  python pep2rna.py \\
      --peptides  peptides.tsv \\
      --nuorf     /path/PA_nuORFdb_v1.2_protein.fasta \\
      --metamorf  /path/MetamORF_Hsapiens_aa.fasta \\
      --bam-dir   /mnt/pi_raid/.../ANALYSIS/BAM/ \\
      --workers   28 \\
      --output    hits.tsv
"""

import argparse
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
import ahocorasick
import pysam
from Bio.Data import CodonTable


# --------------------------------------------------------------------------- #
# Defaults                                                                    #
# --------------------------------------------------------------------------- #
DEFAULT_SHORT_IDS = {
    "AltP43s", "AltP11s", "AltP59s", "AltP42s", "AltP13s",
    "AltP53s", "AltP19s", "AltP39s", "AltP23s", "AltP34s",
}


# --------------------------------------------------------------------------- #
# Translation helpers                                                         #
# --------------------------------------------------------------------------- #
_T = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = {**_T.forward_table, **{s: "*" for s in _T.stop_codons}}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(nt):
    """Translate a nucleotide string, ignoring a partial trailing codon."""
    return "".join(CODON_TO_AA.get(nt[i:i + 3], "X")
                   for i in range(0, len(nt) - 2, 3))


def read_mrna(read):
    """Transcript-space nucleotide sequence of the read, in mRNA orientation."""
    qs = read.query_sequence
    if qs is None:
        return None
    pairs = [(q, r) for q, r in read.get_aligned_pairs()
             if q is not None and r is not None]
    if read.is_reverse:
        pairs = pairs[::-1]           # walk 5'->3' along the mRNA
        return "".join(qs[q].translate(_COMP) for q, r in pairs)
    return "".join(qs[q] for q, r in pairs)


# --------------------------------------------------------------------------- #
# Building the peptide -> loci map                                            #
# --------------------------------------------------------------------------- #
def build_pattern_to_ids(peptides_tsv, keep_ids):
    """Read the peptide TSV, return {peptide: {id, ...}} for kept ids."""
    pattern_to_ids = {}
    with open(peptides_tsv) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            peptide, id_ = fields[0].strip(), fields[1].strip()
            if id_ not in keep_ids:
                continue
            pattern_to_ids.setdefault(peptide, set()).add(id_)
    return pattern_to_ids


def build_automaton(pattern_to_ids):
    A = ahocorasick.Automaton()
    for p in pattern_to_ids:
        A.add_word(p, p)
    A.make_automaton()
    return A


def scan_metamorf(path, A, pattern_to_ids):
    results = {p: {"ids": sorted(pattern_to_ids[p]), "loci": []}
               for p in pattern_to_ids}
    seen = {p: set() for p in pattern_to_ids}

    with open(path) as f:
        header = ""
        for line in f:
            if line.startswith(">"):
                header = line.rstrip("\n")[1:]
                continue
            seq = line.rstrip("\n")
            for _, p in A.iter(seq):
                if header in seen[p]:
                    continue
                seen[p].add(header)
                chrom = re.search(r"\bchr=(\S+)", header).group(1)
                strand = re.search(r"\bstrand=(\S+)", header).group(1)
                start = re.search(r"\bstart_pos=(\d+)", header).group(1)
                end = re.search(r"\bstop_pos=(\d+)", header).group(1)
                results[p]["loci"].append((chrom, start, end, strand))
    return results


def scan_nuorf(path, A, pattern_to_ids):
    results = {p: {"ids": sorted(pattern_to_ids[p]), "loci": []}
               for p in pattern_to_ids}
    seen = {p: set() for p in pattern_to_ids}

    with open(path) as f:
        header = ""
        for line in f:
            if line.startswith(">"):
                header = line.rstrip("\n")[1:]
                continue
            seq = line.rstrip("\n")
            for _, p in A.iter(seq):
                if header in seen[p]:
                    continue
                seen[p].add(header)
                core = header.split("|")[0]
                loc = core.rsplit(":", 2)
                chrom = loc[0].rsplit("_", 1)[-1]
                start, end = loc[1].split("-")
                strand = loc[2]
                results[p]["loci"].append((chrom, start, end, strand))
    return results


def merge_databases(*dicts):
    """Merge {peptide: {ids, loci}} dicts, deduplicating ids and loci."""
    merged_ids = defaultdict(set)
    merged_loci = defaultdict(dict)          # peptide -> {tuple: locus}
    for d in dicts:
        for peptide, info in d.items():
            merged_ids[peptide].update(info.get("ids", []))
            for locus in info.get("loci", []):
                merged_loci[peptide].setdefault(tuple(locus), locus)
    return {
        p: {"ids": sorted(merged_ids[p]),
            "loci": list(merged_loci[p].values())}
        for p in merged_ids
    }


def merge_overlapping_loci(peptides):
    """In-place: collapse overlapping intervals per (chrom, strand)."""
    for info in peptides.values():
        by_key = defaultdict(list)
        for chrom, start, end, strand in info["loci"]:
            by_key[(chrom, strand)].append((int(start), int(end)))

        merged = []
        for (chrom, strand), intervals in by_key.items():
            intervals.sort()
            s0, e0 = intervals[0]
            for s, e in intervals[1:]:
                if s <= e0:
                    e0 = max(e0, e)
                else:
                    merged.append((chrom, str(s0), str(e0), strand))
                    s0, e0 = s, e
            merged.append((chrom, str(s0), str(e0), strand))
        info["loci"] = merged


# --------------------------------------------------------------------------- #
# Per-BAM scan                                                                #
# --------------------------------------------------------------------------- #
def process_bam(bam_file, bam_dir, peptides, min_mapq):
    """Scan one BAM against all peptides/loci. Returns (bam_file, hits, rows)."""
    bam_path = os.path.join(bam_dir, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")

    local_hits = 0
    rows = []
    try:
        for peptide, info in peptides.items():
            if not info["loci"]:
                continue
            for chrom, start, end, strand in info["loci"]:
                chrom_name = "chr" + chrom
                start_i, end_i = int(start), int(end)
                n_reads = n_hits = 0

                for read in bam.fetch(chrom_name, start_i - 1, end_i):
                    if read.mapping_quality < min_mapq:
                        continue
                    if strand == "-" and not read.is_reverse:
                        continue
                    if strand == "+" and read.is_reverse:
                        continue

                    mrna = read_mrna(read)
                    if mrna is None or len(mrna) < 3 * len(peptide):
                        continue

                    n_reads += 1
                    for frame in range(3):
                        if peptide in translate(mrna[frame:]):
                            n_hits += 1
                            break

                if n_hits > 0:
                    local_hits += n_hits
                    rows.append((bam_file, peptide, chrom_name,
                                 start_i, end_i, strand, n_reads, n_hits))
    finally:
        bam.close()

    return bam_file, local_hits, rows


# --------------------------------------------------------------------------- #
# CLI                                                                         #
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Search RNA-seq BAMs for reads whose translation "
                    "contains any of a set of target peptides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--peptides", required=True,
                   help="TSV file: peptide <tab> identifier (with header).")
    p.add_argument("--nuorf", required=True,
                   help="nuORFdb protein FASTA.")
    p.add_argument("--metamorf", required=True,
                   help="MetamORF protein FASTA.")
    p.add_argument("--bam-dir", required=True,
                   help="Directory containing .bam files.")
    p.add_argument("--pattern-ids", default=None,
                   help="Comma-separated subset of identifiers to keep. "
                        "If omitted, the built-in 'short' set is used.")
    p.add_argument("--output", default="hits.tsv",
                   help="Output TSV path.")
    p.add_argument("--workers", type=int, default=28,
                   help="Number of parallel worker processes.")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Minimum mapping quality for a read to be kept.")
    p.add_argument("--quiet", action="store_true",
                   help="Do not echo the report to stdout.")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #
def main(argv=None):
    args = parse_args(argv)

    keep_ids = (set(args.pattern_ids.split(","))
                if args.pattern_ids else DEFAULT_SHORT_IDS)

    # 1) peptides -> ids
    pattern_to_ids = build_pattern_to_ids(args.peptides, keep_ids)
    print(f"[info] {len(pattern_to_ids)} distinct peptide patterns "
          f"for {len(set().union(*pattern_to_ids.values()))} ids",
          file=sys.stderr)

    # 2) automaton
    A = build_automaton(pattern_to_ids)

    # 3) scan references
    metamorf_results = scan_metamorf(args.metamorf, A, pattern_to_ids)
    nuorf_results = scan_nuorf(args.nuorf, A, pattern_to_ids)

    # 4) merge + collapse
    peptides = merge_databases(nuorf_results, metamorf_results)
    merge_overlapping_loci(peptides)

    n_with_loci = sum(1 for info in peptides.values() if info["loci"])
    print(f"[info] {n_with_loci} / {len(peptides)} peptides have loci",
          file=sys.stderr)

    # 5) BAM files
    bam_files = [x for x in os.listdir(args.bam_dir) if x.endswith("bam")]
    print(f"[info] {len(bam_files)} BAM files to scan",
          file=sys.stderr)

    # 6) parallel scan
    total_hits = 0
    out_rows = []
    tasks = [(bf, args.bam_dir, peptides, args.min_mapq) for bf in bam_files]
    with Pool(processes=args.workers) as pool:
        for bam_file, hits, rows in pool.starmap(process_bam, tasks,
                                                 chunksize=1):
            total_hits += hits
            out_rows.extend(rows)

    # 7) output
    with open(args.output, "w") as f:
        f.write("bam_file\tpeptide\tchrom\tstart\tend\tstrand\treads\thits\n")
        for row in out_rows:
            f.write("\t".join(map(str, row)) + "\n")

    if not args.quiet:
        for row in out_rows:
            print("\t".join(map(str, row)))

    print(f"Total Hits: {total_hits}", file=sys.stderr)
    print(f"[info] wrote {len(out_rows)} rows to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
