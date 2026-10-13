import os
import re
import pysam
import ahocorasick
from Bio.Data import CodonTable
from multiprocessing import Pool
from collections import Counter
from collections import defaultdict

short = {"AltP43s", "AltP11s", "AltP59s", "AltP42s", "AltP13s",
         "AltP53s", "AltP19s", "AltP39s", "AltP23s", "AltP34s"}

long  = {"AltP26L", "AltP25L", "AltP43L", "AltP59L", "AltP42L",
         "AltP13L", "AltP53L", "AltP19L", "AltP39L", "AltP40L",
         "AltP23L", "AltP34L"}

nuORF = "/home/dimitris/Workspace/reference/PA_nuORFdb_v1.2_protein.fasta"
MetamORF = "/home/dimitris/Workspace/reference/MetamORF_Hsapiens_aa.fasta"

# ---- 1. build the pattern -> ids map ---------------------------------------
pattern_to_ids = {}   # p -> set of ids that produced p
with open('peptides.tsv') as f:
    next(f)
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 2:
            continue
        peptide, id_ = fields[0].strip(), fields[1].strip()
        if id_ in short:
            pats = [peptide]
#        elif id_ in long:
#            pats = [peptide[i:i+8] for i in range(len(peptide) - 7)]
        else:
            continue

        for p in pats:
            pattern_to_ids.setdefault(p, set()).add(id_)

# ---- 2. build the automaton ------------------------------------------------
A = ahocorasick.Automaton()
for p in pattern_to_ids:
    A.add_word(p, p)
A.make_automaton()

# ---- 3. scan the nuORF, filling the results dict ---------------------------
MetamORF_results = {p: {"ids": sorted(pattern_to_ids[p]), "loci": []}
           for p in pattern_to_ids}

seen = {p: set() for p in pattern_to_ids}   # dedupe loci per peptide

with open(MetamORF) as f:
    header = ''
    for line in f:
        if line.startswith('>'):
            header = line.rstrip('\n')[1:]   # drop the leading '>'
            continue
        seq = line.rstrip('\n')
        for _, p in A.iter(seq):
            if header not in seen[p]:
                seen[p].add(header)
                chrom  = re.search(r'\bchr=(\S+)',       header).group(1)
                strand = re.search(r'\bstrand=(\S+)',    header).group(1)
                start  = re.search(r'\bstart_pos=(\d+)', header).group(1)
                end    = re.search(r'\bstop_pos=(\d+)',  header).group(1)
                MetamORF_results[p]["loci"].append((chrom, start,end,strand))

# ---- 4. scan the nuORF, filling the results dict ---------------------------
nuORF_results = {p: {"ids": sorted(pattern_to_ids[p]), "loci": []}
           for p in pattern_to_ids}

seen = {p: set() for p in pattern_to_ids}   # dedupe loci per peptide

with open(nuORF) as f:
    header = ''
    for line in f:
        if line.startswith('>'):
            header = line.rstrip('\n')[1:]   # drop the leading '>'
            continue
        seq = line.rstrip('\n')
        for _, p in A.iter(seq):
            if header not in seen[p]:
                seen[p].add(header)
                core   = header.split('|')[0]
                loc    = core.rsplit(':', 2)
                chrom  = loc[0].rsplit('_', 1)[-1]
                start, end = loc[1].split('-')
                strand = loc[2]
                nuORF_results[p]["loci"].append((chrom, start,end,strand))

# -------- merge databases ------------------------
_merged_ids  = defaultdict(set)
_merged_loci = defaultdict(dict)   # peptide -> {tuple(locus): locus}

for d in (nuORF_results, MetamORF_results):
    for peptide, info in d.items():
        _merged_ids[peptide].update(info.get("ids", []))
        for locus in info.get("loci", []):
            _merged_loci[peptide].setdefault(tuple(locus), locus)

peptides = {
    p: {"ids": sorted(_merged_ids[p]),
        "loci": list(_merged_loci[p].values())}
    for p in _merged_ids
}
# ---- merge overlaping loci --------------
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

#for peptide,info in peptides.items():
#    print(peptide,info['loci'])
# ---- translation -----------------------------------------------------------
_T = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = {**_T.forward_table, **{s: "*" for s in _T.stop_codons}}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(nt):
    return "".join(CODON_TO_AA.get(nt[i:i+3], "X")
                   for i in range(0, len(nt) - 2, 3))


def read_mrna(read):
    """Transcript-space nucleotide sequence of the read, in mRNA orientation."""
    qs = read.query_sequence
    if qs is None:
        return None
    pairs = [(q, r) for q, r in read.get_aligned_pairs()
             if q is not None and r is not None]
    if read.is_reverse:
        pairs = pairs[::-1]      # walk 5'->3' along the mRNA
        return "".join(qs[q].translate(_COMP) for q, r in pairs)
    return "".join(qs[q] for q, r in pairs)

def get_bam(path):
    if path not in _bam_cache:
        _bam_cache[path] = pysam.AlignmentFile(path, "rb")
    return _bam_cache[path]


def process_bam(bam_file, peptides):
    """Process a single BAM file across all peptides/loci.

    `peptides` is a plain dict passed in (picklable)."""

    bam_path = os.path.join(BASE_PATH, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")

    local_hits = 0
    per_bam_report = []

    for PEPTIDE, info in peptides.items():
        if not info["loci"]:
            continue
        for locus in info["loci"]:
            CHROM = "chr" + locus[0]
            START = int(locus[1])
            END = int(locus[2])
            STRAND = locus[3]
            MIN_MAPQ = 0
            n_reads = 0
            n_hits = 0

            for read in bam.fetch(CHROM, START - 1, END):
                if read.mapping_quality < MIN_MAPQ:
                    continue
                if STRAND == "-" and not read.is_reverse:
                    continue
                if STRAND == "+" and read.is_reverse:
                    continue

                mrna = read_mrna(read)
                if mrna is None or len(mrna) < 3 * len(PEPTIDE):
                    continue

                n_reads += 1

                for frame in range(3):
                    if PEPTIDE in translate(mrna[frame:]):
                        n_hits += 1
                        break

            if n_hits > 0:
                local_hits += n_hits
                per_bam_report.append(
                    f"{bam_file}  {CHROM}:{START}-{END} ({STRAND})  "
                    f"reads={n_reads}  hits={n_hits}  peptide={PEPTIDE}"
                )

    bam.close()
    return bam_file, local_hits, per_bam_report

# --- Globals that need to be visible in workers ---
#BASE_PATH = "/mnt/pi_raid/raid/PJ2410297_20250115/ANALYSIS/BAM/" # 142 samples
BASE_PATH = "/mnt/pi_raid/raid/PJ2406160_20240807_Integragen_Tessa/ANALYSIS/BAM/" # 43 samples 
#BASE_PATH = "/mnt/pi_raid/raid/PJ2406160_20240909_Integragen_Tessa/ANALYSIS/BAM/" # 5 samples
_bam_cache = {}

def main():
    contents = os.listdir(BASE_PATH)
    bam_files = [x for x in contents if x.endswith("bam")]
    total_hits = 0
    with Pool(processes=28) as pool:
        for bam_file, hits, report in pool.starmap(
            process_bam,
            [(bf, peptides) for bf in bam_files],
            chunksize=1,
        ):
            total_hits += hits
            for line in report:
                print(line)

    print("Total Hits:", total_hits)


if __name__ == "__main__":
    main()
