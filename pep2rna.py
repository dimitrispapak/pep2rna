import multiprocessing as mp
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, List, Tuple

import pysam
from Bio.Data import CodonTable
from Bio.Seq import Seq

"""
pep2rna.py — Find genomic loci encoding a peptide and count supporting reads.

Parallel, load-balanced version:
  * Genome scan: fixed-size windows with peptide-length overlap.
  * Read counting: loci clustered so overlapping ones share a worker.
  * Per-window mapped-read counting for a matched CPM denominator.
"""



# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------
_TABLE = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = dict(_TABLE.forward_table)
for _s in _TABLE.stop_codons:
    CODON_TO_AA[_s] = "*"

AA_TO_CODONS = {}
for _codon, _aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def _translate_nt(nt: str) -> str:
    """Fast translation without Bio.Seq overhead."""
    out = []
    for i in range(0, len(nt) - 2, 3):
        out.append(CODON_TO_AA.get(nt[i:i + 3], "X"))
    return "".join(out)


# ---------------------------------------------------------------------------
# Worker state (set once per process via Pool initializer)
# ---------------------------------------------------------------------------
_WORKER: dict = {}


def _init_genome_worker(reference_fasta: str, peptide: str) -> None:
    _WORKER["ref"] = pysam.FastaFile(reference_fasta)
    _WORKER["peptide"] = peptide
    _WORKER["L"] = len(peptide) * 3
    _WORKER["chrom_len"] = {}


def _init_count_worker(cram_path: str, reference_fasta: str,
                       peptide: str, min_mapq: int,
                       require_full_span: bool) -> None:
    _WORKER["bam"] = pysam.AlignmentFile(
        cram_path, "rc", reference_filename=reference_fasta
    )
    _WORKER["peptide"] = peptide
    _WORKER["L"] = len(peptide) * 3
    _WORKER["min_mapq"] = min_mapq
    _WORKER["require_full_span"] = require_full_span


def _init_mapped_worker(cram_path: str, reference_fasta: str,
                        min_mapq: int) -> None:
    _WORKER["bam"] = pysam.AlignmentFile(
        cram_path, "rc", reference_filename=reference_fasta
    )
    _WORKER["min_mapq"] = min_mapq


# ---------------------------------------------------------------------------
# Locus finding
# ---------------------------------------------------------------------------
def find_peptide_loci_in_sequence(seq: str, peptide: str,
                                  chrom: str) -> List[Tuple[str, int, int, str]]:
    """
    Find all loci (both strands, all frames) encoding `peptide`.
    Uses translation + str.find — much faster than regex.
    """
    hits: List[Tuple[str, int, int, str]] = []
    L = len(peptide) * 3
    n = len(seq)

    for frame in range(3):
        sub = seq[frame: frame + ((n - frame) // 3) * 3]

        # --- forward strand ---
        prot = _translate_nt(sub)
        start = 0
        while True:
            j = prot.find(peptide, start)
            if j < 0:
                break
            gs = frame + j * 3
            hits.append((chrom, gs, gs + L, "+"))
            start = j + 1

        # --- reverse strand ---
        rc = str(Seq(sub).reverse_complement())
        prot_rc = _translate_nt(rc)
        sub_len = len(sub)
        start = 0
        while True:
            j = prot_rc.find(peptide, start)
            if j < 0:
                break
            gs = frame + sub_len - (j * 3 + L)
            hits.append((chrom, gs, gs + L, "-"))
            start = j + 1

    return hits


def _find_loci_in_window(task: Tuple[str, int, int, int]
                         ) -> List[Tuple[str, int, int, str]]:
    """
    Scan [start, end) of `chrom` with `pad` nt of overlap on each side so
    peptides straddling window boundaries are still found. Only hits
    *starting* inside [start, end) are returned, so windows tile without
    duplicates or gaps.
    """
    chrom, start, end, pad = task
    ref = _WORKER["ref"]
    peptide = _WORKER["peptide"]

    chrom_len = _WORKER["chrom_len"].get(chrom)
    if chrom_len is None:
        chrom_len = ref.get_reference_length(chrom)
        _WORKER["chrom_len"][chrom] = chrom_len

    s = max(0, start - pad)
    e = min(chrom_len, end + pad)
    seq = ref.fetch(chrom, s, e).upper()

    raw = find_peptide_loci_in_sequence(seq, peptide, chrom)

    out = []
    for (c, hs, he, strand) in raw:
        gs, ge = hs + s, he + s
        if start <= gs < end:
            out.append((c, gs, ge, strand))
    return out


def _plan_windows(reference_fasta: str,
                  window: int) -> List[Tuple[str, int, int, int]]:
    tasks = []
    with pysam.FastaFile(reference_fasta) as ref:
        for chrom in ref.references:
            length = ref.get_reference_length(chrom)
            for start in range(0, length, window):
                end = min(start + window, length)
                tasks.append((chrom, start, end, 0))
    return tasks


def find_all_loci_parallel(reference_fasta: str, peptide: str,
                           n_jobs: int = None,
                           window: int = 5_000_000
                           ) -> List[Tuple[str, int, int, str]]:
    """
    Scan the whole genome for peptide-encoding loci, load-balanced.
    `window` controls task granularity (5 Mb ≈ 600 tasks for human).
    """
    if n_jobs is None:
        n_jobs = min(mp.cpu_count(), 16)

    pad = 3 * len(peptide)
    tasks = [(c, s, e, pad)
             for (c, s, e, _) in _plan_windows(reference_fasta, window)]

    loci: List[Tuple[str, int, int, str]] = []
    with mp.Pool(
        processes=n_jobs,
        initializer=_init_genome_worker,
        initargs=(reference_fasta, peptide),
    ) as pool:
        for chunk in pool.imap_unordered(_find_loci_in_window, tasks,
                                         chunksize=1):
            loci.extend(chunk)

    return loci


# ---------------------------------------------------------------------------
# Locus clustering for safe parallel counting
# ---------------------------------------------------------------------------
@dataclass
class _Cluster:
    loci: List[Tuple[str, int, int, str]]
    start: int
    end: int
    chrom: str


def _cluster_loci(loci: Iterable[Tuple[str, int, int, str]],
                  gap: int) -> List[_Cluster]:
    """Group loci within `gap` nt on the same chromosome."""
    by_chrom: dict = defaultdict(list)
    for loc in loci:
        by_chrom[loc[0]].append(loc)

    clusters: List[_Cluster] = []
    for chrom, locs in by_chrom.items():
        locs.sort(key=lambda x: x[1])
        cur: List[Tuple[str, int, int, str]] = []
        cur_start = cur_end = None
        for loc in locs:
            if not cur:
                cur = [loc]
                cur_start, cur_end = loc[1], loc[2]
                continue
            if loc[1] - cur_end <= gap:
                cur.append(loc)
                cur_end = max(cur_end, loc[2])
            else:
                clusters.append(_Cluster(cur, cur_start, cur_end, chrom))
                cur = [loc]
                cur_start, cur_end = loc[1], loc[2]
        if cur:
            clusters.append(_Cluster(cur, cur_start, cur_end, chrom))
    return clusters


# ---------------------------------------------------------------------------
# Peptide-support read counting
# ---------------------------------------------------------------------------
def _count_reads_for_loci(loci: List[Tuple[str, int, int, str]]
                          ) -> Tuple[int, dict]:
    bam = _WORKER["bam"]
    peptide = _WORKER["peptide"]
    min_mapq = _WORKER["min_mapq"]
    require_full_span = _WORKER["require_full_span"]

    per_locus = defaultdict(int)
    seen_reads = set()

    for (chrom, start, end, strand) in loci:
        try:
            for read in bam.fetch(chrom, start, end):
                if (read.is_unmapped or read.is_secondary
                        or read.is_supplementary):
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                qseq = read.query_sequence
                if qseq is None:
                    continue

                if require_full_span:
                    if (read.reference_start > start
                            or read.reference_end < end):
                        continue

                ref_to_query = {}
                for q, r in read.get_aligned_pairs(matches_only=True):
                    if start <= r < end:
                        ref_to_query[r] = q

                nt = []
                ok = True
                for r in range(start, end):
                    q = ref_to_query.get(r)
                    if q is None:
                        ok = False
                        break
                    nt.append(qseq[q])
                if not ok:
                    continue
                nt = "".join(nt)

                if strand == "-":
                    nt = str(Seq(nt).reverse_complement())

                if _translate_nt(nt) != peptide:
                    continue

                key = (read.query_name, read.reference_start, strand)
                if key not in seen_reads:
                    seen_reads.add(key)
                    per_locus[(chrom, start, end, strand)] += 1
        except ValueError:
            continue

    return sum(per_locus.values()), dict(per_locus)


def count_reads_parallel(cram_path: str, reference_fasta: str,
                         loci: List[Tuple[str, int, int, str]],
                         peptide: str, min_mapq: int = 0,
                         require_full_span: bool = True,
                         n_jobs: int = None) -> dict:
    """
    Count supporting reads in parallel with greedy load balancing.
    """
    if n_jobs is None:
        n_jobs = min(mp.cpu_count(), 16)

    if not loci:
        return {"total": 0, "per_locus": {}}

    # Cluster so reads spanning multiple loci go to one worker.
    clusters = _cluster_loci(loci, gap=1000)

    # Greedy LPT: largest cluster first, always to lightest bucket.
    clusters.sort(key=lambda c: len(c.loci), reverse=True)
    buckets: List[List[_Cluster]] = [[] for _ in range(n_jobs)]
    weights = [0] * n_jobs
    for cl in clusters:
        i = weights.index(min(weights))
        buckets[i].append(cl)
        weights[i] += len(cl.loci)

    tasks = [
        [loc for cl in bucket for loc in cl.loci]
        for bucket in buckets
        if bucket
    ]

    with mp.Pool(
        processes=n_jobs,
        initializer=_init_count_worker,
        initargs=(cram_path, reference_fasta, peptide,
                  min_mapq, require_full_span),
    ) as pool:
        results = pool.map(_count_reads_for_loci, tasks)

    per_locus: dict = defaultdict(int)
    total = 0
    for t, pl in results:
        total += t
        for k, v in pl.items():
            per_locus[k] += v

    return {"total": total, "per_locus": dict(per_locus)}


# ---------------------------------------------------------------------------
# Per-window mapped-read counting
# ---------------------------------------------------------------------------
def _count_mapped_in_window(task: Tuple[str, int, int]
                            ) -> Tuple[Tuple[str, int, int], int]:
    """
    Count reads whose leftmost aligned position falls in [start, end).
    Filters mirror the peptide-support counter so the denominator matches.
    """
    chrom, start, end = task
    bam = _WORKER["bam"]
    min_mapq = _WORKER["min_mapq"]

    count = 0
    try:
        for read in bam.fetch(chrom, start, end):
            if (read.is_unmapped or read.is_secondary
                    or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            rs = read.reference_start
            if start <= rs < end:
                count += 1
    except ValueError:
        pass
    return (chrom, start, end), count


def _plan_mapped_windows(reference_fasta: str,
                         window: int) -> List[Tuple[str, int, int]]:
    tasks = []
    with pysam.FastaFile(reference_fasta) as ref:
        for chrom in ref.references:
            length = ref.get_reference_length(chrom)
            for start in range(0, length, window):
                end = min(start + window, length)
                tasks.append((chrom, start, end))
    return tasks


def count_mapped_per_window(cram_path: str, reference_fasta: str,
                            min_mapq: int = 0,
                            n_jobs: int = None,
                            window: int = 1_000_000) -> dict:
    """
    Return per-window mapped-read counts and the total.
    Each read is counted once, at its leftmost aligned position, so the
    per-window counts sum to the total primary mapped reads passing MAPQ.
    """
    if n_jobs is None:
        n_jobs = min(mp.cpu_count(), 16)

    tasks = _plan_mapped_windows(reference_fasta, window)

    per_window: dict = {}
    with mp.Pool(
        processes=n_jobs,
        initializer=_init_mapped_worker,
        initargs=(cram_path, reference_fasta, min_mapq),
    ) as pool:
        for key, cnt in pool.imap_unordered(_count_mapped_in_window, tasks,
                                            chunksize=4):
            per_window[key] = cnt

    total = sum(per_window.values())
    return {"per_window": per_window, "total": total}


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------
def quantify_peptide(cram_path: str, reference_fasta: str, peptide: str,
                     min_mapq: int = 0,
                     require_full_span: bool = True,
                     n_jobs: int = None,
                     window: int = 5_000_000,
                     mapped_window: int = 1_000_000) -> dict:
    loci = find_all_loci_parallel(
        reference_fasta, peptide, n_jobs=n_jobs, window=window
    )

    counts = count_reads_parallel(
        cram_path, reference_fasta, loci, peptide,
        min_mapq=min_mapq, require_full_span=require_full_span,
        n_jobs=n_jobs,
    )

    mapped = count_mapped_per_window(
        cram_path, reference_fasta, min_mapq=min_mapq,
        n_jobs=n_jobs, window=mapped_window,
    )
    counts["total_mapped"] = mapped["total"]
    counts["mapped_per_window"] = mapped["per_window"]

    counts["cpm"] = (counts["total"] / mapped["total"] * 1e6
                     if mapped["total"] else 0.0)
    counts["n_loci"] = len(loci)
    return counts


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    peptide = "RQFEPRIQY"
    cram_path = ("/home/dimitris/Workspace/postdoc/Darvac/wes_pipeline/"
                 "results/preprocessing/recalibrated/16H1303-2_LGD/"
                 "16H1303-2_LGD.recal.cram")
    reference_fasta = ("/home/dimitris/Workspace/reference/"
                       "Homo_sapiens_assembly38.fasta")

    result = quantify_peptide(
        cram_path, reference_fasta, peptide,
        n_jobs=16,
        window=5_000_000,
        mapped_window=1_000_000,
    )

    print("total peptide-supporting reads:", result["total"])
    print("total mapped reads (MAPQ-filtered, primary):", result["total_mapped"])
    print("CPM:", result["cpm"])
    print("n_loci:", result["n_loci"])
