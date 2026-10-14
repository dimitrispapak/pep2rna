"""Per-BAM scan: fetch reads at peptide loci, translate, match."""
import os
from typing import Dict, Iterable, List, Set, Tuple

import pysam

from .patterns import FLEX_MIN_LEN
from .references import Locus, PeptideInfo
from .translate import read_mrna, translate

Row = Tuple[str, str, str, int, int, str, int, int]


def _peptide_match(peptide: str, aa: str, flexible: bool) -> bool:
    if not flexible:
        return peptide in aa
    for k in range(len(peptide), FLEX_MIN_LEN - 1, -1):
        for i in range(len(peptide) - k + 1):
            if peptide[i:i + k] in aa:
                return True
    return False


def process_bam(
    bam_file: str,
    bam_dir: str,
    peptides: Dict[str, PeptideInfo],
    min_mapq: int,
    flexible_set: Set[str],
) -> Tuple[str, int, List[Row]]:
    """Scan one BAM against all peptide loci."""
    bam_path = os.path.join(bam_dir, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")
    local_hits, rows = 0, []

    try:
        for peptide, info in peptides.items():
            if not info.loci:
                continue
            flexible = peptide in flexible_set

            for chrom, start, end, strand in info.loci:
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
                    if mrna is None:
                        continue

                    n_reads += 1
                    for frame in range(3):
                        aa = translate(mrna[frame:])
                        if _peptide_match(peptide, aa, flexible):
                            n_hits += 1
                            break

                if n_hits:
                    local_hits += n_hits
                    rows.append((bam_file, peptide, chrom_name,
                                 start_i, end_i, strand, n_reads, n_hits))
    finally:
        bam.close()

    return bam_file, local_hits, rows
