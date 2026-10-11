from itertools import product
from Bio.Data import CodonTable
from Bio.Seq import Seq
import re
import pysam
from collections import defaultdict


peptide = "RQFEPRIQY"

# Standard genetic code, codon -> amino acid
CODON_TABLE = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = {**CODON_TABLE.forward_table}
for stop in CODON_TABLE.stop_codons:
    CODON_TO_AA[stop] = '*'

# Reverse: amino acid -> list of codons
AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


def peptide_to_codon_sets(peptide: str):
    """Return list of codon lists, one per peptide position."""
    return [AA_TO_CODONS[aa] for aa in peptide]

def peptide_regex(peptide: str) -> str:
    """Build a regex matching any nucleotide encoding of the peptide."""
    parts = []
    for aa in peptide:
        codons = AA_TO_CODONS[aa]
        # Group codons: e.g. (TTT|TTC) for Phe
        parts.append("(" + "|".join(codons) + ")")
    return "".join(parts)

def find_peptide_loci_in_sequence(seq: str, peptide: str, chrom: str, 
                                  frame: int = 0):
    """Scan one frame (0, 1, or 2) of a sequence."""
    sub = seq[frame:]
    pattern = peptide_regex(peptide)
    L = len(peptide) * 3
    hits = []
    
    for m in re.finditer(f"(?=({pattern}))", sub):
        hits.append((chrom, frame + m.start(), frame + m.start() + L, "+"))
    
    rc = str(Seq(sub).reverse_complement())
    sub_len = len(sub)
    for m in re.finditer(f"(?=({pattern}))", rc):
        start_fwd = frame + sub_len - (m.start() + L)
        hits.append((chrom, start_fwd, start_fwd + L, "-"))
    
    return hits

def count_reads_for_peptide(bam_path: str,
                            reference_fasta : str,
                            loci: list, peptide: str, 
                             min_mapq: int = 0, 
                             require_full_span: bool = True) -> dict:
    """
    Count reads supporting a peptide at each candidate locus.
    
    loci: list of (chrom, start, end, strand)
    peptide: the peptide sequence (1-letter AA)
    
    Returns: dict with 'total' count and per-locus counts.
    """
    bam = pysam.AlignmentFile(
        cram_path, "rc",
        reference_filename=reference_fasta
        )
    per_locus = defaultdict(int)
    seen_reads = set()  # avoid double-counting a read across overlapping loci
    
    pep_len_nt = len(peptide) * 3

    for (chrom, start, end, strand) in loci:
        try:
            for read in bam.fetch(chrom, start, end):
                # ... (rest of your filtering logic is identical)
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.mapping_quality < min_mapq:
                    continue
                if read.query_sequence is None:
                    continue

                if require_full_span:
                    if read.reference_start > start:
                        continue
                    if read.reference_end < end:
                        continue 
                
                # Extract the read's aligned sequence in reference frame
                # get_aligned_pairs gives (query_pos, ref_pos) tuples
                pairs = read.get_aligned_pairs(matches_only=False)
                ref_to_query = {r: q for q, r in pairs 
                                if q is not None and r is not None}
                
                # Build the reference-frame nucleotide sequence of the read
                seq_chars = []
                for ref_pos in range(start, end):
                    qpos = ref_to_query.get(ref_pos)
                    if qpos is None:
                        seq_chars.append("N")
                    else:
                        seq_chars.append(read.query_sequence[qpos])
                nt = "".join(seq_chars)
                
                if "N" in nt:
                    continue
                
                # Account for strand
                if strand == "-":
                    nt = str(Seq(nt).reverse_complement())
                
                # Translate and check
                try:
                    translated = str(Seq(nt).translate())
                except Exception:
                    continue
                
                if translated == peptide:
                    read_key = (read.query_name, read.reference_start, strand)
                    if read_key not in seen_reads:
                        seen_reads.add(read_key)
                        per_locus[(chrom, start, end, strand)] += 1
        except ValueError:
            # chromosome not in BAM header
            continue
    
    bam.close()
    return {
        "total": sum(per_locus.values()),
        "per_locus": dict(per_locus),
    }

def normalize_count(count: int, total_mapped_reads: int, 
                    peptide_nt_len: int, 
                    total_transcriptome_nt: int = None) -> dict:
    """
    Return CPM and (if total_transcriptome_nt provided) a TPM-like value.
    """
    cpm = (count / total_mapped_reads) * 1e6
    
    result = {"raw": count, "cpm": cpm}
    
    if total_transcriptome_nt is not None:
        # TPM-like: rate normalized by peptide length, then scaled
        rate = count / (peptide_nt_len / 1000)  # reads per kb
        # This is only meaningful if you compute it for all features
        # and sum them; here we report the rate directly.
        result["rate_per_kb"] = rate
    
    return result

def quantify_peptide(bam_path: str, reference_fasta: str, peptide: str,
                     min_mapq: int = 0):
    """
    Full pipeline: find loci in the reference, then count reads.
    """
    # Step 1+2: find candidate loci
    loci = []
    with pysam.FastaFile(reference_fasta) as ref:
        for chrom in ref.references:
            seq = ref.fetch(chrom).upper()
            for frame in range(3):
                loci.extend(find_peptide_loci_in_sequence(
                    seq, peptide, chrom, frame))
    
    # Step 3: count
    counts = count_reads_for_peptide(bam_path, loci, peptide, min_mapq)
    
    # Normalization
    bam = pysam.AlignmentFile(bam_path, "rc")
    total_mapped = bam.mapped
    bam.close()
    
    counts["cpm"] = counts["total"] / total_mapped * 1e6
    counts["n_loci"] = len(loci)
    return counts

cram_path = "/home/dimitris/Workspace/postdoc/Darvac/wes_pipeline/results/preprocessing/recalibrated/16H1303-2_LGD/16H1303-2_LGD.recal.cram"
reference_fasta = "/home/dimitris/Workspace/reference/Homo_sapiens_assembly38.fasta"

quantify_peptide(cram_path , reference_fasta, peptide)
