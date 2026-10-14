"""Codon translation and read-mRNA reconstruction."""
from Bio.Data import CodonTable

_T = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = {**_T.forward_table, **{s: "*" for s in _T.stop_codons}}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def translate(nt: str) -> str:
    """Translate a nucleotide string, ignoring a partial trailing codon."""
    return "".join(CODON_TO_AA.get(nt[i:i + 3], "X")
                   for i in range(0, len(nt) - 2, 3))


def read_mrna(read) -> str | None:
    """Transcript-space nucleotide sequence of the read, in mRNA orientation."""
    qs = read.query_sequence
    if qs is None:
        return None
    pairs = [(q, r) for q, r in read.get_aligned_pairs()
             if q is not None and r is not None]
    if read.is_reverse:
        pairs = pairs[::-1]
        return "".join(qs[q].translate(_COMP) for q, r in pairs)
    return "".join(qs[q] for q, r in pairs)
