"""Peptide TSV parsing and short/long split."""
from dataclasses import dataclass
from typing import Dict, Iterable, Set

# type alias: peptide -> set of target ids
PatternMap = Dict[str, Set[str]]

# Minimum length for a flexible (substring) match.
FLEX_MIN_LEN = 8


@dataclass
class PatternSets:
    """Peptides split into exact (strict) and flexible (substring) groups."""
    exact: PatternMap
    flexible: PatternMap


def read_peptides_tsv(path: str) -> PatternMap:
    """Read TSV (col1=peptide, col2=id, with header). No filtering."""
    out: PatternMap = {}
    with open(path) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            peptide, id_ = fields[0].strip(), fields[1].strip()
            if not peptide or not id_:
                continue
            out.setdefault(peptide, set()).add(id_)
    return out


def read_id_list(path: str) -> Set[str]:
    """One id per line; blank lines and '#' comments ignored."""
    ids: Set[str] = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line)
    return ids

def split_patterns(pattern_map: PatternMap,
                   flexible_ids: Set[str] | None = None) -> PatternSets:
    """
    Split into exact vs flexible based on whether any of a peptide's ids is
    in `flexible_ids`. If `flexible_ids` is None/empty, everything is exact.
    """
    flexible_ids = flexible_ids or set()
    exact, flexible = {}, {}
    for peptide, ids in pattern_map.items():
        (flexible if ids & flexible_ids else exact)[peptide] = ids
    return PatternSets(exact=exact, flexible=flexible)
