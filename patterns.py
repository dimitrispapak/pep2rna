"""Peptide TSV parsing and the short/long split."""
from dataclasses import dataclass
from typing import Dict, Set


DEFAULT_SHORT_IDS = {
    "AltP43s", "AltP11s", "AltP59s", "AltP42s", "AltP13s",
    "AltP53s", "AltP19s", "AltP39s", "AltP23s", "AltP34s",
}

LONG_IDS = {
    "AltP26L", "AltP25L", "AltP43L", "AltP59L", "AltP42L", "AltP13L",
    "AltP53L", "AltP19L", "AltP39L", "AltP40L", "AltP23L", "AltP34L",
}

FLEX_MIN_LEN = 8

# type alias: peptide -> set of target ids
PatternMap = Dict[str, Set[str]]


@dataclass
class PatternSets:
    """Peptides split into exact (short) and flexible (long) groups."""
    exact: PatternMap
    flexible: PatternMap

    @property
    def all_peptides(self):
        return {**self.exact, **self.flexible}


def read_peptides_tsv(path: str, keep_ids: Set[str]) -> PatternMap:
    """Read TSV (col1=peptide, col2=id, with header); keep only kept ids."""
    out: PatternMap = {}
    with open(path) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            peptide, id_ = fields[0].strip(), fields[1].strip()
            if id_ not in keep_ids:
                continue
            out.setdefault(peptide, set()).add(id_)
    return out


def split_short_long(pattern_map: PatternMap,
                     long_ids: Set[str] = LONG_IDS) -> PatternSets:
    exact, flexible = {}, {}
    for peptide, ids in pattern_map.items():
        (flexible if ids & long_ids else exact)[peptide] = ids
    return PatternSets(exact=exact, flexible=flexible)
