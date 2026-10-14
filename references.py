"""Scan nuORFdb and MetamORF FASTA files for peptide loci."""
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Set, Tuple

from .patterns import PatternMap

# (chrom, start, end, strand)
Locus = Tuple[str, str, str, str]


@dataclass
class PeptideInfo:
    ids: Set[str] = field(default_factory=set)
    loci: Set[Locus] = field(default_factory=set)


def _iter_fasta(path: str):
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line.rstrip("\n")[1:], []
            else:
                chunks.append(line.rstrip("\n"))
    if header is not None:
        yield header, "".join(chunks)


class ReferenceScanner:
    """Base class: parse headers via `parse_header`, then Aho-Corasick match."""

    def __init__(self, automaton, pattern_to_ids: PatternMap):
        self.A = automaton
        self.patterns = pattern_to_ids

    def parse_header(self, header: str) -> Locus:
        raise NotImplementedError

    def scan(self, path: str) -> Dict[str, PeptideInfo]:
        results: Dict[str, PeptideInfo] = {}
        if self.A is None or not self.patterns:
            return results
        for header, seq in _iter_fasta(path):
            for _, p in self.A.iter(seq):
                info = results.setdefault(
                    p, PeptideInfo(ids=set(self.patterns[p])))
                info.loci.add(self.parse_header(header))
        return results


class MetamorfScanner(ReferenceScanner):
    _CHROM = re.compile(r"\bchr=(\S+)")
    _STRAND = re.compile(r"\bstrand=(\S+)")
    _START = re.compile(r"\bstart_pos=(\d+)")
    _END = re.compile(r"\bstop_pos=(\d+)")

    def parse_header(self, header: str) -> Locus:
        return (
            self._CHROM.search(header).group(1),
            self._START.search(header).group(1),
            self._END.search(header).group(1),
            self._STRAND.search(header).group(1),
        )


class NuorfScanner(ReferenceScanner):
    def parse_header(self, header: str) -> Locus:
        core = header.split("|")[0]
        chrom_part, span, strand = core.rsplit(":", 2)
        chrom = chrom_part.rsplit("_", 1)[-1]
        start, end = span.split("-")
        return chrom, start, end, strand


def merge_reference_results(*dicts) -> Dict[str, PeptideInfo]:
    """Merge several {peptide: PeptideInfo} dicts."""
    merged: Dict[str, PeptideInfo] = {}
    for d in dicts:
        for p, info in d.items():
            slot = merged.setdefault(p, PeptideInfo())
            slot.ids.update(info.ids)
            slot.loci.update(info.loci)
    return merged


def collapse_overlapping_loci(results: Dict[str, PeptideInfo]) -> None:
    """In-place: merge overlapping intervals per (chrom, strand)."""
    for info in results.values():
        by_key = defaultdict(list)
        for chrom, start, end, strand in info.loci:
            by_key[(chrom, strand)].append((int(start), int(end)))

        merged: Set[Locus] = set()
        for (chrom, strand), intervals in by_key.items():
            intervals.sort()
            s0, e0 = intervals[0]
            for s, e in intervals[1:]:
                if s <= e0:
                    e0 = max(e0, e)
                else:
                    merged.add((chrom, str(s0), str(e0), strand))
                    s0, e0 = s, e
            merged.add((chrom, str(s0), str(e0), strand))
        info.loci = merged
