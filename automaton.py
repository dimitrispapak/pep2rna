"""Aho-Corasick automaton builders (exact and flexible)."""
import ahocorasick

from .patterns import FLEX_MIN_LEN, PatternMap


def build_exact(pattern_to_ids: PatternMap) -> ahocorasick.Automaton | None:
    if not pattern_to_ids:
        return None
    A = ahocorasick.Automaton()
    for p in pattern_to_ids:
        A.add_word(p, p)
    A.make_automaton()
    return A


def build_flexible(pattern_to_ids: PatternMap,
                   min_len: int = FLEX_MIN_LEN) -> ahocorasick.Automaton | None:
    if not pattern_to_ids:
        return None
    A = ahocorasick.Automaton()
    for p in pattern_to_ids:
        for k in range(min_len, len(p) + 1):
            for i in range(len(p) - k + 1):
                A.add_word(p[i:i + k], p)
    A.make_automaton()
    return A
