"""Levenshtein (edit) distance with unit costs for substitution, insertion and deletion."""

from rapidfuzz.distance import Levenshtein


def distance(a: str, b: str, *, score_cutoff: int | None = None) -> int:
    """Edit distance; returns ``score_cutoff + 1`` if it exceeds ``score_cutoff``."""
    return Levenshtein.distance(a, b, score_cutoff=score_cutoff)
