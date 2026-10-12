"""Hamming distance: substitutions only.

Suitable only for fixed-position indexes, where indels cannot shift the read.
"""

from rapidfuzz.distance import Hamming


def distance(a: str, b: str, *, score_cutoff: int | None = None) -> int:
    """Mismatch count; returns ``score_cutoff + 1`` if it exceeds ``score_cutoff``."""
    if len(a) != len(b):
        raise ValueError(f"Hamming distance needs equal lengths: {len(a)} != {len(b)}")
    return Hamming.distance(a, b, score_cutoff=score_cutoff)
