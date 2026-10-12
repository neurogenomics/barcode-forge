"""Per-position base composition of a set of equal-length sequences.

Inline barcodes are read in the first sequencing cycles, so each position
should carry a balanced mix of bases for cluster calling.
"""

from collections.abc import Sequence

from barcode_forge.sequence import DNA


def positional_composition(seqs: Sequence[str]) -> list[dict[str, float]]:
    """Fraction of each base at each position."""
    if not seqs:
        raise ValueError("need at least one sequence")
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise ValueError("sequences must all have the same length")
    n = len(seqs)
    return [{b: sum(s[i] == b for s in seqs) / n for b in DNA} for i in range(length)]


def imbalance(seqs: Sequence[str]) -> float:
    """Sum over positions of squared deviation from 25% per base (0 = perfect)."""
    return sum(
        (frac - 0.25) ** 2
        for position in positional_composition(seqs)
        for frac in position.values()
    )
