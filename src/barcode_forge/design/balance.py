"""Swap barcodes between sets to even out per-position base composition.

Inline barcodes are read in the first sequencing cycles; a set dominated by
one base at a position gives poor cluster calling, so each set should be as
close to 25% per base per position as possible.

Hill climbing: a random swap between two sets is kept only if the summed
``composition.imbalance`` of the sets strictly decreases. The change is
computed incrementally from integer base counts, which is exact and O(length)
per swap. Per set, imbalance = sum((4c - n)**2) / (16 n**2) over positions
and bases, where c is the base count and n the set size.
"""

from collections.abc import Sequence
from random import Random

from barcode_forge.sequence import DNA

_INDEX = {b: i for i, b in enumerate(DNA)}


def _counts(seqs: Sequence[str]) -> list[list[int]]:
    counts = [[0] * len(DNA) for _ in range(len(seqs[0]))] if seqs else []
    for s in seqs:
        for p, base in enumerate(s):
            counts[p][_INDEX[base]] += 1
    return counts


def _delta(counts: list[list[int]], n: int, out: str, into: str) -> int:
    """Change in sum((4c - n)**2) when ``out`` is replaced by ``into`` in a set."""
    total = 0
    for p, (a, b) in enumerate(zip(out, into, strict=True)):
        if a != b:
            ca, cb = counts[p][_INDEX[a]], counts[p][_INDEX[b]]
            total += -8 * (4 * ca - n) + 16 + 8 * (4 * cb - n) + 16
    return total


def _apply(counts: list[list[int]], out: str, into: str) -> None:
    for p, (a, b) in enumerate(zip(out, into, strict=True)):
        counts[p][_INDEX[a]] -= 1
        counts[p][_INDEX[b]] += 1


def balance(sets: Sequence[Sequence[str]], rng: Random, rounds: int = 20_000) -> list[list[str]]:
    """Return a copy of ``sets`` with barcodes swapped between sets to reduce imbalance."""
    out = [list(s) for s in sets]
    if len(out) < 2:
        return out
    counts = [_counts(s) for s in out]
    sizes = [len(s) for s in out]
    for _ in range(rounds):
        i, j = rng.sample(range(len(out)), 2)
        if not sizes[i] or not sizes[j]:
            continue
        x, y = rng.randrange(sizes[i]), rng.randrange(sizes[j])
        bx, by = out[i][x], out[j][y]
        di = _delta(counts[i], sizes[i], bx, by)
        dj = _delta(counts[j], sizes[j], by, bx)
        # Strict decrease of di / n_i**2 + dj / n_j**2, in exact integers.
        if di * sizes[j] ** 2 + dj * sizes[i] ** 2 < 0:
            _apply(counts[i], bx, by)
            _apply(counts[j], by, bx)
            out[i][x], out[j][y] = by, bx
    return out
