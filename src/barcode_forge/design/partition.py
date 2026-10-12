"""Split a pool into equally sized sets."""

from collections.abc import Sequence
from random import Random


def partition(pool: Sequence[str], n_sets: int, set_size: int, rng: Random) -> list[list[str]]:
    """Shuffle a copy of ``pool`` and slice it into ``n_sets`` sets of ``set_size``."""
    if n_sets < 1 or set_size < 1 or len(pool) != n_sets * set_size:
        raise ValueError(f"cannot split {len(pool)} barcodes into {n_sets} sets of {set_size}")
    shuffled = list(pool)
    rng.shuffle(shuffled)
    return [shuffled[i * set_size : (i + 1) * set_size] for i in range(n_sets)]
