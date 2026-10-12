"""GC content window.

40-60% GC gives even PCR amplification across barcodes and a stable Tm.
Bounds are applied as whole numbers of G/C bases, rounding each bound to the
nearest count, so short barcodes are not over-restricted (at 8 nt a strict
40-60% window admits only 4 GC; rounding admits 3-5).
"""

import math

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule
from barcode_forge.sequence import gc_fraction

__all__ = ["RULE", "gc_count_bounds", "gc_fraction"]

_EPS = 1e-9  # so bounds like 0.1 + 0.2 still admit exactly 3/10


def gc_count_bounds(length: int, gc_min: float, gc_max: float) -> tuple[int, int]:
    """Allowed G/C base counts: the strict window widened to the nearest whole base."""
    low, high = gc_min * length, gc_max * length
    lo = min(math.ceil(low - _EPS), math.floor(low + 0.5))
    hi = max(math.floor(high + _EPS), math.floor(high + 0.5))
    return lo, hi


def check(seq: str, cfg: Config) -> bool:
    lo, hi = gc_count_bounds(len(seq), cfg.gc_min, cfg.gc_max)
    return lo <= seq.count("G") + seq.count("C") <= hi


RULE = Rule(name="gc", cost=1, check=check)
