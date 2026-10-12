"""GC content window.

40-60% GC gives even PCR amplification across barcodes and a stable Tm.
"""

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule
from barcode_forge.sequence import gc_fraction

__all__ = ["RULE", "gc_fraction"]

_EPS = 1e-9  # so bounds like 0.1 + 0.2 still admit exactly 3/10


def check(seq: str, cfg: Config) -> bool:
    return cfg.gc_min - _EPS <= gc_fraction(seq) <= cfg.gc_max + _EPS


RULE = Rule(name="gc", cost=1, check=check)
