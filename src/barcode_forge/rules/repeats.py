"""Dinucleotide (e.g. ATATAT) repeats.

Low-complexity repeats cause polymerase slippage and out-of-register mispriming.
"""

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule


def longest_dinucleotide_repeat(seq: str) -> int:
    """Most consecutive copies of a two-different-base unit; 0 if none occurs."""
    best = 0
    for i in range(len(seq) - 1):
        unit = seq[i : i + 2]
        if unit[0] == unit[1]:
            continue
        copies = 1
        while seq[i + 2 * copies : i + 2 * copies + 2] == unit:
            copies += 1
        best = max(best, copies)
    return best


def check(seq: str, cfg: Config) -> bool:
    limit = cfg.max_dinucleotide_repeat
    return limit is None or longest_dinucleotide_repeat(seq) <= limit


RULE = Rule(name="repeats", cost=2, check=check)
