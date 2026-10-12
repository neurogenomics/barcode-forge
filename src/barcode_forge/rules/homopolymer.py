"""Maximum homopolymer run.

Homopolymers cause polymerase slippage and are the main source of
sequencing indel errors (especially on Nanopore/Ion Torrent).
"""

from itertools import groupby

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule


def longest_homopolymer(seq: str) -> int:
    return max((sum(1 for _ in run) for _, run in groupby(seq)), default=0)


def check(seq: str, cfg: Config) -> bool:
    return cfg.max_homopolymer is None or longest_homopolymer(seq) <= cfg.max_homopolymer


RULE = Rule(name="homopolymer", cost=1, check=check)
