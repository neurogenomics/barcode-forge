"""Avoided motifs (restriction sites, adapter fragments), on both strands.

A barcode carrying a cut site or adapter substring is cleaved or misprimed.
"""

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule
from barcode_forge.sequence import reverse_complement


def contains_motif(seq: str, motifs: tuple[str, ...]) -> str | None:
    """Return the first motif present on either strand, else None."""
    rc = reverse_complement(seq)
    return next((m for m in motifs if m in seq or m in rc), None)


def check(seq: str, cfg: Config) -> bool:
    return contains_motif(seq, cfg.avoid_motifs) is None


RULE = Rule(name="motifs", cost=1, check=check)
