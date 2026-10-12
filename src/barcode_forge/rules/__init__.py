"""Biological design rules; one rule per module.

``ALL_RULES`` is the registry; ``ruleset(cfg)`` builds the RuleSet the design uses.
"""

from collections.abc import Collection

from barcode_forge.config import Config
from barcode_forge.rules import context, gc, homopolymer, motifs, palindrome, repeats, thermo
from barcode_forge.rules.ruleset import Rule, RuleSet

ALL_RULES: tuple[Rule, ...] = (
    gc.RULE,
    homopolymer.RULE,
    repeats.RULE,
    palindrome.RULE,
    motifs.RULE,
    *thermo.RULES,
    context.CONTEXT_RULE,
)


def ruleset(cfg: Config, exclude: Collection[str] = ()) -> RuleSet:
    return RuleSet((r for r in ALL_RULES if r.name not in exclude), cfg)
