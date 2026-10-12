"""Leading dark cycles on Illumina 2-colour SBS (NextSeq, MiniSeq, NovaSeq).

In 2-colour chemistry G gives no signal. Inline barcodes are read in the first
cycles, where the instrument registers clusters and estimates phasing, so two
or more dark cycles at read start hurt cluster detection; empty or failed
clusters also read as poly-G. Illumina's low-diversity and index-pooling
guidance therefore advises against starting reads with GG.
"""

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule


def leading_dark_bases(seq: str) -> int:
    """Number of leading Gs (dark cycles on 2-colour SBS)."""
    return len(seq) - len(seq.lstrip("G"))


def check(seq: str, cfg: Config) -> bool:
    if cfg.chemistry != "2colour" or cfg.max_dark_start is None:
        return True
    return leading_dark_bases(seq) <= cfg.max_dark_start


RULE = Rule(name="dark_start", cost=1, check=check)
