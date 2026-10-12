"""Rules evaluated on the barcode inside its construct (flank5 + barcode + flank3).

A barcode that is clean on its own can still create problems at the junctions
with the adapter/linker it is cloned between. This rule only blames the
barcode for problems that involve at least one barcode base:

(a) a homopolymer run spanning a junction longer than ``max_homopolymer``
    (long runs raise indel error rates in sequencing-by-synthesis and nanopore);
(b) an avoided motif (either strand) overlapping the barcode, e.g. a
    restriction site formed across a junction;
(c) a hairpin in the local window (barcode plus up to 20 nt of each flank)
    below ``min_hairpin_dg`` *and* more than ``HAIRPIN_MARGIN`` more stable
    than the same window with the barcode replaced by an inert spacer (the
    least-structured of poly-A/C/G/T, since any one homopolymer can pair with
    a flank rich in its complement).
    The baseline captures structure the flanks form by themselves or with each
    other (e.g. the complementary ends of Illumina TruSeq adapters), which no
    barcode choice can fix; the margin absorbs barcode-dependent loop and
    terminal-mismatch terms (~1 nearest-neighbour stack), so only barcodes
    that actually pair into a stem are rejected.

Bare-barcode checks belong to the other rule modules; with no flanks this
rule passes.
"""

import re
from functools import cache

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule
from barcode_forge.rules.thermo import hairpin_dg
from barcode_forge.sequence import DNA, reverse_complement

FLANK_WINDOW = 20  # nt of each flank included in the hairpin window
MAX_THERMO_LEN = 60  # primer3's thermodynamic alignment limit
HAIRPIN_MARGIN = 1.0  # kcal/mol the barcode may add to the flanks' own structure


def local_window(seq: str, flank5: str, flank3: str) -> tuple[str, int, int]:
    """Barcode plus nearby flank, capped at primer3's 60 nt limit.

    Returns the window and the barcode's [start, end) within it.
    """
    room = max(0, (MAX_THERMO_LEN - len(seq)) // 2)
    n = min(FLANK_WINDOW, room)
    left = flank5[-n:] if n else ""
    right = flank3[:n]
    return left + seq + right, len(left), len(left) + len(seq)


def _junction_runs(seq: str, flank5: str, flank3: str):
    """Yield lengths of homopolymer runs touching a barcode base and a flank base."""
    construct = flank5 + seq + flank3
    start, end = len(flank5), len(flank5) + len(seq)
    for match in re.finditer(r"A+|C+|G+|T+", construct):
        touches_barcode = match.start() < end and match.end() > start
        spans_junction = match.start() < start or match.end() > end
        if touches_barcode and spans_junction:
            yield match.end() - match.start()


def longest_junction_homopolymer(seq: str, flank5: str, flank3: str) -> int:
    """Longest homopolymer run spanning a barcode/flank junction (0 if none)."""
    return max(_junction_runs(seq, flank5, flank3), default=0)


def _motif_overlaps_barcode(seq: str, cfg: Config) -> bool:
    construct = cfg.flank5 + seq + cfg.flank3
    start, end = len(cfg.flank5), len(cfg.flank5) + len(seq)
    for motif in cfg.avoid_motifs:
        for site in {motif, reverse_complement(motif)}:
            for match in re.finditer(f"(?={site})", construct):
                if match.start() < end and match.start() + len(site) > start:
                    return True
    return False


@cache  # the baseline depends only on the flanks and barcode length
def baseline_hairpin_dg(length: int, cfg: Config) -> float:
    """Hairpin ΔG of the window with the least-structured homopolymer spacer in place."""
    return max(
        hairpin_dg(local_window(base * length, cfg.flank5, cfg.flank3)[0], cfg) for base in DNA
    )


def _barcode_creates_hairpin(seq: str, cfg: Config) -> bool:
    window, _, _ = local_window(seq, cfg.flank5, cfg.flank3)
    dg = hairpin_dg(window, cfg)
    baseline = baseline_hairpin_dg(len(seq), cfg)
    return dg < cfg.min_hairpin_dg and dg < baseline - HAIRPIN_MARGIN


def _check(seq: str, cfg: Config) -> bool:
    if not cfg.flank5 and not cfg.flank3:
        return True
    if (
        cfg.max_homopolymer is not None
        and longest_junction_homopolymer(seq, cfg.flank5, cfg.flank3) > cfg.max_homopolymer
    ):
        return False
    if _motif_overlaps_barcode(seq, cfg):
        return False
    return cfg.min_hairpin_dg is None or not _barcode_creates_hairpin(seq, cfg)


CONTEXT_RULE = Rule(name="context", cost=20, check=_check)
