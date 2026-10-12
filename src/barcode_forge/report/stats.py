"""Summary statistics of a finished design (pure data, no I/O)."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations
from statistics import fmean

from rapidfuzz import process

from barcode_forge.composition import imbalance, positional_composition
from barcode_forge.config import PRIMER3_MAX_LENGTH, Config
from barcode_forge.design.channels import channel_violations
from barcode_forge.distance import Metric, get_metric
from barcode_forge.rules.homopolymer import longest_homopolymer
from barcode_forge.rules.palindrome import longest_palindrome
from barcode_forge.rules.repeats import longest_dinucleotide_repeat
from barcode_forge.rules.thermo import hairpin_dg, homodimer_dg, tm
from barcode_forge.sequence import gc_fraction


@dataclass(frozen=True)
class Stats:
    min_distance: int
    metric: str
    n_sets: int
    set_size: int
    length: int
    within_min: int | None  # None when sets have one barcode
    within_mean: float | None
    cross_min: int | None  # None when there is one set
    cross_mean: float | None
    overall_min: int | None
    gc_min: float
    gc_max: float
    tm_min: float
    tm_max: float
    max_homopolymer: int
    max_dinucleotide_repeat: int
    max_palindrome: int
    worst_hairpin_dg: float | None  # None above primer3's 60 nt alignment limit
    worst_homodimer_dg: float | None
    composition: list[list[dict[str, float]]]  # per set, per position, per base
    imbalance: list[float]  # per set
    channel_violations: list[list[str]]  # per set, under cfg.chemistry
    rejections: dict[str, int]
    candidates_screened: int
    seed: int


def _row(query: str, choices: Sequence[str], metric: Metric) -> list[int]:
    """Exact distances from ``query`` to each choice, in order."""
    if hasattr(metric, "_RF_Scorer"):
        # rapidfuzz's C batch path; cdist would need numpy, which is not a dependency.
        hits = process.extract(query, choices, scorer=metric, limit=None)
        return [d for _, d, _ in sorted(hits, key=lambda hit: hit[2])]
    return [metric(query, c) for c in choices]


def _within(group: Sequence[str], metric: Metric) -> list[int]:
    return [d for i, a in enumerate(group[:-1]) for d in _row(a, group[i + 1 :], metric)]


def _cross(a: Sequence[str], b: Sequence[str], metric: Metric) -> list[int]:
    return [d for query in a for d in _row(query, b, metric)]


def _min_mean(ds: list[int]) -> tuple[int | None, float | None]:
    return (min(ds), fmean(ds)) if ds else (None, None)


def compute(
    sets: Sequence[Sequence[str]],
    cfg: Config,
    *,
    min_distance: int,
    rejections: Mapping[str, int],
    candidates_screened: int,
    seed: int,
) -> Stats:
    """Exact (no cutoff) distance and per-barcode property summary of ``sets``."""
    metric = get_metric(cfg.metric)
    within = [d for group in sets for d in _within(group, metric)]
    cross = [d for a, b in combinations(sets, 2) for d in _cross(a, b, metric)]
    within_min, within_mean = _min_mean(within)
    cross_min, cross_mean = _min_mean(cross)
    seqs = [s for group in sets for s in group]
    gcs = [gc_fraction(s) for s in seqs]
    tms = [tm(s, cfg) for s in seqs]
    aligns = len(seqs[0]) <= PRIMER3_MAX_LENGTH
    return Stats(
        min_distance=min_distance,
        metric=cfg.metric,
        n_sets=len(sets),
        set_size=len(sets[0]),
        length=len(seqs[0]),
        within_min=within_min,
        within_mean=within_mean,
        cross_min=cross_min,
        cross_mean=cross_mean,
        overall_min=min(within + cross, default=None),
        gc_min=min(gcs),
        gc_max=max(gcs),
        tm_min=min(tms),
        tm_max=max(tms),
        max_homopolymer=max(map(longest_homopolymer, seqs)),
        max_dinucleotide_repeat=max(map(longest_dinucleotide_repeat, seqs)),
        max_palindrome=max(map(longest_palindrome, seqs)),
        worst_hairpin_dg=min(hairpin_dg(s, cfg) for s in seqs) if aligns else None,
        worst_homodimer_dg=min(homodimer_dg(s, cfg) for s in seqs) if aligns else None,
        composition=[positional_composition(g) for g in sets],
        imbalance=[imbalance(g) for g in sets],
        channel_violations=[channel_violations(g, cfg.chemistry) for g in sets],
        rejections=dict(sorted(rejections.items())),
        candidates_screened=candidates_screened,
        seed=seed,
    )
