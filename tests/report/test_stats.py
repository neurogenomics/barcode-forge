from collections import Counter

import pytest

from barcode_forge.composition import imbalance, positional_composition
from barcode_forge.config import Config
from barcode_forge.distance import levenshtein
from barcode_forge.report import stats
from barcode_forge.rules.thermo import hairpin_dg, homodimer_dg, tm

# Hand-checkable Levenshtein distances (equal lengths, so d >= 6 - LCS):
#   within set 1: AAAAAA-AAAACC 2;  within set 2: CCCCCC-GAATTC 5
#   cross: AAAAAA-CCCCCC 6, AAAAAA-GAATTC 4, AAAACC-CCCCCC 4, AAAACC-GAATTC 3
SETS = [["AAAAAA", "AAAACC"], ["CCCCCC", "GAATTC"]]
ALL = [s for group in SETS for s in group]
CFG = Config(n_sets=2, set_size=2, length=6, metric="levenshtein")  # hand-computed values


def _compute(sets=SETS, cfg=CFG, **kwargs) -> stats.Stats:
    args = {"min_distance": 2, "rejections": {}, "candidates_screened": 10, "seed": 7} | kwargs
    return stats.compute(sets, cfg, **args)


def test_distances_match_hand_computed_values():
    s = _compute()
    assert (s.within_min, s.within_mean) == (2, 3.5)
    assert (s.cross_min, s.cross_mean) == (3, 4.25)
    assert s.overall_min == 2


def test_echoes_design_parameters():
    s = _compute(rejections=Counter(gc=3, homopolymer=1), candidates_screened=99, seed=42)
    assert (s.min_distance, s.metric, s.n_sets, s.set_size, s.length) == (
        2,
        "levenshtein",
        2,
        2,
        6,
    )
    assert s.rejections == {"gc": 3, "homopolymer": 1}
    assert (s.candidates_screened, s.seed) == (99, 42)


def test_sequence_property_extremes():
    s = _compute()
    assert (s.gc_min, s.gc_max) == (0.0, 1.0)
    assert s.max_homopolymer == 6
    assert s.max_dinucleotide_repeat == 1
    assert s.max_palindrome == 6  # GAATTC (EcoRI)


def test_thermodynamics_extremes():
    s = _compute()
    tms = [tm(x, CFG) for x in ALL]
    assert (s.tm_min, s.tm_max) == (min(tms), max(tms))
    assert s.worst_hairpin_dg == min(hairpin_dg(x, CFG) for x in ALL)
    assert s.worst_homodimer_dg == min(homodimer_dg(x, CFG) for x in ALL)


def test_per_set_composition_and_imbalance():
    s = _compute()
    assert s.composition == [positional_composition(g) for g in SETS]
    assert s.imbalance == [imbalance(g) for g in SETS]


def test_cross_set_is_none_for_one_set():
    s = _compute(sets=[ALL], cfg=Config(n_sets=1, set_size=4, length=6))
    assert s.cross_min is None and s.cross_mean is None
    assert s.within_min == s.overall_min == 2


def test_set_size_one_has_no_within_pairs():
    s = _compute(sets=[["AAAAAA"], ["GAATTC"]], cfg=Config(n_sets=2, set_size=1, length=6))
    assert s.within_min is None and s.within_mean is None
    assert s.cross_min == s.overall_min == 4


def test_single_barcode_has_no_pairs():
    s = _compute(sets=[["GAATTC"]], cfg=Config(n_sets=1, set_size=1, length=6))
    assert s.overall_min is None


def test_fast_path_agrees_with_loop(monkeypatch):
    fast = _compute()
    # A lambda hides rapidfuzz's _RF_Scorer marker, forcing the pure-Python loop.
    monkeypatch.setattr(
        stats, "get_metric", lambda name: lambda a, b, **kw: levenshtein.distance(a, b, **kw)
    )
    assert _compute() == fast


@pytest.mark.parametrize("metric", ["hamming", "seqlev"])
def test_other_metrics_use_loop(metric):
    s = _compute(cfg=Config(n_sets=2, set_size=2, length=6, metric=metric))
    assert s.metric == metric
    assert s.overall_min is not None and s.overall_min <= 2


def test_channel_violations_two_colour_hand_checked():
    # 2colour: red = A+C, green = A+T, dark = G.
    # Set 1 (AAAAAA, AAAACC): red 100% everywhere; green 100% (pos 1-4), 50% (pos 5-6).
    # Set 2 (CCCCCC, GAATTC): pos 1 C/G red 50 green 0 dark 50 (<= 50 OK);
    #   pos 2-3 C/A red 100 green 50; pos 4-5 C/T red 50 green 50; pos 6 C/C green 0.
    assert _compute().channel_violations == [
        [],
        ["position 1: green 0% < 25%", "position 6: green 0% < 25%"],
    ]


def test_channel_violations_four_colour_hand_checked():
    # 4colour: red = A+C, green = G+T.
    # Set 1: only A and C, so green 0% at all six positions.
    # Set 2: pos 1 C/G and pos 4-5 C/T are 50/50; pos 2-3 C/A and pos 6 C/C have green 0%.
    s = _compute(cfg=Config(n_sets=2, set_size=2, length=6, chemistry="4colour"))
    assert s.channel_violations == [
        [f"position {p}: green 0% < 25%" for p in range(1, 7)],
        [f"position {p}: green 0% < 25%" for p in (2, 3, 6)],
    ]


def test_no_channel_checks_without_chemistry():
    s = _compute(cfg=Config(n_sets=2, set_size=2, length=6, chemistry=None))
    assert s.channel_violations == [[], []]


def test_long_barcodes_skip_primer3_structure_stats():
    # primer3 cannot align > 60 nt, so hairpin/homodimer dG are not reported.
    cfg = Config(
        n_sets=1,
        set_size=2,
        length=70,
        min_hairpin_dg=None,
        min_homodimer_dg=None,
        chemistry=None,
    )
    s = _compute(sets=[["ACGT" * 17 + "AC", "TGCA" * 17 + "TG"]], cfg=cfg)
    assert s.worst_hairpin_dg is None and s.worst_homodimer_dg is None
    assert s.tm_min > 0
