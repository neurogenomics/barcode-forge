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
CFG = Config(n_sets=2, set_size=2, length=6)


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
