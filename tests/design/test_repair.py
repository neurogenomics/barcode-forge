from itertools import combinations
from random import Random

import pytest

from barcode_forge.config import Config
from barcode_forge.design import pool
from barcode_forge.design.channels import channel_violations
from barcode_forge.design.repair import deficit, repair
from barcode_forge.distance import get_metric
from barcode_forge.sequence import reverse_complement

TWO = Config(length=4, chemistry="2colour")
HAMMING = get_metric("hamming")


def test_deficit_is_zero_for_a_balanced_set():
    assert deficit(["ACGT", "CATG", "GTAC", "TGCA"], "2colour") == 0


def test_deficit_sums_channel_shortfalls_and_dark_excess():
    # Positions 1-2 are G,G: red 0 (short 0.25), green 0 (short 0.25), dark 1 (excess 0.5).
    # Positions 3-4 are A,C: red 1, green 0.5, dark 0. Total 2 * (0.25 + 0.25 + 0.5).
    assert deficit(["GGAA", "GGCC"], "2colour") == pytest.approx(2.0)


def test_deficit_is_zero_exactly_when_there_are_no_violations():
    # 4-colour reads G as green, so C,G at every position is balanced.
    assert deficit(["CCCC", "GGGG"], "4colour") == 0
    # 2-colour: green 0 at all 4 positions -> 4 * 0.25.
    assert deficit(["CCCC", "GGGG"], "2colour") == pytest.approx(1.0)
    assert deficit(["GGGG"], None) == 0


def test_threshold_is_inclusive_despite_float_rounding():
    # red = 1/15 + 5/15 is 0.4 exactly, though not as a float sum (cf. channel_violations).
    seqs = ["AAAA"] + ["CCCC"] * 5 + ["GGGG"] * 9
    expected = 4 * (0.4 - 1 / 15)  # only green falls short
    assert deficit(seqs, "2colour", min_signal=0.4, max_dark=0.6) == pytest.approx(expected)


def test_repair_fixes_an_unbalanced_set_with_a_helpful_spare():
    # C,G everywhere: no green on 2-colour. Replacing G-rich GGGG with TTTT fixes it.
    fixed = repair([["CCCC", "GGGG"]], ["TTTT"], TWO, HAMMING, 2, Random(0))
    assert fixed == [["CCCC", "TTTT"]]
    assert channel_violations(fixed[0], "2colour") == []


def test_no_helpful_spare_returns_sets_unchanged():
    # GGGC in place of GGGG still leaves green 0 at positions 1-3 and 4 (deficit 1 -> 1).
    sets = [["CCCC", "GGGG"]]
    assert repair(sets, ["GGGC"], TWO, HAMMING, 2, Random(0)) == sets


def test_no_chemistry_means_nothing_to_repair():
    sets = [["CCCC", "GGGG"]]
    cfg = Config(length=4, chemistry=None)
    assert repair(sets, ["TTTT"], cfg, HAMMING, 2, Random(0)) == sets


def test_spare_too_close_to_a_barcode_in_another_set_is_not_used():
    # TTTT is 1 from TTTA in set 2; set 2 is itself balanced (T and A light green and red).
    sets = [["CCCC", "GGGG"], ["TTTA", "ACAT"]]
    assert repair(sets, ["TTTT"], TWO, HAMMING, 2, Random(0)) == sets


def test_spare_may_be_close_to_the_barcode_it_replaces():
    # TGGG is 1 from GGGG, which leaves the design. Position 1 becomes C,T (fixed).
    fixed = repair([["CCCC", "GGGG"]], ["TGGG"], TWO, HAMMING, 2, Random(0))
    assert fixed == [["CCCC", "TGGG"]]


def test_spare_too_close_to_a_reserved_sequence_is_not_used():
    sets = [["CCCC", "GGGG"]]
    out = repair(sets, ["TTTG"], TWO, HAMMING, 2, Random(0), reserved=("TTTT",))
    assert out == sets


def test_spare_close_to_its_own_reverse_complement_is_not_used_with_revcomp():
    # ATAT is its own reverse complement.
    sets = [["CCCC", "GGGG"]]
    assert repair(sets, ["ATAT"], TWO, HAMMING, 2, Random(0)) == [["CCCC", "ATAT"]]
    assert repair(sets, ["ATAT"], TWO, HAMMING, 2, Random(0), revcomp=True) == sets


def _random_design(seed: int, metric_name: str, revcomp: bool):
    """Small unbalanced sets from a greedy pool, plus the rejected-or-unused spares."""
    rng = Random(seed)
    cands = list(dict.fromkeys("".join(rng.choice("ACGT") for _ in range(6)) for _ in range(600)))
    metric = get_metric(metric_name)
    reserved = ("GGGGGG",)
    chosen = pool.greedy(cands, 9, 3, metric, reserved=reserved, revcomp=revcomp)
    sets = [chosen[0:3], chosen[3:6], chosen[6:9]]
    spares = [c for c in cands if c not in chosen]
    return sets, spares, metric, reserved


def _assert_invariants(sets, d, metric, reserved, revcomp):
    flat = [b for s in sets for b in s]
    assert len(set(flat)) == len(flat)
    assert all(metric(a, b) >= d for a, b in combinations(flat, 2))
    assert all(metric(b, r) >= d for b in flat for r in reserved)
    if revcomp:
        assert all(metric(a, reverse_complement(b)) >= d for a in flat for b in flat)
        assert all(metric(reverse_complement(b), r) >= d for b in flat for r in reserved)


@pytest.mark.parametrize("metric_name", ["levenshtein", "seqlev"])
@pytest.mark.parametrize("revcomp", [False, True])
@pytest.mark.parametrize("seed", range(4))
def test_repair_keeps_every_global_invariant(seed, metric_name, revcomp):
    sets, spares, metric, reserved = _random_design(seed, metric_name, revcomp)
    cfg = Config(length=6, chemistry="2colour")
    before = [deficit(s, "2colour") for s in sets]
    assert sum(before) > 0  # the fixture really needs repair
    out = repair(sets, spares, cfg, metric, 3, Random(seed), reserved=reserved, revcomp=revcomp)
    _assert_invariants(out, 3, metric, reserved, revcomp)
    assert [len(s) for s in out] == [len(s) for s in sets]
    assert set(b for s in out for b in s) <= set(b for s in sets for b in s) | set(spares)
    after = [deficit(s, "2colour") for s in out]
    assert all(a <= b + 1e-9 for a, b in zip(after, before, strict=True))
    assert sum(after) < sum(before)


def test_repair_is_deterministic_for_a_seed():
    sets, spares, metric, reserved = _random_design(1, "seqlev", True)
    cfg = Config(length=6, chemistry="2colour")

    def go(seed):
        return repair(sets, spares, cfg, metric, 3, Random(seed), reserved=reserved, revcomp=True)

    assert go(5) == go(5)


def test_repair_does_not_mutate_its_input():
    sets = [["CCCC", "GGGG"]]
    repair(sets, ["TTTT"], TWO, HAMMING, 2, Random(0))
    assert sets == [["CCCC", "GGGG"]]
