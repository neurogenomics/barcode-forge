from collections.abc import Callable
from itertools import combinations
from random import Random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from rapidfuzz.distance import Hamming, Levenshtein

from barcode_forge.design.candidates import generate
from barcode_forge.design.pool import InfeasibleError, greedy
from barcode_forge.distance import seqlev
from barcode_forge.sequence import reverse_complement

Metric = Callable[..., int]


def _py_levenshtein(a: str, b: str, *, score_cutoff: int | None = None) -> int:
    """Pure-Python metric to exercise the generic (non-rapidfuzz) path."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    d = prev[-1]
    if score_cutoff is not None and d > score_cutoff:
        return score_cutoff + 1
    return d


def _min_pairwise(seqs: list[str], metric: Metric) -> int:
    return min(metric(a, b) for a, b in combinations(seqs, 2))


METRICS = [Levenshtein.distance, Hamming.distance, _py_levenshtein]


@pytest.mark.parametrize("metric", METRICS)
@settings(max_examples=15, deadline=None)
@given(seed=st.integers(0, 2**32), d=st.integers(1, 4), target=st.integers(2, 25))
def test_output_pairwise_at_least_d(metric: Metric, seed: int, d: int, target: int) -> None:
    candidates = generate(7, 2000, Random(seed))
    pool = greedy(candidates, target, d, metric)
    assert len(pool) == target
    assert len(set(pool)) == target
    assert set(pool) <= set(candidates)
    assert _min_pairwise(pool, metric) >= d


@pytest.mark.parametrize("metric", METRICS)
def test_matches_brute_force_greedy(metric: Metric) -> None:
    candidates = generate(6, 800, Random(11))
    d, target = 3, 40
    expected: list[str] = []
    for c in candidates:
        if all(metric(c, b) >= d for b in expected):
            expected.append(c)
        if len(expected) == target:
            break
    assert greedy(candidates, target, d, metric) == expected


def test_metric_called_with_score_cutoff() -> None:
    seen: list[int | None] = []

    def spy(a: str, b: str, *, score_cutoff: int | None = None) -> int:
        seen.append(score_cutoff)
        return _py_levenshtein(a, b, score_cutoff=score_cutoff)

    greedy(generate(6, 200, Random(0)), 5, 3, spy)
    assert seen and all(c == 2 for c in seen)


def test_stops_at_target() -> None:
    assert len(greedy(generate(8, 1000, Random(0)), 10, 2, Levenshtein.distance)) == 10


def test_deterministic_given_candidate_order() -> None:
    cands = generate(8, 2000, Random(5))
    assert greedy(cands, 50, 3, Levenshtein.distance) == greedy(cands, 50, 3, Levenshtein.distance)


def test_d1_accepts_all_unique() -> None:
    cands = ["AAAA", "AAAC", "AAAG"]
    assert greedy(cands, 3, 1, Hamming.distance) == cands


@pytest.mark.parametrize("metric", METRICS)
def test_infeasible_more_than_space(metric: Metric) -> None:
    cands = generate(2, 16, Random(0))
    with pytest.raises(InfeasibleError, match=r"(?i)length|distance|candidates"):
        greedy(cands, 17, 1, metric)


def test_infeasible_singleton_bound_hamming() -> None:
    # Singleton bound: a length-4 code with Hamming distance 4 has at most
    # 4**(4-4+1) = 4 words, so 5 barcodes are impossible.
    cands = generate(4, 256, Random(0))
    assert len(greedy(cands, 4, 4, Hamming.distance)) <= 4
    with pytest.raises(InfeasibleError):
        greedy(cands, 5, 4, Hamming.distance)


def test_infeasible_is_value_error() -> None:
    assert issubclass(InfeasibleError, ValueError)


def test_message_reports_progress() -> None:
    with pytest.raises(InfeasibleError, match="2"):
        greedy(["AAAA", "AAAC", "CCCC"], 3, 4, Hamming.distance)


def test_invalid_arguments() -> None:
    with pytest.raises(ValueError):
        greedy(["AAAA"], 1, 0, Hamming.distance)
    with pytest.raises(ValueError):
        greedy(["AAAA"], 0, 1, Hamming.distance)


SEQLEV_METRICS = [*METRICS, seqlev.distance]


@pytest.mark.parametrize("metric", SEQLEV_METRICS)
def test_reserved_respected_and_not_returned(metric: Metric) -> None:
    reserved = ["GGGGGGG"]
    cands = ["GGGGGGG", "GGGGGGA", *generate(7, 2000, Random(3))]
    pool = greedy(cands, 30, 3, metric, reserved=reserved)
    assert len(pool) == 30
    assert "GGGGGGG" not in pool and "GGGGGGA" not in pool
    assert all(metric(p, r) >= 3 for p in pool for r in reserved)


def test_reserved_makes_infeasible() -> None:
    with pytest.raises(InfeasibleError, match="reserved"):
        greedy(["AAAA", "GGGG"], 2, 1, Hamming.distance, reserved=["GGGG"])


@pytest.mark.parametrize("metric", SEQLEV_METRICS)
@settings(max_examples=10, deadline=None)
@given(seed=st.integers(0, 2**32), d=st.integers(1, 4))
def test_revcomp_pairs_at_least_d(metric: Metric, seed: int, d: int) -> None:
    pool = greedy(generate(9, 3000, Random(seed)), 12, d, metric, revcomp=True)
    assert len(pool) == 12
    for a in pool:
        for b in pool:
            assert metric(a, reverse_complement(b)) >= d, (a, b)
        assert pool.count(a) == 1


def test_revcomp_rejects_rc_and_palindromes() -> None:
    # ACGT is its own reverse complement; AACC's is GGTT.
    cands = ["ACGT", "AACC", "GGTT", "AAAA"]
    assert greedy(cands, 2, 1, Hamming.distance, revcomp=True) == ["AACC", "AAAA"]


@pytest.mark.parametrize("metric", SEQLEV_METRICS)
def test_defaults_unchanged(metric: Metric) -> None:
    cands = generate(7, 1500, Random(4))
    assert greedy(cands, 25, 3, metric) == greedy(cands, 25, 3, metric, reserved=(), revcomp=False)
    expected: list[str] = []
    for c in cands:
        if all(metric(c, b) >= 3 for b in expected):
            expected.append(c)
        if len(expected) == 25:
            break
    assert greedy(cands, 25, 3, metric) == expected
