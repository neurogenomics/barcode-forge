from random import Random

import pytest
from rapidfuzz.distance import Hamming, Levenshtein

from barcode_forge.design.candidates import generate
from barcode_forge.design.pool import InfeasibleError, greedy
from barcode_forge.design.search import maximise


def _feasible(cands: list[str], target: int, d: int, metric) -> bool:
    try:
        greedy(cands, target, d, metric)
    except InfeasibleError:
        return False
    return True


@pytest.mark.parametrize("metric", [Levenshtein.distance, Hamming.distance])
@pytest.mark.parametrize(("length", "target", "seed"), [(6, 20, 0), (7, 30, 1), (5, 8, 2)])
def test_beats_any_feasible_fixed_d(metric, length: int, target: int, seed: int) -> None:
    cands = generate(length, 3000, Random(seed))
    d, pool = maximise(cands, target, metric, length)
    best = max(x for x in range(1, length + 1) if _feasible(cands, target, x, metric))
    assert d >= best
    assert pool == greedy(cands, target, d, metric)
    assert len(pool) == target


def test_returns_full_length_when_trivial() -> None:
    # Two barcodes that differ at every position satisfy d == length.
    d, pool = maximise(["AAAA", "CCCC", "AAAC"], 2, Hamming.distance, 4)
    assert (d, pool) == (4, ["AAAA", "CCCC"])


def test_infeasible_even_at_d1() -> None:
    with pytest.raises(InfeasibleError):
        maximise(["AAAA", "CCCC"], 3, Hamming.distance, 4)


def test_deterministic() -> None:
    cands = generate(8, 2000, Random(9))
    assert maximise(cands, 40, Levenshtein.distance, 8) == maximise(
        cands, 40, Levenshtein.distance, 8
    )
