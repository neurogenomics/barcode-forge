from random import Random

import pytest
from rapidfuzz.distance import Hamming, Levenshtein

from barcode_forge.design import search
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


def test_reserved_and_revcomp_pass_through() -> None:
    cands = generate(7, 3000, Random(6))
    reserved = ["GGGGGGG"]
    d, pool = maximise(cands, 20, Levenshtein.distance, 7, reserved=reserved, revcomp=True)
    assert pool == greedy(cands, 20, d, Levenshtein.distance, reserved=reserved, revcomp=True)
    assert all(Levenshtein.distance(p, "GGGGGGG") >= d for p in pool)


def test_reserved_infeasible_even_at_d1() -> None:
    with pytest.raises(InfeasibleError):
        maximise(["AAAA", "CCCC"], 2, Hamming.distance, 4, reserved=["CCCC"])


@pytest.mark.parametrize(("length", "target", "seed"), [(6, 20, 0), (7, 30, 1), (8, 60, 3)])
def test_each_d_evaluated_once(monkeypatch, length: int, target: int, seed: int) -> None:
    calls: list[int] = []

    def counting(cands, target, d, metric, **kw):
        calls.append(d)
        return greedy(cands, target, d, metric, **kw)

    monkeypatch.setattr(search, "greedy", counting)
    cands = generate(length, 3000, Random(seed))
    maximise(cands, target, Levenshtein.distance, length)
    assert len(calls) == len(set(calls)), calls
