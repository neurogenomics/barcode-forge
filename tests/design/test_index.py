import time
from collections.abc import Callable
from random import Random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from rapidfuzz.distance import Hamming, Levenshtein

from barcode_forge.design.candidates import generate
from barcode_forge.design.index import CloseIndex
from barcode_forge.distance import hamming, seqlev
from barcode_forge.sequence import DNA

Metric = Callable[..., int]


def _py_levenshtein(a: str, b: str, *, score_cutoff: int | None = None) -> int:
    """Pure-Python metric to exercise the generic path."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    d = prev[-1]
    return d if score_cutoff is None or d <= score_cutoff else score_cutoff + 1


def _naive(pool: list[str], seq: str, d: int, metric: Metric) -> bool:
    return any(metric(seq, b) < d for b in pool)


@st.composite
def pools(draw):
    """A pool and queries of one length; small alphabets make close pairs common."""
    length = draw(st.integers(5, 10))
    alphabet = draw(st.sampled_from(["AC", "ACG", DNA]))
    word = st.text(alphabet=alphabet, min_size=length, max_size=length)
    return draw(st.lists(word, max_size=12)), draw(st.lists(word, min_size=1, max_size=12))


def _check(pool: list[str], queries: list[str], d: int, metric: Metric) -> None:
    index = CloseIndex(metric, d)
    for s in pool:
        index.add(s)
    for q in queries + pool:
        assert index.has_close(q) == _naive(pool, q, d, metric), (q, pool, d)


@pytest.mark.parametrize(
    "metric", [Levenshtein.distance, Hamming.distance, hamming.distance, _py_levenshtein]
)
@settings(max_examples=100, deadline=None)
@given(data=pools(), d=st.integers(1, 5))
def test_matches_naive(metric: Metric, data, d: int) -> None:
    _check(*data, d, metric)


@settings(max_examples=1500, deadline=None)
@given(data=pools(), d=st.integers(1, 5))
def test_seqlev_matches_naive(data, d: int) -> None:
    _check(*data, d, seqlev.distance)


@settings(max_examples=300, deadline=None)
@given(
    pool=st.lists(st.text(alphabet="ACG", min_size=1, max_size=8), max_size=10),
    queries=st.lists(st.text(alphabet="ACG", min_size=1, max_size=8), min_size=1, max_size=10),
    d=st.integers(1, 5),
)
def test_seqlev_mixed_lengths_match_naive(pool: list[str], queries: list[str], d: int) -> None:
    _check(pool, queries, d, seqlev.distance)


def test_seqlev_known_cases() -> None:
    index = CloseIndex(seqlev.distance, 2)
    index.add("ACGTAC")
    assert index.has_close("ACGTAC")  # identical
    assert index.has_close("CGTACG")  # 5' deletion: free 3' shift, seqlev 1
    assert not index.has_close("TTTTTT")


def test_empty_index() -> None:
    assert not CloseIndex(seqlev.distance, 3).has_close("ACGTACGT")
    assert not CloseIndex(Levenshtein.distance, 3).has_close("ACGTACGT")


def test_seqlev_fast() -> None:
    rng = Random(0)
    index = CloseIndex(seqlev.distance, 4)
    for s in generate(10, 384, rng):
        index.add(s)
    queries = generate(10, 5000, rng)
    start = time.perf_counter()
    for q in queries:
        index.has_close(q)
    assert time.perf_counter() - start < 2.0


def test_hamming_wrapper_keeps_equal_length_check() -> None:
    index = CloseIndex(hamming.distance, 2)
    index.add("ACGT")
    with pytest.raises(ValueError):
        index.has_close("ACG")
