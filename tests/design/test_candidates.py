from random import Random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from barcode_forge.design.candidates import generate
from barcode_forge.sequence import DNA


def _ok(seqs: list[str], length: int) -> bool:
    return all(len(s) == length and set(s) <= set(DNA) for s in seqs)


@settings(max_examples=40, deadline=None)
@given(length=st.integers(1, 12), limit=st.integers(0, 300), seed=st.integers(0, 2**32))
def test_unique_correct_length_and_alphabet(length: int, limit: int, seed: int) -> None:
    out = generate(length, limit, Random(seed))
    assert len(out) == len(set(out))
    assert len(out) <= limit
    assert _ok(out, length)


@pytest.mark.parametrize("length", [3, 8, 14])
def test_same_seed_same_output(length: int) -> None:
    assert generate(length, 200, Random(7)) == generate(length, 200, Random(7))


@pytest.mark.parametrize("length", [6, 14])
def test_different_seeds_differ(length: int) -> None:
    assert generate(length, 50, Random(1)) != generate(length, 50, Random(2))


@pytest.mark.parametrize("length", [5, 14])
def test_accept_filter_honoured(length: int) -> None:
    def no_a(s: str) -> bool:
        return "A" not in s

    out = generate(length, 100, Random(0), accept=no_a)
    assert out
    assert all(no_a(s) for s in out)


def test_exhaustive_covers_whole_space_when_limit_large() -> None:
    out = generate(3, 10_000, Random(0))
    assert sorted(out) == sorted(a + b + c for a in DNA for b in DNA for c in DNA)


def test_exhaustive_order_is_shuffled() -> None:
    out = generate(4, 256, Random(0))
    assert out != sorted(out)


def test_random_mode_reaches_limit() -> None:
    out = generate(14, 1000, Random(3))
    assert len(out) == 1000


def test_random_mode_terminates_when_filter_rejects_everything() -> None:
    assert generate(14, 100, Random(0), accept=lambda s: False) == []


def test_exhaustive_mode_filter_can_reduce_below_limit() -> None:
    out = generate(2, 16, Random(0), accept=lambda s: s[0] == "A")
    assert sorted(out) == ["AA", "AC", "AG", "AT"]


def test_invalid_arguments() -> None:
    with pytest.raises(ValueError):
        generate(0, 10, Random(0))
    with pytest.raises(ValueError):
        generate(5, -1, Random(0))
