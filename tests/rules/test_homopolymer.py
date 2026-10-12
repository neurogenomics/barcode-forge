import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.config import Config
from barcode_forge.rules.homopolymer import RULE, longest_homopolymer
from barcode_forge.sequence import DNA


@pytest.mark.parametrize(
    ("seq", "expected"),
    [("A", 1), ("ACGT", 1), ("AACGT", 2), ("ACGTTT", 3), ("GGGGACCC", 4), ("ACCA", 2)],
)
def test_longest_homopolymer(seq, expected):
    assert longest_homopolymer(seq) == expected


@pytest.mark.parametrize(
    ("seq", "ok"),
    [("ACGTACGTAC", True), ("AACCGGTTAC", True), ("ACGTTTACGA", False), ("AAAAAAAAAA", False)],
)
def test_default_max_two(seq, ok):
    assert RULE.check(seq, Config()) is ok


def test_boundary():
    assert RULE.check("ACGTTTAC", Config(max_homopolymer=3))
    assert not RULE.check("ACGTTTTA", Config(max_homopolymer=3))


def test_disabled_when_none():
    assert RULE.check("AAAAAAAAAA", Config(max_homopolymer=None))


@given(st.sampled_from(DNA), st.integers(1, 30))
def test_single_run(base, k):
    assert longest_homopolymer(base * k) == k


@given(st.text(alphabet=DNA, min_size=1, max_size=40))
def test_bounded_by_length(seq):
    assert 1 <= longest_homopolymer(seq) <= len(seq)


def test_metadata():
    assert RULE.name == "homopolymer"
    assert 1 <= RULE.cost <= 2
