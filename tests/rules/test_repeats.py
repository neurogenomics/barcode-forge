import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.config import Config
from barcode_forge.rules.repeats import RULE, longest_dinucleotide_repeat
from barcode_forge.sequence import DNA


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        ("ATATAT", 3),
        ("GCGC", 2),
        ("ACGT", 1),
        ("ACGTCA", 1),
        ("CATATAG", 2),  # ATATA: two full copies
        ("GGACACACACTT", 4),
        ("ATATGCGCGC", 3),  # the longer of two repeats wins
        ("AAAA", 0),  # homopolymer, not a dinucleotide repeat
        ("A", 0),
    ],
)
def test_longest_dinucleotide_repeat(seq, expected):
    assert longest_dinucleotide_repeat(seq) == expected


@pytest.mark.parametrize(
    ("seq", "ok"),
    [("ACGTACGTAC", True), ("ATATCGGCAC", True), ("ATATATCGGC", False), ("GCAGAGAGTC", False)],
)
def test_default_max_two(seq, ok):
    assert RULE.check(seq, Config()) is ok


def test_boundary():
    assert RULE.check("ATATATCG", Config(max_dinucleotide_repeat=3))
    assert not RULE.check("ATATATAT", Config(max_dinucleotide_repeat=3))


def test_disabled_when_none():
    assert RULE.check("ATATATATAT", Config(max_dinucleotide_repeat=None))


@given(st.sampled_from([a + b for a in DNA for b in DNA if a != b]), st.integers(1, 15))
def test_pure_repeat(di, k):
    assert longest_dinucleotide_repeat(di * k) == k


@given(st.text(alphabet=DNA, min_size=2, max_size=40))
def test_bounded_by_half_length(seq):
    assert longest_dinucleotide_repeat(seq) <= len(seq) // 2


def test_metadata():
    assert RULE.name == "repeats"
    assert 1 <= RULE.cost <= 2
