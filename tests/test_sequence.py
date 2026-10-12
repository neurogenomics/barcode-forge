import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.sequence import DNA, gc_fraction, reverse_complement, validate

dna = st.text(alphabet=DNA, min_size=1, max_size=40)


def test_alphabet_is_acgt():
    assert DNA == "ACGT"


@pytest.mark.parametrize(
    ("seq", "expected"),
    [("ACGT", "ACGT"), ("AAAA", "TTTT"), ("GATTACA", "TGTAATC"), ("", "")],
)
def test_reverse_complement(seq, expected):
    assert reverse_complement(seq) == expected


@given(dna)
def test_reverse_complement_is_involution(seq):
    assert reverse_complement(reverse_complement(seq)) == seq


@pytest.mark.parametrize(
    ("seq", "expected"),
    [("GGCC", 1.0), ("AATT", 0.0), ("ACGT", 0.5), ("GAAAA", 0.2)],
)
def test_gc_fraction(seq, expected):
    assert gc_fraction(seq) == pytest.approx(expected)


def test_gc_fraction_empty_raises():
    with pytest.raises(ValueError):
        gc_fraction("")


def test_validate_uppercases():
    assert validate("acgt") == "ACGT"


@pytest.mark.parametrize("bad", ["ACGN", "ACGU", "AC GT", ""])
def test_validate_rejects_non_dna(bad):
    with pytest.raises(ValueError):
        validate(bad)
