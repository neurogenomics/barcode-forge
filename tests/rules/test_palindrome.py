import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.config import Config
from barcode_forge.rules.palindrome import RULE, longest_palindrome
from barcode_forge.sequence import DNA, reverse_complement

dna = st.text(alphabet=DNA, min_size=1, max_size=40)


@pytest.mark.parametrize(
    ("seq", "expected"),
    [
        ("GAATTC", 6),  # EcoRI
        ("AATT", 4),
        ("AT", 2),
        ("AA", 0),
        ("A", 0),
        ("AAAA", 0),
        ("CCGAATTCGG", 10),
        ("TTGAATTCAA", 10),
        ("AAGCAA", 2),  # GC
        ("AACCCAAGCTTCC", 6),  # AAGCTT (HindIII) embedded
    ],
)
def test_longest_palindrome(seq, expected):
    assert longest_palindrome(seq) == expected


@pytest.mark.parametrize(
    ("seq", "ok"),
    [("ACGTCAGGTC", True), ("AAGCTTCAGC", False), ("CAGAATTCAG", False), ("ACGTACGTAC", False)],
)
def test_default_max_four(seq, ok):
    assert RULE.check(seq, Config()) is ok


def test_boundary():
    assert RULE.check("GAATTCAC", Config(max_palindrome=6))
    assert not RULE.check("GGAATTCC", Config(max_palindrome=6))


def test_disabled_when_none():
    assert RULE.check("GAATTC", Config(max_palindrome=None))


@given(dna)
def test_length_is_even(seq):
    assert longest_palindrome(seq) % 2 == 0


@given(dna)
def test_strand_symmetric(seq):
    assert longest_palindrome(reverse_complement(seq)) == longest_palindrome(seq)


@given(dna)
def test_self_reverse_complement_is_whole_length(seq):
    pal = seq + reverse_complement(seq)
    assert longest_palindrome(pal) == len(pal)


@given(dna)
def test_matches_brute_force(seq):
    expected = max(
        (
            j - i
            for i in range(len(seq))
            for j in range(i + 1, len(seq) + 1)
            if seq[i:j] == reverse_complement(seq[i:j])
        ),
        default=0,
    )
    assert longest_palindrome(seq) == expected


def test_metadata():
    assert RULE.name == "palindrome"
    assert 1 <= RULE.cost <= 2
