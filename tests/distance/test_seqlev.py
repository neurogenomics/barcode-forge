import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.distance import levenshtein
from barcode_forge.distance.seqlev import distance
from barcode_forge.sequence import DNA

dna = st.text(alphabet=DNA, max_size=12)


@st.composite
def same_length(draw):
    size = draw(st.integers(0, 12))
    return draw(st.text(alphabet=DNA, min_size=size, max_size=size)), draw(
        st.text(alphabet=DNA, min_size=size, max_size=size)
    )


def by_definition(a: str, b: str) -> int:
    """Levenshtein after freely truncating either sequence at its 3' end.

    This is the minimum over the last row and last column of the Levenshtein DP
    matrix, written with an independent Levenshtein implementation.
    """
    return min(
        min(levenshtein.distance(a, b[:j]) for j in range(len(b) + 1)),
        min(levenshtein.distance(a[:i], b) for i in range(len(a) + 1)),
    )


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("", "", 0),
        ("ACGT", "ACGT", 0),
        # A substitution costs 1, as in Levenshtein.
        ("ACGT", "ACGA", 1),
        # Read lost its leading A: ACGT -> CGT, and the next downstream base (A)
        # fills position 4. One deletion; Levenshtein needs a second edit (2).
        ("ACGT", "CGTA", 1),
        # Extra A inserted after position 1: AACGT, read through a 4-nt window
        # gives AACG; the pushed-out T is free. Levenshtein = 2.
        ("ACGT", "AACG", 1),
        # Two leading bases deleted: ACGTAC -> GTAC, then two downstream bases
        # (TT). Two deletions; Levenshtein = 4 (also re-insert the tail).
        ("ACGTAC", "GTACTT", 2),
        # No shift helps when no bases are shared.
        ("AAAA", "TTTT", 4),
    ],
)
def test_known_values(a, b, expected):
    assert distance(a, b) == expected
    assert distance(b, a) == expected
    assert by_definition(a, b) == expected


@pytest.mark.parametrize(("a", "b"), [("ACGT", "CGTA"), ("ACGT", "AACG"), ("ACGTAC", "GTACTT")])
def test_smaller_than_levenshtein_for_shifts(a, b):
    assert distance(a, b) < levenshtein.distance(a, b)


def test_prefix_costs_nothing():
    # A sequence is 0 away from any extension of itself: the extra 3' bases are
    # treated as downstream sequence. So d == 0 iff equal only for equal lengths.
    assert distance("ACG", "ACGTT") == 0


def test_not_a_metric():
    # Triangle inequality fails via a shorter intermediate: AC and AG are each
    # extensions of A (distance 0), yet differ from each other by one substitution.
    assert distance("AC", "AG") > distance("AC", "A") + distance("A", "AG")


@given(dna, dna)
def test_matches_definition(a, b):
    assert distance(a, b) == by_definition(a, b)


@given(dna)
def test_identity(a):
    assert distance(a, a) == 0


@given(dna, dna)
def test_symmetry(a, b):
    assert distance(a, b) == distance(b, a)


@given(same_length())
def test_zero_iff_equal_for_equal_lengths(pair):
    a, b = pair
    assert (distance(a, b) == 0) == (a == b)


@given(same_length())
def test_at_most_levenshtein_for_equal_lengths(pair):
    a, b = pair
    assert distance(a, b) <= levenshtein.distance(a, b)


@given(dna, dna, st.integers(0, 12))
def test_score_cutoff(a, b, cutoff):
    true = distance(a, b)
    assert distance(a, b, score_cutoff=cutoff) == (true if true <= cutoff else cutoff + 1)


@pytest.mark.parametrize(("cutoff", "expected"), [(0, 1), (1, 2), (2, 2), (5, 2)])
def test_score_cutoff_known(cutoff, expected):
    assert distance("ACGTAC", "GTACTT", score_cutoff=cutoff) == expected
