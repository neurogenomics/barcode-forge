import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.distance.hamming import distance
from barcode_forge.sequence import DNA


@st.composite
def same_length(draw, n=2):
    size = draw(st.integers(0, 12))
    return tuple(draw(st.text(alphabet=DNA, min_size=size, max_size=size)) for _ in range(n))


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("", "", 0),
        ("ACGT", "ACGT", 0),
        ("ACGT", "ACGA", 1),  # last base substituted
        ("ACGT", "CGTA", 4),  # a one-base shift mismatches every position
        ("AAAA", "TTTT", 4),
        ("GATTACA", "GACTATA", 2),  # T->C at 3, C->T at 5
    ],
)
def test_known_values(a, b, expected):
    assert distance(a, b) == expected


def test_unequal_lengths_raise():
    with pytest.raises(ValueError):
        distance("ACGT", "ACG")


@pytest.mark.parametrize(("cutoff", "expected"), [(0, 1), (1, 2), (2, 3), (3, 4), (4, 4), (9, 4)])
def test_score_cutoff_known(cutoff, expected):
    assert distance("AAAA", "TTTT", score_cutoff=cutoff) == expected


@given(st.text(alphabet=DNA, max_size=12))
def test_identity(a):
    assert distance(a, a) == 0


@given(same_length())
def test_symmetry(pair):
    a, b = pair
    assert distance(a, b) == distance(b, a)


@given(same_length())
def test_zero_iff_equal(pair):
    a, b = pair
    assert (distance(a, b) == 0) == (a == b)


@given(same_length(3))
def test_triangle_inequality(triple):
    a, b, c = triple
    assert distance(a, c) <= distance(a, b) + distance(b, c)


@given(same_length(), st.integers(0, 12))
def test_score_cutoff(pair, cutoff):
    a, b = pair
    true = distance(a, b)
    assert distance(a, b, score_cutoff=cutoff) == (true if true <= cutoff else cutoff + 1)
