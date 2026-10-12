import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.distance import hamming
from barcode_forge.distance.levenshtein import distance
from barcode_forge.sequence import DNA

dna = st.text(alphabet=DNA, max_size=12)


@st.composite
def same_length(draw):
    size = draw(st.integers(0, 12))
    return draw(st.text(alphabet=DNA, min_size=size, max_size=size)), draw(
        st.text(alphabet=DNA, min_size=size, max_size=size)
    )


def reference(a: str, b: str) -> int:
    """Textbook Wagner-Fischer DP with unit costs."""
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("", "", 0),
        ("", "ACG", 3),  # three insertions
        ("ACGT", "ACGT", 0),
        ("ACGT", "ACGA", 1),  # one substitution
        ("ACGT", "ACT", 1),  # delete G
        ("ACGT", "CGTA", 2),  # delete leading A, append A (Hamming would say 4)
        ("AAAA", "TTTT", 4),  # four substitutions; no indel helps
        ("GATTACA", "GTTACAA", 2),  # delete A at 2, insert A at end
    ],
)
def test_known_values(a, b, expected):
    assert distance(a, b) == expected
    assert reference(a, b) == expected


@given(dna, dna)
def test_matches_reference_dp(a, b):
    assert distance(a, b) == reference(a, b)


@given(dna)
def test_identity(a):
    assert distance(a, a) == 0


@given(dna, dna)
def test_symmetry(a, b):
    assert distance(a, b) == distance(b, a)


@given(dna, dna)
def test_zero_iff_equal(a, b):
    assert (distance(a, b) == 0) == (a == b)


@given(dna, dna, dna)
def test_triangle_inequality(a, b, c):
    assert distance(a, c) <= distance(a, b) + distance(b, c)


@given(same_length())
def test_at_most_hamming(pair):
    a, b = pair
    assert distance(a, b) <= hamming.distance(a, b)


@given(dna, dna, st.integers(0, 12))
def test_score_cutoff(a, b, cutoff):
    true = distance(a, b)
    assert distance(a, b, score_cutoff=cutoff) == (true if true <= cutoff else cutoff + 1)


def test_is_a_rapidfuzz_scorer_so_pool_search_runs_in_c():
    # design.pool uses rapidfuzz.process for scorers carrying this attribute.
    assert hasattr(distance, "_RF_Scorer")
