from collections import Counter
from random import Random

from hypothesis import given, settings
from hypothesis import strategies as st

from barcode_forge.composition import imbalance
from barcode_forge.design.balance import balance
from barcode_forge.design.candidates import generate
from barcode_forge.design.partition import partition


def _total(sets: list[list[str]]) -> float:
    return sum(imbalance(s) for s in sets)


@settings(max_examples=30, deadline=None)
@given(n_sets=st.integers(1, 4), set_size=st.integers(1, 12), seed=st.integers(0, 2**32))
def test_preserves_membership_sizes_and_never_worsens(
    n_sets: int, set_size: int, seed: int
) -> None:
    pool = generate(6, n_sets * set_size, Random(seed))
    sets = partition(pool, n_sets, set_size, Random(seed))
    snapshot = [list(s) for s in sets]
    out = balance(sets, Random(seed), rounds=200)
    assert sets == snapshot  # input not mutated
    assert [len(s) for s in out] == [len(s) for s in sets]
    assert Counter(b for s in out for b in s) == Counter(pool)
    assert _total(out) <= _total(sets) + 1e-12


def test_strictly_improves_obviously_unbalanced() -> None:
    # Each set is a single base: maximally unbalanced. Mixing them is ideal.
    sets = [["AAAA", "AAAA", "AAAA", "AAAA"], ["CCCC"] * 4, ["GGGG"] * 4, ["TTTT"] * 4]
    out = balance(sets, Random(0), rounds=2000)
    assert _total(out) < _total(sets)
    assert _total(out) == 0.0  # every set ends up with one of each base


def test_single_set_unchanged() -> None:
    sets = [["AAAA", "CCCC", "ACGT"]]
    assert balance(sets, Random(0)) == sets


def test_deterministic() -> None:
    pool = generate(8, 96, Random(3))
    sets = partition(pool, 4, 24, Random(3))
    assert balance(sets, Random(5), rounds=500) == balance(sets, Random(5), rounds=500)


def test_improves_random_partition() -> None:
    pool = generate(8, 4 * 24, Random(4))
    sets = partition(pool, 4, 24, Random(4))
    assert _total(balance(sets, Random(4))) < _total(sets)
