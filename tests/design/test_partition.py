from random import Random

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from barcode_forge.design.partition import partition


@settings(max_examples=50, deadline=None)
@given(n_sets=st.integers(1, 6), set_size=st.integers(1, 20), seed=st.integers(0, 2**32))
def test_sizes_disjoint_and_union(n_sets: int, set_size: int, seed: int) -> None:
    pool = [f"S{i}" for i in range(n_sets * set_size)]
    original = list(pool)
    sets = partition(pool, n_sets, set_size, Random(seed))
    assert pool == original  # input not mutated
    assert len(sets) == n_sets
    assert all(len(s) == set_size for s in sets)
    flat = [b for s in sets for b in s]
    assert len(flat) == len(set(flat))
    assert set(flat) == set(pool)


def test_deterministic_and_shuffled() -> None:
    pool = [f"S{i:03}" for i in range(200)]
    a = partition(pool, 4, 50, Random(1))
    assert a == partition(pool, 4, 50, Random(1))
    assert [b for s in a for b in s] != pool


@pytest.mark.parametrize(("n", "k"), [(2, 3), (3, 3), (0, 5)])
def test_wrong_size_raises(n: int, k: int) -> None:
    with pytest.raises(ValueError):
        partition([f"S{i}" for i in range(5)], n, k, Random(0))
