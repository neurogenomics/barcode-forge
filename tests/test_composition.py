import pytest

from barcode_forge.composition import imbalance, positional_composition


def test_positional_composition_counts_fractions():
    comp = positional_composition(["AC", "AG", "TC", "GC"])
    assert comp[0] == {"A": 0.5, "C": 0.0, "G": 0.25, "T": 0.25}
    assert comp[1] == {"A": 0.0, "C": 0.75, "G": 0.25, "T": 0.0}


def test_positional_composition_requires_equal_lengths():
    with pytest.raises(ValueError):
        positional_composition(["AC", "A"])


def test_positional_composition_requires_sequences():
    with pytest.raises(ValueError):
        positional_composition([])


def test_imbalance_is_zero_for_perfectly_balanced_set():
    assert imbalance(["ACGT", "CGTA", "GTAC", "TACG"]) == pytest.approx(0.0)


def test_imbalance_is_worse_for_monotonous_set():
    assert imbalance(["AAAA", "AAAA"]) > imbalance(["AAAA", "CCCC"]) > 0.0
