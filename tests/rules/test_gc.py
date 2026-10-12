import pytest

from barcode_forge.config import Config
from barcode_forge.rules.gc import RULE, gc_fraction


@pytest.mark.parametrize(
    ("seq", "ok"),
    [
        ("ACGTACGTAC", True),  # 50%
        ("AACGTACGTA", True),  # 40% lower bound
        ("GACGTACGTC", True),  # 60% upper bound
        ("AATGTACATA", False),  # 20%
        ("GCCGTGCGTC", False),  # 80%
        ("AAAAAAAAAA", False),
        ("GGGGGGGGGG", False),
    ],
)
def test_default_window(seq, ok):
    assert RULE.check(seq, Config()) is ok


def test_float_tolerance_at_bounds():
    # 0.1 + 0.2 == 0.30000000000000004 > 3/10; an exact bound must still pass
    assert RULE.check("GCGAAAAAAA", Config(gc_min=0.1 + 0.2, gc_max=0.5))
    assert RULE.check("GCGAAAAAAA", Config(gc_min=0.0, gc_max=0.7 - 0.4))


def test_just_outside_bounds_fails():
    assert not RULE.check("GCGAAAAAAA", Config(gc_min=0.31, gc_max=0.5))
    assert not RULE.check("GCGAAAAAAA", Config(gc_min=0.0, gc_max=0.29))


def test_reexports_measurement():
    assert gc_fraction("GCAT") == 0.5


def test_metadata():
    assert RULE.name == "gc"
    assert 1 <= RULE.cost <= 2
