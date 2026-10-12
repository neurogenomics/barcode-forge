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
    # more than half a base outside the window (3 GC of 10 vs 3.6 / 2.4 bases)
    assert not RULE.check("GCGAAAAAAA", Config(gc_min=0.36, gc_max=0.5))
    assert not RULE.check("GCGAAAAAAA", Config(gc_min=0.0, gc_max=0.24))


def test_reexports_measurement():
    assert gc_fraction("GCAT") == 0.5


def test_metadata():
    assert RULE.name == "gc"
    assert 1 <= RULE.cost <= 2


@pytest.mark.parametrize(
    ("seq", "ok"),
    [
        ("ACGTACTA", True),  # 3/8 = 37.5%: 0.4 * 8 = 3.2 rounds to 3 GC bases
        ("ACGCACGA", True),  # 5/8 = 62.5%: 0.6 * 8 = 4.8 rounds to 5 GC bases
        ("ACATACTA", False),  # 2/8
        ("GCGCACGA", False),  # 6/8
    ],
)
def test_bounds_round_to_whole_bases_for_short_barcodes(seq, ok):
    # Without rounding, a 40-60% window admits only exactly 4/8 GC at 8 nt.
    assert RULE.check(seq, Config(length=8)) is ok


def test_rounding_never_narrows_the_window():
    # 12 nt: 0.4 * 12 = 4.8 -> 5 GC (41.7%); 4 GC (33%) stays out.
    assert not RULE.check("ACGTAATCTAAC", Config())  # 4 GC
    assert RULE.check("ACGTAGTCTAAC", Config())  # 5 GC
