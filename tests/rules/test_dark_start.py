import pytest
from hypothesis import given
from hypothesis import strategies as st

from barcode_forge.config import Config
from barcode_forge.rules.dark_start import RULE, leading_dark_bases
from barcode_forge.sequence import DNA


@pytest.mark.parametrize(
    ("seq", "expected"),
    [("ACGT", 0), ("GACT", 1), ("GGACTCATGC", 2), ("AGGCTCATGC", 0), ("GGGG", 4)],
)
def test_leading_dark_bases(seq, expected):
    assert leading_dark_bases(seq) == expected


def test_two_dark_cycles_fail_on_two_colour_default():
    assert not RULE.check("GGACTCATGC", Config())


def test_internal_gs_are_not_dark_start():
    assert RULE.check("AGGCTCATGC", Config())


def test_one_leading_g_allowed_by_default():
    assert RULE.check("GACTCATGCA", Config())


def test_max_zero_forbids_leading_g():
    assert not RULE.check("GACTCATGCA", Config(max_dark_start=0))


def test_not_applied_on_four_colour():
    # 4-colour: G is imaged in the green channel, so it is not dark.
    assert RULE.check("GGACTCATGC", Config(chemistry="4colour"))


def test_not_applied_without_chemistry():
    assert RULE.check("GGACTCATGC", Config(chemistry=None))


def test_disabled_via_without():
    assert RULE.check("GGACTCATGC", Config().without("dark_start"))


@given(st.text(alphabet=DNA, max_size=30))
def test_bounded_by_length(seq):
    n = leading_dark_bases(seq)
    assert 0 <= n <= len(seq)
    assert seq[:n] == "G" * n


def test_metadata():
    assert RULE.name == "dark_start"
    assert RULE.cost == 1
