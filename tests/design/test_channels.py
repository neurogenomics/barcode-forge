import pytest

from barcode_forge.design.channels import channel_signal, channel_violations

BALANCED = ["ACGT", "CATG", "GTAC", "TGCA"]  # each position has one of each base


def test_two_colour_signal_per_position():
    # GGAA, GGCC: positions 1-2 are G,G -> red 0, green 0, dark 1.
    # position 3-4 are A,C -> red f(A)+f(C) = 1, green f(A)+f(T) = 0.5, dark 0.
    assert channel_signal(["GGAA", "GGCC"], "2colour") == [
        {"red": 0.0, "green": 0.0, "dark": 1.0},
        {"red": 0.0, "green": 0.0, "dark": 1.0},
        {"red": 1.0, "green": 0.5, "dark": 0.0},
        {"red": 1.0, "green": 0.5, "dark": 0.0},
    ]


def test_four_colour_signal_per_position():
    # CCCC, GGGG: every position is C,G -> red f(A)+f(C) = 0.5, green f(G)+f(T) = 0.5.
    assert channel_signal(["CCCC", "GGGG"], "4colour") == [{"red": 0.5, "green": 0.5}] * 4


def test_all_dark_start_fails_two_colour():
    # Positions 1-2: red 0%, green 0%, dark 100% (> 50%); positions 3-4 OK
    # (red 100%, green 50%, dark 0%).
    assert channel_violations(["GGAA", "GGCC"], "2colour") == [
        "position 1: red 0% < 25%",
        "position 1: green 0% < 25%",
        "position 1: dark 100% > 50%",
        "position 2: red 0% < 25%",
        "position 2: green 0% < 25%",
        "position 2: dark 100% > 50%",
    ]


@pytest.mark.parametrize("chemistry", ["2colour", "4colour"])
def test_balanced_set_passes(chemistry):
    # Each position has one A, C, G, T: 2colour red 50%, green 50%, dark 25%;
    # 4colour red 50%, green 50%.
    assert channel_violations(BALANCED, chemistry) == []


def test_c_and_g_only_fails_two_colour_but_passes_four_colour():
    # 2colour, every position C,G: red f(C) = 50%, green f(A)+f(T) = 0%, dark 50% (<= 50%).
    assert channel_violations(["CCCC", "GGGG"], "2colour") == [
        f"position {p}: green 0% < 25%" for p in range(1, 5)
    ]
    # 4colour: red f(C) = 50%, green f(G) = 50%.
    assert channel_violations(["CCCC", "GGGG"], "4colour") == []


def test_no_chemistry_means_no_checks():
    assert channel_violations(["GGGG", "GGGG"], None) == []


def test_thresholds_are_adjustable():
    # Balanced set: dark 25% at every position, so max_dark=0.2 flags all four.
    assert len(channel_violations(BALANCED, "2colour", max_dark=0.2)) == 4
    # red/green 50% each, so min_signal=0.6 flags both channels at all 4 positions.
    assert len(channel_violations(BALANCED, "4colour", min_signal=0.6)) == 8


def test_unknown_chemistry_raises():
    with pytest.raises(ValueError, match="chemistry"):
        channel_signal(BALANCED, "3colour")


def test_threshold_is_inclusive_despite_float_rounding():
    # 15 barcodes, every position: 1 A, 5 C, 9 G -> red 1/15 + 5/15 = 40% exactly
    # (as floats 0.3999...), green 1/15 = 7% < 40%, dark 60%.
    seqs = ["AAAA"] + ["CCCC"] * 5 + ["GGGG"] * 9
    problems = channel_violations(seqs, "2colour", min_signal=0.4, max_dark=0.6)
    assert problems == [f"position {p}: green 7% < 40%" for p in range(1, 5)]
