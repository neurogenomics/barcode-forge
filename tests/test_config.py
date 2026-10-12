import dataclasses

import pytest

from barcode_forge.config import Config, from_layers, load_toml


def test_defaults_follow_inline_barcode_guidelines():
    cfg = Config()
    assert (cfg.gc_min, cfg.gc_max) == (0.4, 0.6)
    assert cfg.max_homopolymer == 2
    assert cfg.min_distance is None  # maximise
    # Inline barcodes: an indel shifts downstream bases into the barcode window,
    # so plain Levenshtein's correction guarantee does not hold (Buschmann 2013).
    assert cfg.metric == "seqlev"


def test_tm_rule_is_opt_in():
    # An embedded barcode never anneals alone (primers bind the constant flanks),
    # so the Tm of the free 8-12-mer is not a PCR criterion; GC evens amplification.
    assert Config().tm_tolerance is None


def test_config_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        Config().length = 5  # type: ignore[misc]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"n_sets": 0},
        {"set_size": 0},
        {"length": 3},
        {"gc_min": 0.7, "gc_max": 0.6},
        {"gc_max": 1.5},
        {"min_distance": 0},
        {"min_distance": 11, "length": 10},
        {"metric": "euclid"},
        {"max_homopolymer": 0},
        {"chemistry": "3colour"},
        {"max_dark_start": -1},
        {"tm_tolerance": -0.5},
        {"seed": -5},  # Random(-5) gives the same stream as Random(5)
        {"avoid_motifs": ("GANTTC",)},
        {"flank5": "ACGN"},
    ],
)
def test_invalid_values_raise(kwargs):
    with pytest.raises(ValueError):
        Config(**kwargs)


def test_motifs_and_flanks_are_normalised_to_uppercase():
    cfg = Config(avoid_motifs=("gaattc",), flank5="acg", flank3="tt")
    assert cfg.avoid_motifs == ("GAATTC",)
    assert (cfg.flank5, cfg.flank3) == ("ACG", "TT")


def test_total_is_sets_times_size():
    assert Config(n_sets=3, set_size=8).total == 24


def test_with_overrides_ignores_none():
    cfg = Config().with_overrides(length=12, min_distance=None)
    assert cfg.length == 12
    assert cfg.min_distance is None


def test_load_toml(tmp_path):
    path = tmp_path / "cfg.toml"
    path.write_text('n_sets = 2\nlength = 12\navoid_motifs = ["GAATTC"]\n')
    assert load_toml(path) == {"n_sets": 2, "length": 12, "avoid_motifs": ("GAATTC",)}


def test_load_toml_rejects_unknown_keys(tmp_path):
    path = tmp_path / "cfg.toml"
    path.write_text("colour = 'blue'\n")
    with pytest.raises(ValueError, match="colour"):
        load_toml(path)


@pytest.mark.parametrize(
    ("rule", "field"),
    [
        ("homopolymer", "max_homopolymer"),
        ("repeats", "max_dinucleotide_repeat"),
        ("palindrome", "max_palindrome"),
        ("hairpin", "min_hairpin_dg"),
        ("homodimer", "min_homodimer_dg"),
        ("tm", "tm_tolerance"),
        ("dark_start", "max_dark_start"),
    ],
)
def test_without_disables_rule(rule, field):
    assert getattr(Config().without(rule), field) is None


def test_without_unknown_rule_raises():
    with pytest.raises(ValueError, match="gc"):
        Config().without("gc")


def test_load_toml_accepts_disable_list(tmp_path):
    path = tmp_path / "cfg.toml"
    path.write_text('disable = ["tm", "hairpin"]\n')
    assert load_toml(path) == {"disable": ("tm", "hairpin")}


def test_sequencing_defaults_target_two_colour_sbs():
    # NextSeq/NovaSeq 2-colour SBS is the most common platform; G is dark there,
    # so at most one leading G by default.
    cfg = Config()
    assert cfg.chemistry == "2colour"
    assert cfg.max_dark_start == 1
    assert cfg.check_revcomp is False


@pytest.mark.parametrize("chemistry", ["2colour", "4colour", None])
def test_valid_chemistries(chemistry):
    assert Config(chemistry=chemistry).chemistry == chemistry


def test_max_dark_start_zero_forbids_any_leading_g():
    assert Config(max_dark_start=0).max_dark_start == 0


def test_check_revcomp_is_settable():
    assert Config(check_revcomp=True).check_revcomp is True


NO_PRIMER3 = {"min_hairpin_dg": None, "min_homodimer_dg": None, "tm_tolerance": None}


@pytest.mark.parametrize(
    "active",
    [
        {"min_hairpin_dg": -3.0, "min_homodimer_dg": -6.0},  # the defaults
        {"min_hairpin_dg": -3.0},
        {"min_homodimer_dg": -6.0},
        {"tm_tolerance": 2.0},
    ],
)
def test_long_barcodes_need_primer3_rules_off(active):
    # primer3 refuses thermodynamic alignment of sequences over 60 nt.
    with pytest.raises(ValueError, match="hairpin.*homodimer.*tm"):
        Config(length=70, **(NO_PRIMER3 | active))


def test_long_barcodes_allowed_without_primer3_rules():
    assert Config(length=70, **NO_PRIMER3).length == 70
    assert Config(length=60).length == 60  # primer3 limit is inclusive


def test_zero_tm_tolerance_and_seed_are_valid():
    cfg = Config(tm_tolerance=0.0, seed=0)
    assert (cfg.tm_tolerance, cfg.seed) == (0.0, 0)


def test_chemistry_none_string_means_no_checks():
    # TOML and the CLI have no null, so "none" switches channel checks off.
    assert Config(chemistry="none").chemistry is None


def test_from_layers_applies_later_layers_and_disables_in_one_validation():
    # A long barcode is only valid once the primer3 rules are off, so disabling
    # must happen in the same step as setting the length.
    cfg = from_layers(
        {"length": 9},
        {"length": 70, "set_size": None},
        disable=("hairpin", "homodimer"),
    )
    assert cfg.length == 70
    assert cfg.set_size == Config().set_size  # None in a layer means "not given"
    assert cfg.min_hairpin_dg is None and cfg.min_homodimer_dg is None


def test_from_layers_rejects_unknown_rule():
    with pytest.raises(ValueError, match="gc"):
        from_layers(disable=("gc",))
