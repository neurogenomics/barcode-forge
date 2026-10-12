import dataclasses

import pytest

from barcode_forge.config import Config, load_toml


def test_defaults_follow_inline_barcode_guidelines():
    cfg = Config()
    assert (cfg.gc_min, cfg.gc_max) == (0.4, 0.6)
    assert cfg.max_homopolymer == 2
    assert cfg.min_distance is None  # maximise
    assert cfg.metric == "levenshtein"


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
