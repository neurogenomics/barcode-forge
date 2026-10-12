from itertools import combinations

import pytest
from rapidfuzz.distance import Levenshtein
from typer.testing import CliRunner

from barcode_forge import __version__
from barcode_forge.cli import app

runner = CliRunner()
SMALL = ["generate", "-n", "2", "-s", "5", "-l", "8", "--seed", "3", "--max-candidates", "2000"]


def _invoke(*args):
    return runner.invoke(app, [str(a) for a in args])


def _read_sets(outdir):
    return [p.read_text().splitlines() for p in sorted(outdir.glob("set_*.tsv"))]


def test_version():
    result = _invoke("--version")
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help():
    result = _invoke("--help")
    assert result.exit_code == 0
    assert "barcodes" in result.output


def test_generate_writes_one_column_sets_and_stats(tmp_path):
    result = _invoke(*SMALL, "-o", tmp_path)
    assert result.exit_code == 0, result.output
    sets = _read_sets(tmp_path)
    assert [len(s) for s in sets] == [5, 5]
    assert all(len(b) == 8 and set(b) <= set("ACGT") for s in sets for b in s)
    assert (tmp_path / "stats.txt").exists()


def test_generate_respects_min_distance_across_sets(tmp_path):
    result = _invoke(*SMALL, "-d", "3", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    flat = [b for s in _read_sets(tmp_path) for b in s]
    assert len(set(flat)) == len(flat)
    assert min(Levenshtein.distance(a, b) for a, b in combinations(flat, 2)) >= 3


def test_same_seed_gives_identical_files(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    _invoke(*SMALL, "-o", a)
    _invoke(*SMALL, "-o", b)
    for name in ("set_01.tsv", "set_02.tsv", "stats.txt"):
        assert (a / name).read_bytes() == (b / name).read_bytes()


def test_config_file_is_used_and_flags_override_it(tmp_path):
    cfg = tmp_path / "cfg.toml"
    cfg.write_text('length = 9\nset_size = 4\ndisable = ["tm"]\n')
    result = _invoke(*SMALL, "--config", cfg, "-l", "10", "-o", tmp_path / "out")
    assert result.exit_code == 0, result.output
    sets = _read_sets(tmp_path / "out")
    assert [len(s) for s in sets] == [5, 5]  # -s 5 beats set_size = 4
    assert all(len(b) == 10 for s in sets for b in s)  # -l 10 beats length = 9
    assert "tm_tolerance = None" in (tmp_path / "out" / "stats.txt").read_text()


def test_disable_flag(tmp_path):
    result = _invoke(*SMALL, "--disable", "tm", "--disable", "homodimer", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    stats = (tmp_path / "stats.txt").read_text()
    assert "tm_tolerance = None" in stats
    assert "min_homodimer_dg = None" in stats


def test_avoid_motif_is_absent_on_both_strands(tmp_path):
    result = _invoke(*SMALL, "--avoid", "ACGT", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    assert not any("ACGT" in b for s in _read_sets(tmp_path) for b in s)


@pytest.mark.parametrize(
    "args",
    [
        ["generate", "-n", "10", "-s", "96", "-l", "6", "-d", "3", "--seed", "1"],
        ["generate", "-l", "2"],
        ["generate", "--metric", "euclid"],
        ["generate", "--disable", "gc"],
    ],
)
def test_bad_input_exits_nonzero_with_message(tmp_path, args):
    result = _invoke(*args, "-o", tmp_path)
    assert result.exit_code != 0
    assert "Error" in result.output
    assert "Traceback" not in result.output


def test_long_barcodes_work_when_primer3_rules_are_disabled(tmp_path):
    args = ["generate", "-n", "1", "-s", "3", "-l", "70", "-d", "10", "--seed", "1"]
    disable = ["--disable", "hairpin", "--disable", "homodimer"]
    # Three 70-mers cannot light both channels at all 70 positions, and ~50 spares
    # cannot repair that; this test is about primer3 limits, not channel balance.
    disable += ["--chemistry", "none"]
    result = _invoke(*args, *disable, "--max-candidates", "200", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    assert all(len(b) == 70 for s in _read_sets(tmp_path) for b in s)


def test_chemistry_option(tmp_path):
    result = _invoke(*SMALL, "--chemistry", "none", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    assert "chemistry = None" in (tmp_path / "stats.txt").read_text()


def test_revcomp_option_is_recorded(tmp_path):
    result = _invoke(*SMALL, "--revcomp", "-o", tmp_path)
    assert result.exit_code == 0, result.output
    assert "check_revcomp = True" in (tmp_path / "stats.txt").read_text()
