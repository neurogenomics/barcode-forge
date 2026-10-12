import dataclasses

from barcode_forge import __version__
from barcode_forge.config import Config
from barcode_forge.report import stats, writer

SETS = [["AAAAAA", "AAAACC"], ["CCCCCC", "GAATTC"]]
CFG = Config(n_sets=2, set_size=2, length=6)


def _stats(seed: int = 7) -> stats.Stats:
    return stats.compute(
        SETS, CFG, min_distance=2, rejections={"gc": 5}, candidates_screened=12, seed=seed
    )


def test_write_sets_names_and_format(tmp_path):
    outdir = tmp_path / "new" / "out"
    paths = writer.write_sets(outdir, SETS)
    assert [p.name for p in paths] == ["set_01.tsv", "set_02.tsv"]
    for path, group in zip(paths, SETS, strict=True):
        text = path.read_text()
        assert text.endswith("\n")
        lines = text.splitlines()
        assert lines == group  # no header, one barcode per line
        assert len(lines) == CFG.set_size
        assert all("\t" not in line for line in lines)  # one column


def test_write_sets_zero_pads_to_width(tmp_path):
    sets = [[f"ACGT{a}{b}"] for a in "ACGT" for b in "ACG"]  # 12 sets
    names = [p.name for p in writer.write_sets(tmp_path, sets)]
    assert names[0] == "set_01.tsv" and names[-1] == "set_12.tsv"
    assert sorted(names) == names


def test_write_sets_pads_to_three_digits_for_100_sets(tmp_path):
    sets = [["ACGTAC"]] * 100
    names = [p.name for p in writer.write_sets(tmp_path, sets)]
    assert names[0] == "set_001.tsv" and names[-1] == "set_100.tsv"


def test_format_stats_header_and_parameters():
    text = writer.format_stats(_stats(seed=1234), CFG)
    assert text.splitlines()[0].startswith(f"barcode-forge {__version__}")
    assert "Parameters" in text
    for field in dataclasses.fields(Config):
        assert f"{field.name} = " in text
    assert "seed = 1234" in text  # resolved seed, not Config's None


def test_format_stats_sections_and_values():
    text = writer.format_stats(_stats(), CFG)
    for section in ("Distances", "Sequence properties", "Thermodynamics", "Rule rejections"):
        assert section in text
    assert "gc = 5" in text
    assert "set 01" in text and "set 02" in text
    assert "GAATTC" not in text  # barcodes live in the set files only


def test_format_stats_is_deterministic():
    assert writer.format_stats(_stats(), CFG) == writer.format_stats(_stats(), CFG)


def test_format_stats_shows_none_for_single_set():
    cfg = Config(n_sets=1, set_size=4, length=6)
    s = stats.compute(
        [SETS[0] + SETS[1]], cfg, min_distance=2, rejections={}, candidates_screened=4, seed=0
    )
    text = writer.format_stats(s, cfg)
    assert "cross-set" in text and "n/a" in text


def test_write_stats(tmp_path):
    path = writer.write_stats(tmp_path / "stats.txt", _stats(), CFG)
    assert path.read_text() == writer.format_stats(_stats(), CFG)
