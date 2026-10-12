import pytest

from barcode_forge.config import Config
from barcode_forge.rules.motifs import RULE, contains_motif

ECORI, BSAI = "GAATTC", "GGTCTC"


@pytest.mark.parametrize(
    ("seq", "motifs", "expected"),
    [
        ("ACGGAATTCA", (ECORI,), ECORI),
        ("ACGTACGTAC", (ECORI,), None),
        ("AAGGTCTCAA", (BSAI,), BSAI),  # forward strand
        ("AAGAGACCAA", (BSAI,), BSAI),  # reverse strand only
        ("AAGAGACCAA", (ECORI, BSAI), BSAI),
        ("GAATTCTCTC", (ECORI, BSAI), ECORI),  # first motif in the given order
        ("ACGTACGTAC", (), None),
    ],
)
def test_contains_motif(seq, motifs, expected):
    assert contains_motif(seq, motifs) == expected


def test_rule_fails_on_motif_either_strand():
    cfg = Config(avoid_motifs=(BSAI,))
    assert not RULE.check("AAGGTCTCAA", cfg)
    assert not RULE.check("AAGAGACCAA", cfg)
    assert RULE.check("ACGTACGTAC", cfg)


def test_no_motifs_always_passes():
    assert RULE.check("GAATTCGAATTC", Config())


def test_metadata():
    assert RULE.name == "motifs"
    assert 1 <= RULE.cost <= 2
