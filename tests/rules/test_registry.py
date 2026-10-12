from barcode_forge.config import Config
from barcode_forge.rules import ALL_RULES, ruleset


def test_registry_contains_every_rule_once():
    names = [r.name for r in ALL_RULES]
    assert sorted(names) == sorted(
        [
            "gc",
            "homopolymer",
            "repeats",
            "palindrome",
            "motifs",
            "tm",
            "hairpin",
            "homodimer",
            "context",
            "dark_start",
        ]
    )


def test_ruleset_runs_string_rules_before_thermodynamic_ones():
    names = ruleset(Config()).names
    assert names.index("homopolymer") < names.index("hairpin")


def test_ruleset_can_exclude_rules():
    assert "tm" not in ruleset(Config(), exclude={"tm"}).names


def test_ruleset_rejects_a_bad_barcode_and_accepts_a_good_one():
    rs = ruleset(Config(tm_tolerance=None))
    assert not rs.passes("AAAAAAAAAA")
    assert rs.passes("ACGTCAGTGA") or rs.passes("ACAGTCTGAC")
