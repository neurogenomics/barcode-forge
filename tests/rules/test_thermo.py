import primer3
import pytest

from barcode_forge.config import Config
from barcode_forge.rules import thermo
from barcode_forge.rules.ruleset import Rule

CFG = Config()
STRONG_HAIRPIN = "GCGCGAAAAGCGCG"  # GC-rich 5 bp stem + AAAA loop
SELF_DIMER = "GCGCGCGCGC"  # fully self-complementary, GC-rich
BENIGN = "ACAGTCATGA"


def _kwargs(cfg: Config) -> dict[str, float]:
    return dict(
        mv_conc=cfg.mv_conc, dv_conc=cfg.dv_conc, dntp_conc=cfg.dntp_conc, dna_conc=cfg.dna_conc
    )


# --- measurements -----------------------------------------------------------


def test_tm_gc_rich_above_at_rich():
    assert thermo.tm("GCGGCCGCGC", CFG) > thermo.tm("ATTAATATTA", CFG)


def test_tm_increases_with_length():
    assert thermo.tm("ACGTGCATGCAGTC", CFG) > thermo.tm("ACGTGCAT", CFG)


def test_tm_uses_config_salt():
    low_salt = Config(mv_conc=10.0, dv_conc=0.0)
    assert thermo.tm(BENIGN, low_salt) < thermo.tm(BENIGN, CFG)


def test_hairpin_dg_is_kcal_per_mol():
    raw = primer3.calc_hairpin(STRONG_HAIRPIN, **_kwargs(CFG)).dg
    assert thermo.hairpin_dg(STRONG_HAIRPIN, CFG) == pytest.approx(raw / 1000.0)


def test_homodimer_dg_is_kcal_per_mol():
    raw = primer3.calc_homodimer(SELF_DIMER, **_kwargs(CFG)).dg
    assert thermo.homodimer_dg(SELF_DIMER, CFG) == pytest.approx(raw / 1000.0)


def test_strong_hairpin_is_negative_and_plausible():
    dg = thermo.hairpin_dg(STRONG_HAIRPIN, CFG)
    assert -10.0 < dg < -3.0


def test_no_structure_gives_zero():
    assert thermo.hairpin_dg("AAAAAAAAAA", CFG) == 0.0
    assert thermo.homodimer_dg("AAAAAAAAAA", CFG) == 0.0


def test_self_complementary_dimer_more_stable_than_benign():
    assert thermo.homodimer_dg(SELF_DIMER, CFG) < thermo.homodimer_dg(BENIGN, CFG)


def test_estimate_tm_target_is_median():
    seqs = ["ATTAATATTA", BENIGN, "GCGGCCGCGC"]
    assert thermo.estimate_tm_target(seqs, CFG) == pytest.approx(thermo.tm(BENIGN, CFG))


def test_estimate_tm_target_even_count_is_mean_of_middle():
    seqs = ["ATTAATATTA", "GCGGCCGCGC"]
    expected = (thermo.tm(seqs[0], CFG) + thermo.tm(seqs[1], CFG)) / 2
    assert thermo.estimate_tm_target(seqs, CFG) == pytest.approx(expected)


def test_estimate_tm_target_empty_raises():
    with pytest.raises(ValueError):
        thermo.estimate_tm_target([], CFG)


# --- rules ------------------------------------------------------------------


def test_rules_exported_with_names_and_high_cost():
    assert thermo.RULES == (thermo.HAIRPIN_RULE, thermo.HOMODIMER_RULE, thermo.TM_RULE)
    assert [r.name for r in thermo.RULES] == ["hairpin", "homodimer", "tm"]
    for rule in thermo.RULES:
        assert isinstance(rule, Rule)
        assert rule.cost >= 10


@pytest.mark.parametrize(
    ("seq", "ok"),
    [(STRONG_HAIRPIN, False), (BENIGN, True), ("ACTGACTGAC", True)],
)
def test_hairpin_rule(seq, ok):
    assert thermo.HAIRPIN_RULE.check(seq, CFG) is ok


@pytest.mark.parametrize(
    ("seq", "ok"),
    [(SELF_DIMER, False), ("GGCCGGCCGGCC", False), (BENIGN, True)],
)
def test_homodimer_rule(seq, ok):
    assert thermo.HOMODIMER_RULE.check(seq, CFG) is ok


def test_hairpin_rule_disabled_when_none():
    assert thermo.HAIRPIN_RULE.check(STRONG_HAIRPIN, Config(min_hairpin_dg=None))


def test_homodimer_rule_disabled_when_none():
    assert thermo.HOMODIMER_RULE.check(SELF_DIMER, Config(min_homodimer_dg=None))


def test_tm_rule_window():
    target = thermo.tm(BENIGN, CFG)
    cfg = Config(tm_target=target, tm_tolerance=5.0)
    assert thermo.TM_RULE.check(BENIGN, cfg)
    assert not thermo.TM_RULE.check("GCGGCCGCGC", cfg)  # far hotter
    assert not thermo.TM_RULE.check("ATTAATATTA", cfg)  # far colder


def test_tm_rule_boundary_is_inclusive():
    target = thermo.tm(BENIGN, CFG) + 5.0
    assert thermo.TM_RULE.check(BENIGN, Config(tm_target=target, tm_tolerance=5.0))


def test_tm_rule_disabled_when_target_or_tolerance_none():
    assert thermo.TM_RULE.check("GCGGCCGCGC", Config(tm_target=None, tm_tolerance=5.0))
    assert thermo.TM_RULE.check("GCGGCCGCGC", Config(tm_target=0.0, tm_tolerance=None))
