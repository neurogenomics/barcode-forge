"""Melting temperature and secondary-structure rules (primer3 nearest-neighbour).

Conditions follow a standard PCR buffer: 50 mM monovalent cation, 1.5 mM Mg2+,
0.2 mM dNTP and 250 nM oligo (``Config.mv_conc`` etc.); Mg2+ and dNTP matter
because dNTPs chelate Mg2+ and free Mg2+ strongly stabilises duplexes.

Thresholds follow common oligo-design guidance (Premier Biosoft "PCR primer
design guidelines"; IDT OligoAnalyzer): hairpins should be weaker than about
-3 kcal/mol and self-dimers weaker than about -6 kcal/mol (for internal
structures) to avoid self-priming and loss of available template.
"""

import statistics
from collections.abc import Iterable

import primer3

from barcode_forge.config import Config
from barcode_forge.rules.ruleset import Rule


def _conditions(cfg: Config) -> dict[str, float]:
    return {
        "mv_conc": cfg.mv_conc,
        "dv_conc": cfg.dv_conc,
        "dntp_conc": cfg.dntp_conc,
        "dna_conc": cfg.dna_conc,
    }


def tm(seq: str, cfg: Config) -> float:
    """Melting temperature (°C) under the configured buffer."""
    return primer3.calc_tm(seq, **_conditions(cfg))


def _dg_kcal(result: primer3.thermoanalysis.ThermoResult) -> float:
    # primer3 reports ΔG in cal/mol.
    return result.dg / 1000.0 if result.structure_found else 0.0


def hairpin_dg(seq: str, cfg: Config) -> float:
    """Most stable hairpin ΔG (kcal/mol); 0.0 if none forms."""
    return _dg_kcal(primer3.calc_hairpin(seq, **_conditions(cfg)))


def homodimer_dg(seq: str, cfg: Config) -> float:
    """Most stable self-dimer ΔG (kcal/mol); 0.0 if none forms."""
    return _dg_kcal(primer3.calc_homodimer(seq, **_conditions(cfg)))


def estimate_tm_target(seqs: Iterable[str], cfg: Config) -> float:
    """Median Tm of ``seqs``; the pipeline centres the Tm window on it."""
    tms = [tm(s, cfg) for s in seqs]
    if not tms:
        raise ValueError("cannot estimate a Tm target from no sequences")
    return statistics.median(tms)


def _hairpin_ok(seq: str, cfg: Config) -> bool:
    return cfg.min_hairpin_dg is None or hairpin_dg(seq, cfg) >= cfg.min_hairpin_dg


def _homodimer_ok(seq: str, cfg: Config) -> bool:
    return cfg.min_homodimer_dg is None or homodimer_dg(seq, cfg) >= cfg.min_homodimer_dg


def _tm_ok(seq: str, cfg: Config) -> bool:
    # Similar Tm across barcodes gives even amplification efficiency.
    if cfg.tm_target is None or cfg.tm_tolerance is None:
        return True
    return abs(tm(seq, cfg) - cfg.tm_target) <= cfg.tm_tolerance


# Tm is cheapest of the three; all run after the string rules.
TM_RULE = Rule(name="tm", cost=10, check=_tm_ok)
HAIRPIN_RULE = Rule(name="hairpin", cost=15, check=_hairpin_ok)
HOMODIMER_RULE = Rule(name="homodimer", cost=20, check=_homodimer_ok)
RULES = (HAIRPIN_RULE, HOMODIMER_RULE, TM_RULE)
