import pytest

from barcode_forge.config import Config
from barcode_forge.rules import context
from barcode_forge.rules.ruleset import Rule

BENIGN = "ACAGTCATGA"
check = context.CONTEXT_RULE.check


def test_rule_exported():
    assert isinstance(context.CONTEXT_RULE, Rule)
    assert context.CONTEXT_RULE.name == "context"
    assert context.CONTEXT_RULE.cost >= 10


def test_no_flanks_always_passes():
    # Rules on the bare barcode are handled elsewhere; context adds nothing.
    assert check("GCGCGAAAAGCGCG", Config(avoid_motifs=("GAATTC",)))


def test_benign_flanks_pass():
    cfg = Config(flank5="ACTGTCAGTC", flank3="TGACTCAGTG")
    assert check(BENIGN, cfg)


# --- (a) junction homopolymers ---------------------------------------------


def test_junction_homopolymer_5prime_fails():
    cfg = Config(flank5="CTGCAA", max_homopolymer=2)
    assert not check("ACGTCAGTCA", cfg)  # ...AA|A... -> AAA


def test_junction_homopolymer_3prime_fails():
    cfg = Config(flank3="GGTCAG", max_homopolymer=2)
    assert not check("ACGTCAGTCG", cfg)  # ...G|GG... -> GGG


def test_junction_run_at_limit_passes():
    cfg = Config(flank5="CTGCTA", max_homopolymer=2)
    assert check("ACGTCAGTCA", cfg)  # ...A|A... -> AA


def test_homopolymer_inside_flank_only_is_ignored():
    cfg = Config(flank5="TTTTTTGC", max_homopolymer=2)
    assert check(BENIGN, cfg)


def test_junction_homopolymer_disabled_when_none():
    cfg = Config(flank5="CTGCAA", max_homopolymer=None, min_hairpin_dg=None)
    assert check("ACGTCAGTCA", cfg)


def test_longest_junction_homopolymer():
    assert context.longest_junction_homopolymer("ACGTCAGTCA", "CTGCAA", "") == 3
    assert context.longest_junction_homopolymer("ACGTCAGTCA", "CTGCTG", "TCAG") == 0


# --- (b) motifs created across a junction ----------------------------------


def test_motif_created_across_5prime_junction_fails():
    cfg = Config(flank5="CTGAGAA", avoid_motifs=("GAATTC",))
    assert not check("TTCGCAGTCA", cfg)  # GAA|TTC


def test_motif_on_reverse_strand_across_junction_fails():
    # GGATCC is a palindrome; use a non-palindromic motif: BsaI GGTCTC (rc GAGACC)
    cfg = Config(flank3="ACCTGA", avoid_motifs=("GGTCTC",), max_homopolymer=None)
    assert not check("ACAGTCGAG", cfg)  # GAG|ACC = rc(GGTCTC)


def test_motif_entirely_in_flank_is_ignored():
    cfg = Config(flank5="GAATTCTG", avoid_motifs=("GAATTC",))
    assert check(BENIGN, cfg)


# --- (c) hairpins created by the barcode in context -------------------------


def test_barcode_completing_stem_with_flank_fails():
    # flank5 ends with GCGCG; barcode starts with a loop then GCGCG -> stem.
    cfg = Config(flank5="ACTGAGCGCG", max_homopolymer=None)
    assert not check("AAAACGCGCA", cfg)


def test_hairpin_check_disabled_when_none():
    cfg = Config(flank5="ACTGAGCGCG", max_homopolymer=None, min_hairpin_dg=None)
    assert check("AAAACGCGCA", cfg)


def test_hairpin_already_in_flank_is_not_blamed_on_barcode():
    # flank5 carries its own strong hairpin; a benign barcode must still pass.
    cfg = Config(flank5="GCGCGAAAAGCGCG")
    assert context.baseline_hairpin_dg(10, cfg) < cfg.min_hairpin_dg
    assert check("TCAGTCATGA", cfg)


# Illumina TruSeq ends: the read-1 adapter's 3' end and the read-2 adapter's
# 5' end are reverse complements, so the construct folds whatever the barcode.
TRUSEQ5 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
TRUSEQ3 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"


def test_flank_to_flank_stem_is_not_blamed_on_barcode():
    cfg = Config(flank5=TRUSEQ5, flank3=TRUSEQ3)
    assert context.baseline_hairpin_dg(10, cfg) < -10
    assert check("ACAGTCATGA", cfg)


def test_barcode_extending_flank_to_flank_stem_fails():
    # GCC...GGC pairs with itself and continues the adapter stem by 3 GC pairs.
    cfg = Config(flank5=TRUSEQ5, flank3=TRUSEQ3)
    assert not check("GCCAGTTGGC", cfg)


@pytest.mark.parametrize("length", [10, 20, 30])
def test_window_stays_within_primer3_limit(length):
    flank = "ACTGTCAGTCACTGTCAGTCACTGTCAGTC"  # 30 nt
    seq = ("ACAGTCATGA" * 3)[:length]
    window, start, end = context.local_window(seq, flank, flank)
    assert len(window) <= 60
    assert window[start:end] == seq
