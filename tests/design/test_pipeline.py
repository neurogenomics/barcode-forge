import dataclasses
from itertools import combinations

import pytest

from barcode_forge.config import Config
from barcode_forge.design.channels import channel_violations
from barcode_forge.design.pipeline import ChannelBalanceError, Design, run
from barcode_forge.design.pool import InfeasibleError
from barcode_forge.distance import get_metric
from barcode_forge.rules import ruleset
from barcode_forge.sequence import reverse_complement

SMALL = Config(n_sets=3, set_size=6, length=8, seed=7, max_candidates=3000)


@pytest.fixture(scope="module")
def design() -> Design:
    return run(SMALL)


def test_shape(design):
    assert len(design.sets) == 3
    assert all(len(s) == 6 for s in design.sets)


def test_barcodes_unique_across_all_sets(design):
    flat = [b for s in design.sets for b in s]
    assert len(set(flat)) == len(flat) == SMALL.total


def test_every_pair_across_all_sets_meets_reported_distance(design):
    metric = get_metric(SMALL.metric)
    flat = [b for s in design.sets for b in s]
    assert min(metric(a, b) for a, b in combinations(flat, 2)) >= design.min_distance


def test_every_barcode_passes_every_rule_with_resolved_config(design):
    rs = ruleset(design.cfg)
    assert all(rs.passes(b) for s in design.sets for b in s)


def test_tm_target_is_resolved_when_tm_rule_enabled():
    assert run(SMALL.with_overrides(tm_tolerance=5.0)).cfg.tm_target is not None


def test_tm_target_left_unset_when_tm_rule_disabled(design):
    assert design.cfg.tm_target is None


def test_maximise_reaches_at_least_fixed_distance():
    fixed = run(SMALL.with_overrides(min_distance=2))
    assert fixed.min_distance == 2
    assert run(SMALL).min_distance >= 2


def test_same_seed_same_design(design):
    assert run(SMALL).sets == design.sets


def test_unset_seed_is_resolved_and_reproducible():
    first = run(dataclasses.replace(SMALL, seed=None))
    assert first.seed is not None
    assert run(SMALL.with_overrides(seed=first.seed)).sets == first.sets


def test_rejections_and_screened_counts_are_reported(design):
    assert design.candidates_screened >= SMALL.total
    assert sum(design.rejections.values()) > 0


def test_infeasible_design_raises():
    with pytest.raises(InfeasibleError):
        run(Config(n_sets=10, set_size=96, length=6, min_distance=3, seed=1))


def test_two_colour_designs_stay_away_from_poly_g(design):
    # Failed/empty clusters read as poly-G on 2-colour instruments.
    metric = get_metric(SMALL.metric)
    poly_g = "G" * SMALL.length
    assert all(metric(b, poly_g) >= design.min_distance for s in design.sets for b in s)


def test_every_set_is_channel_balanced(design):
    assert all(channel_violations(s, design.cfg.chemistry) == [] for s in design.sets)


def test_revcomp_separation_when_requested():
    cfg = SMALL.with_overrides(check_revcomp=True)
    result = run(cfg)
    metric = get_metric(cfg.metric)
    flat = [b for s in result.sets for b in s]
    assert all(metric(a, reverse_complement(b)) >= result.min_distance for a in flat for b in flat)


def test_unbalanced_set_fails_loudly():
    # A single barcode can never light both channels at every position.
    one = dataclasses.replace(SMALL, n_sets=1, set_size=1)
    with pytest.raises(ChannelBalanceError, match="--chemistry none"):
        run(one)
    assert run(dataclasses.replace(one, chemistry=None)).sets


@pytest.mark.parametrize("seed", range(5))
def test_small_sets_are_repaired_to_channel_balance(seed):
    # 5 barcodes per set often miss a channel minimum by chance; repair swaps in spares.
    cfg = Config(n_sets=2, set_size=5, length=8, seed=seed, max_candidates=2000)
    result = run(cfg)
    metric = get_metric(cfg.metric)
    flat = [b for s in result.sets for b in s]
    assert all(channel_violations(s, cfg.chemistry) == [] for s in result.sets)
    assert len(set(flat)) == len(flat) == cfg.total
    assert min(metric(a, b) for a, b in combinations(flat, 2)) >= result.min_distance
    assert all(metric(b, "G" * cfg.length) >= result.min_distance for b in flat)
    rs = ruleset(result.cfg)
    assert all(rs.passes(b) for b in flat)


def test_repair_keeps_revcomp_separation():
    cfg = Config(n_sets=2, set_size=5, length=8, seed=3, max_candidates=2000, check_revcomp=True)
    result = run(cfg)
    metric = get_metric(cfg.metric)
    flat = [b for s in result.sets for b in s]
    assert all(channel_violations(s, cfg.chemistry) == [] for s in result.sets)
    assert all(metric(a, reverse_complement(b)) >= result.min_distance for a in flat for b in flat)
    assert all(metric(reverse_complement(b), "G" * 8) >= result.min_distance for b in flat)
