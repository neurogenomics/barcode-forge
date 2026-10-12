import dataclasses
from itertools import combinations

import pytest

from barcode_forge.config import Config
from barcode_forge.design.pipeline import Design, run
from barcode_forge.design.pool import InfeasibleError
from barcode_forge.distance import get_metric
from barcode_forge.rules import ruleset

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


def test_tm_target_is_resolved(design):
    assert design.cfg.tm_target is not None


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
