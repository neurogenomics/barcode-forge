"""End-to-end design: screen candidates, build one global pool, split into sets.

A single pool in which every pair is ``>= d`` apart guarantees that barcodes
never collide between sets, so sets can also be pooled in one sequencing run.
"""

import dataclasses
from collections import Counter
from dataclasses import dataclass
from random import Random, SystemRandom

from barcode_forge.config import Config
from barcode_forge.design import candidates, partition, pool, search
from barcode_forge.design.balance import balance
from barcode_forge.distance import get_metric
from barcode_forge.rules import ruleset, thermo


@dataclass(frozen=True)
class Design:
    sets: list[list[str]]
    min_distance: int
    cfg: Config  # resolved: seed and tm_target filled in
    rejections: Counter[str]
    candidates_screened: int

    @property
    def seed(self) -> int:
        assert self.cfg.seed is not None
        return self.cfg.seed


def _resolve_tm(cfg: Config, screened: list[str]) -> Config:
    """Centre the Tm window on the median Tm of candidates passing every other rule."""
    if cfg.tm_tolerance is None or cfg.tm_target is not None or not screened:
        return cfg
    return dataclasses.replace(cfg, tm_target=thermo.estimate_tm_target(screened, cfg))


def run(cfg: Config) -> Design:
    if cfg.seed is None:
        cfg = dataclasses.replace(cfg, seed=SystemRandom().randrange(2**32))
    rng = Random(cfg.seed)

    # Tm needs a target derived from the candidates, so it is applied second.
    pre = ruleset(cfg, exclude={"tm"})
    screened = candidates.generate(cfg.length, cfg.max_candidates, rng, accept=pre.passes)
    cfg = _resolve_tm(cfg, screened)
    tm_only = ruleset(cfg, exclude=set(pre.names))
    passed = [s for s in screened if tm_only.passes(s)]
    rejections = pre.rejections + tm_only.rejections

    metric = get_metric(cfg.metric)
    if cfg.min_distance is None:
        d, chosen = search.maximise(passed, cfg.total, metric, cfg.length)
    else:
        d, chosen = cfg.min_distance, pool.greedy(passed, cfg.total, cfg.min_distance, metric)

    sets = partition.partition(chosen, cfg.n_sets, cfg.set_size, rng)
    if cfg.balance:
        sets = balance(sets, rng)
    return Design(
        sets=sets,
        min_distance=d,
        cfg=cfg,
        rejections=rejections,
        candidates_screened=len(screened) + pre.rejections.total(),
    )
