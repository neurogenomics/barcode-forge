"""Find the largest minimum distance at which the greedy pool reaches the target."""

from collections.abc import Sequence

from barcode_forge.design.pool import InfeasibleError, greedy
from barcode_forge.distance import Metric


def maximise(
    candidates: Sequence[str],
    target: int,
    metric: Metric,
    length: int,
    *,
    reserved: Sequence[str] = (),
    revcomp: bool = False,
) -> tuple[int, list[str]]:
    """Return ``(d, pool)`` for the largest d in ``[1, length]`` that greedy reaches.

    Greedy feasibility is close to, but not provably, monotone in d, so this is
    a binary search: the result's d + 1 is known to fail, but gaps further up
    are not probed.
    """
    # Raises an actionable InfeasibleError if even d = 1 (uniqueness) fails.
    pools: dict[int, list[str] | None] = {
        1: greedy(candidates, target, 1, metric, reserved=reserved, revcomp=revcomp)
    }

    def attempt(d: int) -> list[str] | None:
        if d not in pools:  # greedy is deterministic, so each d runs once
            try:
                pools[d] = greedy(candidates, target, d, metric, reserved=reserved, revcomp=revcomp)
            except InfeasibleError:
                pools[d] = None
        return pools[d]

    lo, hi = 1, length  # lo is known feasible
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if attempt(mid) is None:
            hi = mid - 1
        else:
            lo = mid
    best = pools[lo]
    assert best is not None
    return lo, best
