"""Find the largest minimum distance at which the greedy pool reaches the target."""

from collections.abc import Sequence

from barcode_forge.design.pool import InfeasibleError, Metric, greedy


def _try(candidates: Sequence[str], target: int, d: int, metric: Metric) -> list[str] | None:
    try:
        return greedy(candidates, target, d, metric)
    except InfeasibleError:
        return None


def maximise(
    candidates: Sequence[str], target: int, metric: Metric, length: int
) -> tuple[int, list[str]]:
    """Return ``(d, pool)`` for the largest d in ``[1, length]`` that greedy reaches.

    Greedy feasibility is close to, but not provably, monotone in d. After the
    binary search we keep stepping up while d + 1 also succeeds, so the result
    is never beaten by its immediate neighbour; gaps further up are not probed.
    """
    best = greedy(candidates, target, 1, metric)  # raises InfeasibleError if hopeless
    lo, hi = 1, length  # lo is known feasible
    while lo < hi:
        mid = (lo + hi + 1) // 2
        pool = _try(candidates, target, mid, metric)
        if pool is None:
            hi = mid - 1
        else:
            lo, best = mid, pool
    while lo < length and (pool := _try(candidates, target, lo + 1, metric)) is not None:
        lo, best = lo + 1, pool
    return lo, best
