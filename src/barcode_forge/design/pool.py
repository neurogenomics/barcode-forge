"""Greedy selection of a pool in which every pair is at least ``d`` apart."""

from collections.abc import Sequence

from rapidfuzz import process

from barcode_forge.distance import Metric


class InfeasibleError(ValueError):
    """The requested number of barcodes cannot be reached at this distance."""


def _is_rapidfuzz(metric: Metric) -> bool:
    return hasattr(metric, "_RF_Scorer")


def _has_close(candidate: str, accepted: Sequence[str], d: int, metric: Metric) -> bool:
    """True if some accepted barcode is closer than ``d`` to ``candidate``."""
    if not accepted:
        return False
    cutoff = d - 1
    if _is_rapidfuzz(metric):
        # One C-level scan instead of a Python loop over the pool.
        return (
            process.extractOne(candidate, accepted, scorer=metric, score_cutoff=cutoff) is not None
        )
    return any(metric(candidate, b, score_cutoff=cutoff) <= cutoff for b in accepted)


def greedy(candidates: Sequence[str], target: int, min_distance: int, metric: Metric) -> list[str]:
    """Accept candidates in order while they are ``>= min_distance`` from all accepted ones."""
    if target < 1:
        raise ValueError("target must be >= 1")
    if min_distance < 1:
        raise ValueError("min_distance must be >= 1")
    accepted: list[str] = []
    for candidate in candidates:
        if not _has_close(candidate, accepted, min_distance, metric):
            accepted.append(candidate)
            if len(accepted) == target:
                return accepted
    raise InfeasibleError(
        f"only {len(accepted)} of {target} barcodes reach min distance {min_distance} "
        f"from {len(candidates)} candidates; try a longer barcode length, a smaller "
        "min distance, fewer barcodes, more candidates or looser design rules"
    )
