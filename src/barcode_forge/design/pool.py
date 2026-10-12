"""Greedy selection of a pool in which every pair is at least ``d`` apart."""

from collections.abc import Sequence

from barcode_forge.design.index import CloseIndex
from barcode_forge.distance import Metric
from barcode_forge.sequence import reverse_complement


class InfeasibleError(ValueError):
    """The requested number of barcodes cannot be reached at this distance."""


def greedy(
    candidates: Sequence[str],
    target: int,
    min_distance: int,
    metric: Metric,
    *,
    reserved: Sequence[str] = (),
    revcomp: bool = False,
) -> list[str]:
    """Accept candidates in order while they are ``>= min_distance`` from all accepted ones.

    ``reserved`` sequences are kept ``>= min_distance`` away but never returned
    (e.g. poly-G, which failed clusters read as on 2-colour Illumina chemistry).
    With ``revcomp`` a barcode may be read on either strand (non-directional
    libraries), so each one must also be that far from the reverse complement of
    every accepted barcode, itself included.
    """
    if target < 1:
        raise ValueError("target must be >= 1")
    if min_distance < 1:
        raise ValueError("min_distance must be >= 1")
    d = min_distance
    fwd = CloseIndex(metric, d)  # reserved + accepted
    rev = CloseIndex(metric, d)  # rc(accepted)
    for r in reserved:
        fwd.add(r)
    accepted: list[str] = []
    for c in candidates:
        if fwd.has_close(c):
            continue
        if revcomp:
            rc = reverse_complement(c)
            # c vs rc(b), and b vs rc(c): equal for Levenshtein, not for seqlev.
            if rev.has_close(c) or fwd.has_close(rc) or metric(c, rc, score_cutoff=d - 1) < d:
                continue
            rev.add(rc)
        fwd.add(c)
        accepted.append(c)
        if len(accepted) == target:
            return accepted
    also = [f"{len(reserved)} reserved sequences"] * bool(reserved)
    also += ["reverse complements"] * revcomp
    avoiding = f" while avoiding {' and '.join(also)}" if also else ""
    raise InfeasibleError(
        f"only {len(accepted)} of {target} barcodes reach min distance {d} "
        f"from {len(candidates)} candidates{avoiding}; try a longer barcode length, "
        "a smaller min distance, fewer barcodes, more candidates or looser design rules"
    )
