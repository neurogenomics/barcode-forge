"""Closeness index: does any indexed barcode lie closer than ``d`` to a query?"""

from rapidfuzz import process
from rapidfuzz.distance import Hamming, Levenshtein

from barcode_forge.distance import Metric, hamming, seqlev


def _any_close(seq: str, items: list[str], cutoff: int, scorer: Metric, **kwargs: object) -> bool:
    """One C-level scan of ``items`` instead of a Python loop."""
    return (
        process.extractOne(seq, items, scorer=scorer, score_cutoff=cutoff, scorer_kwargs=kwargs)
        is not None
    )


class CloseIndex:
    """Barcodes added so far, queried for any one at distance ``< d`` from a sequence.

    Sequence-Levenshtein of equal-length a, b (length L) is the minimum over the
    last row and column of the Levenshtein DP matrix:
        seqlev(a, b) = min(min_j lev(a, b[:j]), min_j lev(a[:j], b)).
    Only one side is ever truncated: cells truncating both are interior to the
    matrix and not part of the distance (lev("A", "A") = 0 < seqlev("AC", "AG") = 1).
    Since lev(a, b[:j]) >= L - j, only j >= L - d + 1 can give a value below d,
    so a query costs about 2d C scans over stored prefix lists (j = L is plain lev).
    Mixed lengths fall back to the generic per-pair loop.
    """

    def __init__(self, metric: Metric, d: int) -> None:
        self.metric, self.d = metric, d
        self._items: list[str] = []
        self._length: int | None = None  # common length, or None if empty or mixed
        self._prefixes: dict[int, list[str]] = {}  # j -> [b[:j] for b in items]

    def add(self, seq: str) -> None:
        if not self._items:
            self._length = len(seq)
        elif len(seq) != self._length:
            self._length = None
        self._items.append(seq)
        if self.metric is seqlev.distance and self._length is not None:
            for j in range(max(0, self._length - self.d + 1), self._length):
                self._prefixes.setdefault(j, []).append(seq[:j])

    def has_close(self, seq: str) -> bool:
        items, cutoff, metric = self._items, self.d - 1, self.metric
        if not items:
            return False
        if hasattr(metric, "_RF_Scorer"):
            return _any_close(seq, items, cutoff, metric)
        if metric is hamming.distance:
            # pad=False keeps the wrapper's equal-length check (rapidfuzz pads by default).
            return _any_close(seq, items, cutoff, Hamming.distance, pad=False)
        if metric is seqlev.distance and self._length == len(seq):
            lev = Levenshtein.distance
            return _any_close(seq, items, cutoff, lev) or any(
                _any_close(seq, self._prefixes[j], cutoff, lev)
                or _any_close(seq[:j], items, cutoff, lev)
                for j in self._prefixes
            )
        return any(metric(seq, b, score_cutoff=cutoff) <= cutoff for b in items)
