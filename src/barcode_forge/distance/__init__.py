"""Distance metrics between barcodes."""

from collections.abc import Callable

from barcode_forge.distance import hamming, levenshtein, seqlev

# distance(a, b, *, score_cutoff=None) -> int
Metric = Callable[..., int]

_METRICS: dict[str, Metric] = {
    "hamming": hamming.distance,
    "levenshtein": levenshtein.distance,
    "seqlev": seqlev.distance,
}


def get_metric(name: str) -> Metric:
    try:
        return _METRICS[name]
    except KeyError:
        raise ValueError(f"unknown metric {name!r}; choose from {sorted(_METRICS)}") from None
