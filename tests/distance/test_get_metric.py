import pytest

from barcode_forge.config import METRICS
from barcode_forge.distance import get_metric, hamming, levenshtein, seqlev


@pytest.mark.parametrize(
    ("name", "module"),
    [("hamming", hamming), ("levenshtein", levenshtein), ("seqlev", seqlev)],
)
def test_returns_module_distance(name, module):
    assert get_metric(name) is module.distance


def test_every_config_metric_is_available():
    for name in METRICS:
        assert get_metric(name)("ACGT", "ACGT") == 0


def test_metrics_differ_on_a_shift():
    assert [get_metric(n)("ACGT", "CGTA") for n in ("hamming", "levenshtein", "seqlev")] == [
        4,
        2,
        1,
    ]


def test_unknown_metric_raises():
    with pytest.raises(ValueError, match="unknown metric"):
        get_metric("jaccard")
