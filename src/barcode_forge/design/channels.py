"""Imaging-channel balance of a barcode set, position by position.

Inline barcodes are read in the first sequencing cycles, so every cycle needs
signal in every imaging channel for cluster registration and base calling
(Illumina low-diversity and index-pooling guidance).

- 2-colour SBS (NextSeq, MiniSeq, NovaSeq): C = red, T = green, A = red + green,
  G = dark. A cycle that is mostly G is mostly dark.
- 4-colour SBS (MiSeq, HiSeq): A and C are read by the red laser, G and T by
  the green laser.
"""

from collections.abc import Mapping, Sequence

from barcode_forge.composition import positional_composition

# Bases lighting each channel.
CHANNELS = {
    "2colour": {"red": "AC", "green": "AT"},
    "4colour": {"red": "AC", "green": "GT"},
}
MIN_SIGNAL = 0.25  # least fraction of barcodes lighting each channel per cycle
MAX_DARK = 0.5  # most fraction of barcodes dark (G) per cycle on 2-colour
_EPS = 1e-9  # summed fractions such as 1/15 + 5/15 fall just below 0.4 as floats


def channel_signal(seqs: Sequence[str], chemistry: str) -> list[dict[str, float]]:
    """Fraction of barcodes lighting each channel at each position (plus "dark" on 2-colour)."""
    if chemistry not in CHANNELS:
        raise ValueError(f"chemistry must be one of {tuple(CHANNELS)}")
    return [position_signal(comp, chemistry) for comp in positional_composition(seqs)]


def position_signal(comp: Mapping[str, float], chemistry: str) -> dict[str, float]:
    """Channel fractions at one position from its base fractions (as in ``channel_signal``)."""
    pos = {ch: sum(comp[b] for b in bases) for ch, bases in CHANNELS[chemistry].items()}
    if chemistry == "2colour":
        pos["dark"] = comp["G"]
    return pos


def channel_violations(
    seqs: Sequence[str],
    chemistry: str | None,
    min_signal: float = MIN_SIGNAL,
    max_dark: float = MAX_DARK,
) -> list[str]:
    """Human-readable channel-balance problems; empty if balanced or chemistry is None."""
    if chemistry is None:
        return []
    problems = []
    for p, pos in enumerate(channel_signal(seqs, chemistry), 1):
        problems += [
            f"position {p}: {ch} {frac:.0%} < {min_signal:.0%}"
            for ch, frac in pos.items()
            if ch != "dark" and frac < min_signal - _EPS
        ]
        if pos.get("dark", 0.0) > max_dark + _EPS:
            problems.append(f"position {p}: dark {pos['dark']:.0%} > {max_dark:.0%}")
    return problems
