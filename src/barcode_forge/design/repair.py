"""Repair imaging-channel balance by swapping barcodes for unused candidates.

Small sets often miss the per-cycle channel minimum by chance (e.g. one red
base too few at one position), which leaves early sequencing cycles without
signal in a channel (see ``channels``). Unused rule-passing candidates
("spares") can replace barcodes, provided the design's global invariants hold.

Hill climbing: take the set with the largest channel deficit and replace one
of its barcodes with a spare when that strictly lowers the set's deficit and
the spare is ``>= d`` from every other chosen barcode (in all sets), from the
reserved sequences and, with ``revcomp``, from reverse complements exactly as
in ``pool.greedy``. Each move lowers the total deficit, so the climb ends.
"""

from collections import Counter
from collections.abc import Sequence
from random import Random

from barcode_forge.config import Config
from barcode_forge.design.channels import MAX_DARK, MIN_SIGNAL, channel_signal, position_signal
from barcode_forge.design.index import CloseIndex
from barcode_forge.distance import Metric
from barcode_forge.sequence import DNA, reverse_complement

_EPS = 1e-9  # thresholds are inclusive despite float rounding (cf. channel_violations)


def _gap(x: float) -> float:
    return x if x > _EPS else 0.0


def _position_deficit(signal: dict[str, float], min_signal: float, max_dark: float) -> float:
    short = sum(_gap(min_signal - f) for ch, f in signal.items() if ch != "dark")
    return short + _gap(signal.get("dark", 0.0) - max_dark)


def deficit(
    seqs: Sequence[str],
    chemistry: str | None,
    min_signal: float = MIN_SIGNAL,
    max_dark: float = MAX_DARK,
) -> float:
    """Summed shortfall below ``min_signal`` per channel plus dark excess, over positions.

    Zero exactly when ``channel_violations`` reports nothing.
    """
    if chemistry is None or not seqs:
        return 0.0
    return sum(_position_deficit(s, min_signal, max_dark) for s in channel_signal(seqs, chemistry))


class _Set:
    """Per-position base counts of one set, for O(length) move evaluation."""

    def __init__(self, seqs: list[str], chemistry: str) -> None:
        self.seqs, self.chemistry = seqs, chemistry
        self.counts = [Counter(col) for col in zip(*seqs, strict=True)]
        self.pos = [self._deficit(c) for c in self.counts]

    def _deficit(self, counts: Counter[str]) -> float:
        n = len(self.seqs)
        comp = {b: counts[b] / n for b in DNA}
        return _position_deficit(position_signal(comp, self.chemistry), MIN_SIGNAL, MAX_DARK)

    @property
    def total(self) -> float:
        return sum(self.pos)

    def gain(self, out: str, into: str) -> float:
        """Deficit decrease when ``out`` is replaced by ``into``."""
        total = 0.0
        for p, (a, b) in enumerate(zip(out, into, strict=True)):
            if a != b:
                c = self.counts[p].copy()
                c[a] -= 1
                c[b] += 1
                total += self.pos[p] - self._deficit(c)
        return total

    def replace(self, i: int, into: str) -> None:
        out = self.seqs[i]
        for p, (a, b) in enumerate(zip(out, into, strict=True)):
            if a != b:
                self.counts[p][a] -= 1
                self.counts[p][b] += 1
                self.pos[p] = self._deficit(self.counts[p])
        self.seqs[i] = into


class _Fit:
    """Which barcodes of a target set a spare may replace without breaking ``>= d``."""

    def __init__(
        self,
        others: list[str],
        target: list[str],
        metric: Metric,
        d: int,
        reserved: Sequence[str],
        revcomp: bool,
    ) -> None:
        self.metric, self.d, self.revcomp, self.target = metric, d, revcomp, target
        self.fwd = self._index([*reserved, *others])  # never removed: reserved + other sets
        self.rev = self._index([reverse_complement(b) for b in others] if revcomp else [])
        self.tfwd = self._index(target)
        self.trev = self._index([reverse_complement(b) for b in target] if revcomp else [])

    def _index(self, seqs: list[str]) -> CloseIndex:
        index = CloseIndex(self.metric, self.d)
        for s in seqs:
            index.add(s)
        return index

    def _clashes(self, s: str, rc: str, fwd: CloseIndex, rev: CloseIndex) -> bool:
        # As pool.greedy: s vs b, s vs rc(b) and rc(s) vs b (not equal under seqlev).
        if fwd.has_close(s):
            return True
        return self.revcomp and (rev.has_close(s) or fwd.has_close(rc))

    def _close(self, s: str, rc: str, b: str) -> bool:
        m, cut = self.metric, self.d - 1
        if m(s, b, score_cutoff=cut) <= cut:
            return True
        if not self.revcomp:
            return False
        rc_b = reverse_complement(b)
        return m(s, rc_b, score_cutoff=cut) <= cut or m(rc, b, score_cutoff=cut) <= cut

    def replaceable(self, s: str) -> range | list[int]:
        """Indices in the target set that ``s`` may replace (empty if none)."""
        rc = reverse_complement(s) if self.revcomp else s
        if self._clashes(s, rc, self.fwd, self.rev):
            return []
        if self.revcomp and self.metric(s, rc, score_cutoff=self.d - 1) < self.d:
            return []
        if not self._clashes(s, rc, self.tfwd, self.trev):
            return range(len(self.target))
        close = [i for i, b in enumerate(self.target) if self._close(s, rc, b)]
        return close if len(close) == 1 else []  # s may only be close to the barcode it replaces


def repair(
    sets: Sequence[Sequence[str]],
    spares: Sequence[str],
    cfg: Config,
    metric: Metric,
    d: int,
    rng: Random,
    *,
    reserved: Sequence[str] = (),
    revcomp: bool = False,
    max_steps: int = 10_000,
    sample: int = 32,
) -> list[list[str]]:
    """Return a copy of ``sets`` with barcodes replaced by spares to reduce channel deficit.

    Each step scans spares in random order until ``sample`` usable improving
    spares are found and applies the largest deficit reduction. Sets for which
    no improving move exists are returned as they are.
    """
    out = [list(s) for s in sets]
    if cfg.chemistry is None:
        return out
    state = [_Set(s, cfg.chemistry) for s in out]
    pool = list(dict.fromkeys(spares))
    used = {b for s in out for b in s}
    stuck: set[int] = set()
    for _ in range(max_steps):
        open_ = [k for k, s in enumerate(state) if s.total > _EPS and k not in stuck]
        if not open_:
            break
        t = max(open_, key=lambda k: state[k].total)
        target = state[t]
        others = [b for k, s in enumerate(out) if k != t for b in s]
        fit = _Fit(others, out[t], metric, d, reserved, revcomp)
        best: tuple[float, int, str] | None = None
        found = 0
        for j in rng.sample(range(len(pool)), len(pool)):
            s = pool[j]
            if s in used:
                continue
            moves = [(target.gain(out[t][i], s), i) for i in fit.replaceable(s)]
            gain, i = max(moves, key=lambda m: m[0], default=(0.0, -1))
            if gain <= _EPS:
                continue
            if best is None or gain > best[0] + _EPS:
                best = (gain, i, s)
            found += 1
            if found == sample:
                break
        if best is None:
            stuck.add(t)
            continue
        _, i, s = best
        used.discard(out[t][i])
        used.add(s)
        target.replace(i, s)
        stuck.clear()  # a removed barcode may free spares for other sets
    return out
