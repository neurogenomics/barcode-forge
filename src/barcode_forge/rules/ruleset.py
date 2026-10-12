"""Rule type and an ordered, rejection-counting collection of rules."""

from collections import Counter
from collections.abc import Callable, Iterable
from dataclasses import dataclass

from barcode_forge.config import Config


@dataclass(frozen=True)
class Rule:
    name: str
    cost: int  # relative evaluation cost; cheap rules run first
    check: Callable[[str, Config], bool]  # True means the sequence passes


class RuleSet:
    def __init__(self, rules: Iterable[Rule], cfg: Config) -> None:
        self.rules = tuple(sorted(rules, key=lambda r: r.cost))
        self.cfg = cfg
        self.rejections: Counter[str] = Counter()

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(r.name for r in self.rules)

    def passes(self, seq: str) -> bool:
        for rule in self.rules:
            if not rule.check(seq, self.cfg):
                self.rejections[rule.name] += 1
                return False
        return True
