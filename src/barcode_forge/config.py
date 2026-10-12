"""Run configuration: defaults, validation and TOML loading.

Defaults target inline/cell barcodes amplified by PCR and read in the first
cycles of a short-read run. Setting a rule parameter to ``None`` disables it.
"""

import dataclasses
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from barcode_forge.sequence import validate

METRICS = ("hamming", "levenshtein", "seqlev")

# Rules that can be switched off, and the parameter that controls each.
DISABLEABLE = {
    "homopolymer": "max_homopolymer",
    "repeats": "max_dinucleotide_repeat",
    "palindrome": "max_palindrome",
    "hairpin": "min_hairpin_dg",
    "homodimer": "min_homodimer_dg",
    "tm": "tm_tolerance",
}


@dataclass(frozen=True)
class Config:
    # Size of the design
    n_sets: int = 1
    set_size: int = 96
    length: int = 10

    # Distance: None means maximise the minimum distance achievable
    min_distance: int | None = None
    metric: str = "levenshtein"

    # Sequence composition rules
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_homopolymer: int | None = 2
    max_dinucleotide_repeat: int | None = 2  # copies in a row, e.g. 2 allows ATAT not ATATAT
    max_palindrome: int | None = 4  # longest reverse-complement palindrome allowed (nt)
    avoid_motifs: tuple[str, ...] = ()  # checked on both strands

    # Thermodynamics (primer3; kcal/mol, °C, mM / nM)
    min_hairpin_dg: float | None = -3.0
    min_homodimer_dg: float | None = -6.0
    # Opt-in (e.g. when the barcode itself primes): an embedded barcode never anneals
    # alone, so its Tm is not a PCR criterion; GC already evens amplification.
    tm_tolerance: float | None = None  # ± °C around tm_target
    tm_target: float | None = None  # None: estimated from candidates
    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.2
    dna_conc: float = 250.0

    # Sequence context the barcode sits in (adapter / linker)
    flank5: str = ""
    flank3: str = ""

    # Search
    seed: int | None = None
    max_candidates: int = 50_000
    balance: bool = True

    def __post_init__(self) -> None:
        if self.n_sets < 1 or self.set_size < 1:
            raise ValueError("n_sets and set_size must be >= 1")
        if self.length < 4:
            raise ValueError("length must be >= 4")
        if not 0.0 <= self.gc_min <= self.gc_max <= 1.0:
            raise ValueError("need 0 <= gc_min <= gc_max <= 1")
        if self.min_distance is not None and not 1 <= self.min_distance <= self.length:
            raise ValueError("min_distance must be between 1 and length")
        if self.metric not in METRICS:
            raise ValueError(f"metric must be one of {METRICS}")
        for name in ("max_homopolymer", "max_dinucleotide_repeat", "max_palindrome"):
            value = getattr(self, name)
            if value is not None and value < 1:
                raise ValueError(f"{name} must be >= 1")
        if self.max_candidates < 1:
            raise ValueError("max_candidates must be >= 1")
        object.__setattr__(self, "avoid_motifs", tuple(validate(m) for m in self.avoid_motifs))
        for name in ("flank5", "flank3"):
            value = getattr(self, name)
            object.__setattr__(self, name, validate(value) if value else "")

    @property
    def total(self) -> int:
        return self.n_sets * self.set_size

    def with_overrides(self, **overrides: Any) -> "Config":
        """Return a copy with every non-None override applied."""
        return dataclasses.replace(self, **{k: v for k, v in overrides.items() if v is not None})

    def without(self, *rules: str) -> "Config":
        """Return a copy with the named rules disabled."""
        unknown = sorted(set(rules) - set(DISABLEABLE))
        if unknown:
            raise ValueError(f"cannot disable {unknown}; choose from {sorted(DISABLEABLE)}")
        return dataclasses.replace(self, **{DISABLEABLE[r]: None for r in rules})


def load_toml(path: Path) -> dict[str, Any]:
    """Read config overrides from a TOML file, rejecting unknown keys.

    Besides Config fields, ``disable = [...]`` lists rules to switch off
    (TOML has no null); apply it with ``Config.without``.
    """
    data = tomllib.loads(Path(path).read_text())
    known = {f.name for f in dataclasses.fields(Config)} | {"disable"}
    unknown = sorted(set(data) - known)
    if unknown:
        raise ValueError(f"unknown config keys: {', '.join(unknown)}")
    return {k: tuple(v) if isinstance(v, list) else v for k, v in data.items()}
