"""Write barcode sets (``set_NN.tsv``) and the run report (``stats.txt``)."""

import dataclasses
from collections.abc import Sequence
from pathlib import Path

from barcode_forge import __version__
from barcode_forge.config import Config
from barcode_forge.report.stats import Stats
from barcode_forge.sequence import DNA


def _width(n_sets: int) -> int:
    return max(2, len(str(n_sets)))


def write_sets(outdir: Path, sets: Sequence[Sequence[str]]) -> list[Path]:
    """One headerless single-column file per set: set_01.tsv ... set_NN.tsv."""
    outdir.mkdir(parents=True, exist_ok=True)
    width = _width(len(sets))
    paths = []
    for i, group in enumerate(sets, 1):
        path = outdir / f"set_{i:0{width}d}.tsv"
        path.write_text("".join(f"{s}\n" for s in group))
        paths.append(path)
    return paths


def _fmt(value: object, spec: str = "") -> str:
    return "n/a" if value is None else format(value, spec)


def _section(title: str, rows: list[str]) -> list[str]:
    return ["", title, "-" * len(title), *rows]


def _composition(set_no: str, composition: list[dict[str, float]], imbalance: float) -> list[str]:
    header = "      " + "".join(f"{p:>5d}" for p in range(1, len(composition) + 1))
    rows = [f"  {b} % " + "".join(f"{100 * pos[b]:5.0f}" for pos in composition) for b in DNA]
    return [f"set {set_no} (imbalance {imbalance:.4f})", header, *rows]


def _chemistry(chemistry: str | None, violations: list[list[str]], width: int) -> list[str]:
    if chemistry is None:
        return ["chemistry = none (channel balance not checked)"]
    rows = [f"chemistry = {chemistry}"]
    for i, problems in enumerate(violations, 1):
        if problems:
            noun = "violation" if len(problems) == 1 else "violations"
            rows += [f"set {i:0{width}d}: {len(problems)} {noun}", *(f"  {p}" for p in problems)]
        else:
            rows.append(f"set {i:0{width}d}: OK")
    return rows


def format_stats(stats: Stats, cfg: Config) -> str:
    """Deterministic plain-text report: parameters to reproduce the run, then results."""
    params = dataclasses.asdict(cfg) | {"seed": stats.seed}
    lines = [f"barcode-forge {__version__} report"]
    lines += _section("Parameters", [f"{k} = {v!r}" for k, v in params.items()])
    lines += _section(
        "Distances",
        [
            f"metric = {stats.metric}",
            f"achieved min distance = {stats.min_distance}",
            f"sets = {stats.n_sets} x {stats.set_size} barcodes of length {stats.length}",
            f"within-set min / mean = {_fmt(stats.within_min)} / {_fmt(stats.within_mean, '.2f')}",
            f"cross-set min / mean = {_fmt(stats.cross_min)} / {_fmt(stats.cross_mean, '.2f')}",
            f"overall min = {_fmt(stats.overall_min)}",
        ],
    )
    lines += _section(
        "Sequence properties",
        [
            f"GC = {stats.gc_min:.1%} to {stats.gc_max:.1%}",
            f"max homopolymer = {stats.max_homopolymer}",
            f"max dinucleotide repeat = {stats.max_dinucleotide_repeat}",
            f"max palindrome = {stats.max_palindrome}",
        ],
    )
    lines += _section(
        "Thermodynamics",
        [
            f"Tm = {stats.tm_min:.1f} to {stats.tm_max:.1f} C",
            f"worst hairpin dG = {_fmt(stats.worst_hairpin_dg, '.2f')} kcal/mol",
            f"worst homodimer dG = {_fmt(stats.worst_homodimer_dg, '.2f')} kcal/mol",
        ],
    )
    rejected = [f"{name} = {n}" for name, n in stats.rejections.items()] or ["none"]
    lines += _section(
        "Rule rejections", [f"candidates screened = {stats.candidates_screened}", *rejected]
    )
    width = _width(stats.n_sets)
    lines += _section(
        "Sequencing chemistry", _chemistry(cfg.chemistry, stats.channel_violations, width)
    )
    composition = []
    for i, (comp, imb) in enumerate(zip(stats.composition, stats.imbalance, strict=True), 1):
        composition += _composition(f"{i:0{width}d}", comp, imb)
    lines += _section("Base composition per position (%)", composition)
    return "\n".join(lines) + "\n"


def write_stats(path: Path, stats: Stats, cfg: Config) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(format_stats(stats, cfg))
    return path
