"""Command-line interface: parse options, run the design, write outputs."""

from pathlib import Path
from typing import Annotated

import typer

from barcode_forge import __version__
from barcode_forge.config import DISABLEABLE, Config, from_layers, load_toml
from barcode_forge.design.pipeline import run
from barcode_forge.report.stats import compute
from barcode_forge.report.writer import write_sets, write_stats

app = typer.Typer(help="Design sets of DNA barcodes for sequencing.", no_args_is_help=True)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", callback=_version_callback, is_eager=True, help="Show version."),
    ] = False,
) -> None:
    """Design sets of DNA barcodes for sequencing."""


def _build_config(config_file: Path | None, disable: list[str], **flags: object) -> Config:
    """Defaults < config file < command-line flags; disabled rules from both."""
    file_values = dict(load_toml(config_file)) if config_file else {}
    disabled = [*file_values.pop("disable", ()), *disable]
    return from_layers(file_values, flags, disable=disabled)


@app.command()
def generate(
    output: Annotated[Path, typer.Option("--output", "-o", help="Output directory.")],
    n_sets: Annotated[int | None, typer.Option("--sets", "-n", help="Number of sets.")] = None,
    set_size: Annotated[
        int | None, typer.Option("--set-size", "-s", help="Barcodes per set.")
    ] = None,
    length: Annotated[int | None, typer.Option("--length", "-l", help="Barcode length.")] = None,
    min_distance: Annotated[
        int | None,
        typer.Option("--min-distance", "-d", help="Minimum distance [default: maximise]."),
    ] = None,
    metric: Annotated[
        str | None, typer.Option(help="hamming | levenshtein | seqlev [default: seqlev].")
    ] = None,
    gc_min: Annotated[float | None, typer.Option(help="Minimum GC fraction.")] = None,
    gc_max: Annotated[float | None, typer.Option(help="Maximum GC fraction.")] = None,
    max_homopolymer: Annotated[int | None, typer.Option(help="Longest homopolymer.")] = None,
    tm_tolerance: Annotated[float | None, typer.Option(help="Tm window ± °C.")] = None,
    avoid: Annotated[
        list[str] | None, typer.Option(help="Motif to avoid on both strands (repeatable).")
    ] = None,
    flank5: Annotated[str | None, typer.Option(help="Sequence 5' of the barcode.")] = None,
    flank3: Annotated[str | None, typer.Option(help="Sequence 3' of the barcode.")] = None,
    disable: Annotated[
        list[str] | None,
        typer.Option(help=f"Rule to switch off (repeatable): {', '.join(DISABLEABLE)}."),
    ] = None,
    chemistry: Annotated[
        str | None,
        typer.Option(help="Channel checks: 2colour | 4colour | none [default: 2colour]."),
    ] = None,
    revcomp: Annotated[
        bool | None,
        typer.Option(
            "--revcomp/--no-revcomp",
            help="Also keep reverse complements >= d apart (non-directional libraries).",
        ),
    ] = None,
    seed: Annotated[int | None, typer.Option(help="Random seed [default: random].")] = None,
    max_candidates: Annotated[
        int | None, typer.Option(help="Rule-passing candidates to screen.")
    ] = None,
    balance: Annotated[
        bool | None, typer.Option("--balance/--no-balance", help="Balance base composition.")
    ] = None,
    config_file: Annotated[
        Path | None, typer.Option("--config", exists=True, dir_okay=False, help="TOML config.")
    ] = None,
) -> None:
    """Generate barcode sets: every barcode in every set is at least d edits apart."""
    try:
        cfg = _build_config(
            config_file,
            disable or [],
            n_sets=n_sets,
            set_size=set_size,
            length=length,
            min_distance=min_distance,
            metric=metric,
            gc_min=gc_min,
            gc_max=gc_max,
            max_homopolymer=max_homopolymer,
            tm_tolerance=tm_tolerance,
            avoid_motifs=tuple(avoid) if avoid else None,
            flank5=flank5,
            flank3=flank3,
            chemistry=chemistry,
            check_revcomp=revcomp,
            seed=seed,
            max_candidates=max_candidates,
            balance=balance,
        )
        design = run(cfg)
    except ValueError as err:  # includes InfeasibleError and invalid config
        typer.echo(f"Error: {err}", err=True)
        raise typer.Exit(code=1) from None

    stats = compute(
        design.sets,
        design.cfg,
        min_distance=design.min_distance,
        rejections=design.rejections,
        candidates_screened=design.candidates_screened,
        seed=design.seed,
    )
    write_sets(output, design.sets)
    write_stats(output / "stats.txt", stats, design.cfg)
    typer.echo(
        f"Wrote {cfg.n_sets} x {cfg.set_size} barcodes (length {cfg.length}, "
        f"min {cfg.metric} distance {design.min_distance}, seed {design.seed}) to {output}"
    )
