"""Command-line interface."""

import typer

from barcode_forge import __version__

app = typer.Typer(help="Design sets of DNA barcodes for sequencing.", no_args_is_help=True)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(__version__)
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", callback=_version_callback, is_eager=True, help="Show version."
    ),
) -> None:
    """Design sets of DNA barcodes for sequencing."""
