"""
asap-o-matic

A *lightly* altered version of the asap_to_kite_v2.py script
written by Caleb Lareau and
found at https://github.com/caleblareau/asap_to_kite
"""
from importlib.metadata import PackageNotFoundError, version
from typing import Annotated

import typer
from loguru import logger
from rich import print as rp

# from rich.console import Console

try:
    __version__ = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"

logger.disable("asap-o-matic")

app = typer.Typer(
    name="asap_o_matic",
    help="Reformat antibody-derived reads from ASAP-seq to a format expected by CITE-seq-Count or Bustools",
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="markdown",
)

verbosity_level = 0


def version_callback(value: bool) -> None:  # FBT001
    """Prints the version of the package."""
    if value:
        rp(f"[yellow]asap-to-kite[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


@app.callback()
def verbosity(
    verbose: Annotated[
        int,
        typer.Option(
            "-v",
            "--verbose",
            help="Control output verbosity. Pass this argument multiple times to increase the amount of output.",
            count=True,
        ),
    ] = 0
) -> None:
    verbosity_level = verbose  # noqa: F841


if __name__ == "main":
    app()
