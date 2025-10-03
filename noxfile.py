"""Nox sessions."""

from pathlib import Path

import nox
import nox_uv as nu

PACKAGE = "asap_o_matic"
PYTHON_VERSIONS = ["3.11", "3.12", "3.13", "3.14"]
PRIMARY_PYTHON_VERSION = "3.12"

nox.needs_version = ">=2025.5.1"
nox.options.sessions = (
    "typecheck",
    "lint",
    "test",
)

nox.options.default_venv_backend = "uv"

locations = (
    "python",
    "tests",
)

@nu.session(python=PRIMARY_PYTHON_VERSION, uv_groups=["type_check"], uv_sync_locked=True)
def typecheck(s: nox.Session) -> None:
    """Type-check using mypy."""
    _ = s.run(
        "ty",
        "check",
        "python/"
    )

@nu.session(python=PRIMARY_PYTHON_VERSION, uv_groups=["lint"], uv_sync_locked=True)
def lint(s: nox.Session) -> None:
    """Type-check using mypy."""
    _ = s.run(
        "ruff",
        "check",
        "--fix",
        "python"
    )

@nu.session(python=PYTHON_VERSIONS, uv_groups=["test"], uv_sync_locked=True)
def test(s: nox.Session) -> None:
    """Run the test suite."""
    s.install(".")
    _ = s.run(
        "pytest",
        "--cov=asap_o_matic",
        "tests",
    )


@nu.session(python=PRIMARY_PYTHON_VERSION, uv_groups=["coverage"], uv_sync_locked=True)
def coverage(s: nox.Session) -> None:
    """Produce the coverage report."""
    args = s.posargs or ["report"]
    s.install("coverage[toml]", "codecov", external=True)
    s.install(".")
    if not s.posargs and any(Path().glob(".coverage.*")):
        _ = s.run("coverage", "combine")

    _ = s.run("coverage", "json", "--fail-under=0")
    _ = s.run("codecov", *args)
