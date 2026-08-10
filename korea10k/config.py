"""Site-specific filesystem locations used by the Korea10K analysis scripts.

Every absolute path that depends on where the data, tools and references
actually live is resolved here from environment variables, so that no
server-specific path is committed to this repository.

Set the variables either by exporting them in your shell or by copying
``.env.example`` to ``.env`` in the repository root and filling it in::

    cp .env.example .env
    $EDITOR .env

A variable is only required by the scripts that actually use it, so you can
fill in just the roots you need.  Importing a name whose variable is unset
raises :class:`MissingPathError` with an explanatory message rather than
silently falling back to a wrong location.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = [
    "CONDA_ENV_DIR",
    "KOREF_DIR",
    "PROJECT_DIR",
    "REF_DIR",
    "STORE_DIR",
    "TOOL_DIR",
    "WORK_DIR",
    "MissingPathError",
    "get_dir",
]

#: Maps the exported constant to (environment variable, what it points at).
_ROOTS: dict[str, tuple[str, str]] = {
    "WORK_DIR": (
        "KOREA10K_WORK_DIR",
        "personal analysis workspace holding per-analysis working directories "
        "(admixture, shapeit, hla, depthCoverage, ...)",
    ),
    "PROJECT_DIR": (
        "KOREA10K_PROJECT_DIR",
        "Korea10K project root holding Resources/, Results/ and Analysis/",
    ),
    "TOOL_DIR": (
        "KOREA10K_TOOL_DIR",
        "shared bioinformatics tool installations (plink2, shapeit5, "
        "admixture, minimac, GLIMPSE2, ...)",
    ),
    "REF_DIR": (
        "KOREA10K_REF_DIR",
        "shared reference genomes and liftover chain files",
    ),
    "STORE_DIR": (
        "KOREA10K_STORE_DIR",
        "bulk storage for raw sequencing and per-sample intermediate data",
    ),
    "CONDA_ENV_DIR": (
        "KOREA10K_CONDA_ENV_DIR",
        "conda environment root providing the external binaries invoked by "
        "the pipeline scripts (tabix, bcftools, ...)",
    ),
    "KOREF_DIR": (
        "KOREF_DIR",
        "KOREF personal multiomics reference project root",
    ),
}

_DOTENV_FILENAME = ".env"


class MissingPathError(RuntimeError):
    """Raised when a required path environment variable is not configured."""


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _load_dotenv() -> None:
    """Populate ``os.environ`` from the repository-root ``.env``, if present.

    Variables already present in the environment take precedence, so an
    explicit ``export`` always overrides the file.
    """
    dotenv = _repo_root() / _DOTENV_FILENAME
    if not dotenv.is_file():
        return
    for raw in dotenv.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip().strip("\"'"))


_load_dotenv()


def get_dir(name: str) -> str:
    """Return the configured directory for one of the roots in :data:`_ROOTS`.

    Raises :class:`MissingPathError` if the backing environment variable is
    unset or empty.
    """
    try:
        env_var, description = _ROOTS[name]
    except KeyError:
        raise KeyError(
            f"{name!r} is not a known path root; expected one of "
            f"{', '.join(sorted(_ROOTS))}"
        ) from None

    value = os.environ.get(env_var, "").strip()
    if not value:
        raise MissingPathError(
            f"{env_var} is not set. It must point at the {description}.\n"
            f"Set it in {_repo_root() / _DOTENV_FILENAME} (see .env.example) "
            f"or export it:\n"
            f"    export {env_var}=/path/to/directory"
        )
    return value.rstrip("/")


def __getattr__(name: str) -> str:
    """Resolve the root constants lazily.

    Resolving on attribute access rather than at import time means a script
    only needs the roots it actually imports to be configured.
    """
    if name in _ROOTS:
        return get_dir(name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(__all__)
