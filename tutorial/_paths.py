#!/Users/yusong/miniforge3/bin/python
"""Resolve tutorial case paths independently of the notebook working directory."""

from __future__ import annotations

from pathlib import Path


def _candidate_case_directories(case_name: str, cwd: Path):
    """Yield plausible case directories from a root or nested working path."""
    current = cwd.expanduser().resolve()
    for base in (current, *current.parents):
        if base.name == case_name:
            yield base
        if base.name == "tutorial":
            yield base / case_name
        yield base / "tutorial" / case_name


def resolve_case_paths(
    case_name: str,
    *,
    cwd: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    """Return ``(case_dir, input_dir, processed_dir)`` for a tutorial case.

    The case is located from the current directory or any of its parents, so
    notebooks work when launched from the repository root, ``tutorial/``, or
    the individual case directory. The input directory is intentionally not
    created here; the caller can provide a clear download instruction when it
    is absent.
    """
    if not isinstance(case_name, str) or not case_name.strip():
        raise ValueError("case_name must be a non-empty string")

    case_name = case_name.strip()
    working_directory = Path.cwd() if cwd is None else Path(cwd)
    seen: set[Path] = set()
    for case_dir in _candidate_case_directories(case_name, working_directory):
        case_dir = case_dir.resolve()
        if case_dir in seen:
            continue
        seen.add(case_dir)
        notebook = case_dir / f"run_example_{case_name}.ipynb"
        if notebook.is_file():
            return case_dir, case_dir / "input_public", case_dir / "processed"

    raise FileNotFoundError(
        f"could not locate tutorial case {case_name!r} from "
        f"{working_directory.resolve()}"
    )


__all__ = ["resolve_case_paths"]
