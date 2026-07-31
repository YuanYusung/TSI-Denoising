# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Read SAC surface-wave cross-correlations into validated wavefields.

Files are discovered recursively, read one trace at a time, normalized to
ascending station-pair order, and checked through :class:`Wavefield`.

Notes
-----
SAC ``kevnm``/``kstnm`` and ``dist`` metadata provide pair identity and
distance (km); the correlation-time origin is obtained from SAC ``b``.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from obspy import Stream, read

from .._validation import normalize_pair_directions, trace_pair
from ..wavefield import Wavefield


DEFAULT_PATTERNS = ("*.[Ss][Aa][Cc]*",)


def _find_sac_files(directory: Path, pattern: str | Iterable[str] | None):
    """Recursively find unique SAC paths matching one or more wildcards."""
    patterns = DEFAULT_PATTERNS if pattern is None else pattern
    if isinstance(patterns, str):
        patterns = (patterns,)
    else:
        patterns = tuple(patterns)
    if not patterns:
        raise ValueError("pattern must contain at least one wildcard")

    files = {
        path
        for current_pattern in patterns
        for path in directory.rglob(current_pattern)
        if path.is_file()
    }
    return sorted(files)


def read_sac_directory(
    directory: str | Path,
    pattern: str | Iterable[str] | None = None,
    *,
    component: str | None = None,
) -> Wavefield:
    """Read all matching SAC files under a directory into a ``Wavefield``.

    The default pattern is case-insensitive for ``SAC`` and also matches names
    such as ``.SAC_s``. Custom patterns may be passed as one wildcard or an
    iterable of wildcards. Files are searched recursively and de-duplicated.
    Traces are sorted by their resolved source-receiver names. When omitted,
    the component name is inferred from the directory name.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"SAC directory does not exist: {directory}")

    files = _find_sac_files(directory, pattern)
    if not files:
        raise FileNotFoundError(f"no SAC files found in {directory}")

    stream = Stream()
    for path in files:
        try:
            loaded = read(str(path))
        except Exception as exc:
            raise ValueError(f"failed to read SAC file {path}") from exc
        if len(loaded) != 1:
            raise ValueError(f"SAC file must contain exactly one trace: {path}")
        stream += loaded

    stream = normalize_pair_directions(stream)
    stream.traces.sort(key=trace_pair)
    component_name = directory.name if component is None else str(component).strip()
    if not component_name:
        raise ValueError("component must be a non-empty name")
    return Wavefield(stream, component=component_name, copy=False)

__all__ = ["read_sac_directory"]
