# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Core algorithms and persistence objects for iterative TSI denoising.

For each target station pair, candidate correlations from compatible station
triplets are aligned and selected in a velocity-constrained signal window.
The iterative public entry point returns a new wavefield together with
convergence diagnostics; it never changes the input wavefield.

Notes
-----
Distance is km, velocity is km/s, time is s, and filter bounds are Hz.  The
result archive comprises a final-wavefield NPZ plus a diagnostics NPZ.

References
----------
Qiu, H., Niu, F., and Qin, L. (2021). *JGR: Solid Earth*, 126,
e2021JB021712. https://doi.org/10.1029/2021JB021712
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import tempfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve

from .._validation import validate_frequency_band, validate_n_jobs
from ..wavefield import Wavefield
from ._common import (
    DEFAULT_SIGNAL_VMAX,
    DEFAULT_SIGNAL_VMIN,
    DEFAULT_WINDOW_PADDING,
EPSILON,
    _Context,
    _PairDetails,
    _build_context,
    _canonical_pair,
    _coerce_station_pair,
    _peak_normalize,
    _require_bool,
    _sqrt_amplitude_spectrum,
    _taper_and_filter,
    _validate_common_parameters,
)


_DENOISING_RESULT_FORMAT = "surface-wave-denoising-result"
_DENOISING_RESULT_FORMAT_VERSION = 1


def _result_paths(base_path: str | Path) -> tuple[Path, Path]:
    """Return ``_wavefield`` and ``_info`` NPZ paths for an output stem.

    The parent directory of ``base_path`` is the caller-selected output
    directory (normally ``processed``); its final path component is the result
    name shared by both companion files.
    """
    base = Path(base_path)
    return Path(f"{base}_wavefield.npz"), Path(f"{base}_info.npz")


def _scalar(archive, name: str):
    """Extract a scalar NPZ field and reject arrays with an invalid shape."""
    value = np.asarray(archive[name])
    if value.shape != ():
        raise ValueError(f"{name} must be a scalar")
    return value.item()


@dataclass(frozen=True)
class DenoisingResult:
    """Products and convergence diagnostics from iterative denoising.

    ``example_history`` has shape ``(iterations, n_samples)`` and contains the
    selected station pair after each completed iteration. ``relative_changes``
    contains the corresponding whole-wavefield relative L2 changes.

    Parameters
    ----------
    final_wavefield : Wavefield
        New normalized output wavefield from the final iteration.
    example_pair : tuple[str, str]
        Canonically ordered pair represented by ``example_history``.
    example_history : ndarray
        Array of shape ``(iterations, n_samples)``.
    relative_changes : ndarray
        Whole-wavefield relative L2 change after each iteration.
    converged : bool
        Whether the requested change threshold stopped the iteration.
    stop_reason : {"threshold", "max_iterations"}
        Recorded termination condition.

    Notes
    -----
    :meth:`save` writes ``<name>_wavefield.npz`` and ``<name>_info.npz`` to
    the parent directory of its output stem and never overwrites either unless
    requested explicitly.
    """

    final_wavefield: Wavefield
    example_pair: tuple[str, str]
    example_history: np.ndarray
    relative_changes: np.ndarray
    converged: bool
    stop_reason: str

    @property
    def iterations(self) -> int:
        """Number of completed denoising iterations."""
        return int(self.relative_changes.size)

    def save(self, base_path: str | Path, *, overwrite: bool = False) -> tuple[Path, Path]:
        """Persist the final wavefield and diagnostics as two NPZ files.

        Parameters
        ----------
        base_path : path-like
            Output stem.  Its parent is the destination directory (typically
            ``processed``); ``_wavefield.npz`` and ``_info.npz`` are appended
            to its final component.
        overwrite : bool, default=False
            Replace existing destination files only when true.

        Returns
        -------
        tuple[pathlib.Path, pathlib.Path]
            Paths of the wavefield archive and diagnostics archive.
        """
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean")
        wavefield_path, metadata_path = _result_paths(base_path)
        if not metadata_path.parent.is_dir():
            raise FileNotFoundError(
                f"parent directory does not exist: {metadata_path.parent}"
            )
        if not overwrite:
            for path in (wavefield_path, metadata_path):
                if path.exists():
                    raise FileExistsError(f"denoising result file already exists: {path}")

        history = np.asarray(self.example_history, dtype=float)
        changes = np.asarray(self.relative_changes, dtype=float)
        if history.ndim != 2 or history.shape != (
            changes.size,
            self.final_wavefield.n_samples,
        ):
            raise ValueError("example_history has inconsistent shape")
        if changes.ndim != 1 or not np.all(np.isfinite(history)) or not np.all(np.isfinite(changes)):
            raise ValueError("denoising diagnostics must be finite arrays")
        if (
            not isinstance(self.example_pair, tuple)
            or len(self.example_pair) != 2
            or not all(isinstance(name, str) and name for name in self.example_pair)
        ):
            raise ValueError("example_pair must contain two non-empty station names")
        if self.stop_reason not in {"threshold", "max_iterations"}:
            raise ValueError("unsupported denoising stop reason")

        payload = {
            "format_name": np.asarray(_DENOISING_RESULT_FORMAT),
            "format_version": np.asarray(_DENOISING_RESULT_FORMAT_VERSION, dtype=np.int64),
            "example_pair": np.asarray(self.example_pair, dtype=np.str_),
            "example_history": history,
            "relative_changes": changes,
            "converged": np.asarray(self.converged, dtype=np.bool_),
            "stop_reason": np.asarray(self.stop_reason, dtype=np.str_),
        }
        descriptor, temporary_name = tempfile.mkstemp(
            dir=metadata_path.parent,
            prefix=f".{metadata_path.name}.",
            suffix=".tmp",
        )
        os.close(descriptor)
        temporary_path = Path(temporary_name)
        try:
            with temporary_path.open("wb") as file:
                np.savez_compressed(file, **payload)
                file.flush()
                os.fsync(file.fileno())
            self.final_wavefield.save(wavefield_path, overwrite=overwrite)
            if overwrite:
                os.replace(temporary_path, metadata_path)
            else:
                os.link(temporary_path, metadata_path)
                temporary_path.unlink()
        finally:
            if temporary_path.exists():
                temporary_path.unlink()
        return wavefield_path, metadata_path

    @classmethod
    def load(cls, base_path: str | Path) -> "DenoisingResult":
        """Load and validate a two-file result archive written by :meth:`save`.

        Raises
        ------
        FileNotFoundError
            If either companion archive is absent.
        ValueError
            If the schema, array shapes, finite values, or format version are
            invalid.
        """
        wavefield_path, metadata_path = _result_paths(base_path)
        if not wavefield_path.exists():
            raise FileNotFoundError(f"denoised wavefield file does not exist: {wavefield_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"denoising result file does not exist: {metadata_path}")
        required = {
            "format_name", "format_version", "example_pair", "example_history",
            "relative_changes", "converged", "stop_reason",
        }
        try:
            final_wavefield = Wavefield.load(wavefield_path)
            with np.load(metadata_path, allow_pickle=False) as archive:
                missing = sorted(required - set(archive.files))
                if missing:
                    raise ValueError(f"missing required fields: {missing}")
                if _scalar(archive, "format_name") != _DENOISING_RESULT_FORMAT:
                    raise ValueError("unsupported denoising result format")
                if _scalar(archive, "format_version") != _DENOISING_RESULT_FORMAT_VERSION:
                    raise ValueError("unsupported denoising result format version")
                example_pair = np.asarray(archive["example_pair"])
                if example_pair.shape != (2,) or example_pair.dtype.kind not in {"U", "S"}:
                    raise ValueError("example_pair must contain two station names")
                history = np.asarray(archive["example_history"], dtype=float)
                changes = np.asarray(archive["relative_changes"], dtype=float)
                if history.shape != (changes.size, final_wavefield.n_samples):
                    raise ValueError("example_history has inconsistent shape")
                if changes.ndim != 1 or not np.all(np.isfinite(history)) or not np.all(np.isfinite(changes)):
                    raise ValueError("denoising diagnostics must be finite arrays")
                stop_reason = str(_scalar(archive, "stop_reason"))
                if stop_reason not in {"threshold", "max_iterations"}:
                    raise ValueError("unsupported denoising stop reason")
                return cls(
                    final_wavefield=final_wavefield,
                    example_pair=tuple(example_pair.astype(str).tolist()),
                    example_history=history.copy(),
                    relative_changes=changes.copy(),
                    converged=bool(_scalar(archive, "converged")),
                    stop_reason=stop_reason,
                )
        except FileNotFoundError:
            raise
        except Exception as exc:
            raise ValueError(f"invalid denoising result {metadata_path}: {exc}") from exc


def _denoise_pair(
    context: _Context,
    data: np.ndarray,
    station_pair: tuple[str, str],
    *,
    include_convolution: bool,
    sqrt_spectrum: bool,
    taper_output: bool,
    fmin: float,
    fmax: float,
    distance_threshold: float,
    signal_vmin: float,
    signal_vmax: float,
    window_padding: float,
) -> _PairDetails:
    """Denoise one pair and retain diagnostics needed by the demo plot."""
    first, second = station_pair
    target_index = context.pair_indices[station_pair]
    target = _peak_normalize(data[target_index])
    distance = context.distances[target_index]
    tmin = distance / signal_vmax - window_padding
    tmax = distance / signal_vmin + window_padding
    nstations = len(context.station_order)
    rows = np.full((nstations, data.shape[1]), np.nan, dtype=float)
    first_inputs: dict[int, np.ndarray] = {}
    second_inputs: dict[int, np.ndarray] = {}
    mechanisms: dict[int, str] = {}
    candidates: list[np.ndarray] = []
    candidate_stations: list[int] = []

    first_rank, second_rank = context.rank[first], context.rank[second]
    for third_rank, third in enumerate(context.station_order):
        if third == first or third == second:
            continue
        if third_rank < first_rank:
            first_pair = _canonical_pair(third, first, context.rank)
            second_pair = _canonical_pair(third, second, context.rank)
            is_cross_correlation = True
            reverse_first = True
        elif third_rank > second_rank:
            first_pair = _canonical_pair(first, third, context.rank)
            second_pair = _canonical_pair(second, third, context.rank)
            is_cross_correlation = True
            reverse_first = False
        else:
            first_pair = _canonical_pair(first, third, context.rank)
            second_pair = _canonical_pair(third, second, context.rank)
            is_cross_correlation = False
            reverse_first = False

        first_index = context.pair_indices.get(first_pair)
        second_index = context.pair_indices.get(second_pair)
        if first_index is None or second_index is None:
            continue
        if (
            context.distances[first_index] < distance_threshold
            or context.distances[second_index] < distance_threshold
        ):
            continue
        if not is_cross_correlation and not include_convolution:
            continue

        first_input = _peak_normalize(data[first_index])
        second_input = _peak_normalize(data[second_index])
        if reverse_first:
            interferogram = fftconvolve(
                second_input,
                first_input[::-1],
                mode="same",
            )
        elif is_cross_correlation:
            interferogram = fftconvolve(
                first_input,
                second_input[::-1],
                mode="same",
            )
        else:
            interferogram = fftconvolve(
                first_input,
                second_input,
                mode="same",
            )
        if sqrt_spectrum:
            interferogram = _sqrt_amplitude_spectrum(interferogram)
        interferogram = _peak_normalize(interferogram)
        if not np.any(interferogram):
            continue
        rows[third_rank] = interferogram
        first_inputs[third_rank] = first_input
        second_inputs[third_rank] = second_input
        mechanisms[third_rank] = (
            "cross_correlation" if is_cross_correlation else "convolution"
        )
        candidates.append(interferogram)
        candidate_stations.append(third_rank)

    if len(candidates) < 2:
        return _PairDetails(
            output=target,
            target_index=target_index,
            target_pair=station_pair,
            tmin=tmin,
            tmax=tmax,
            rows=rows,
            first_inputs=first_inputs,
            second_inputs=second_inputs,
            mechanisms=mechanisms,
            selected=set(),
        )

    signal_mask = (context.times > tmin) & (context.times < tmax)
    if np.count_nonzero(signal_mask) < 2:
        return _PairDetails(
            output=target,
            target_index=target_index,
            target_pair=station_pair,
            tmin=tmin,
            tmax=tmax,
            rows=rows,
            first_inputs=first_inputs,
            second_inputs=second_inputs,
            mechanisms=mechanisms,
            selected=set(),
        )

    candidate_array = np.asarray(candidates)
    raw_mean = _peak_normalize(np.mean(candidate_array, axis=0))
    reference = raw_mean
    if taper_output:
        reference = _peak_normalize(
            _taper_and_filter(
                reference,
                context.times,
                delta=context.wavefield.delta,
                tmin=tmin,
                tmax=tmax,
                fmin=fmin,
                fmax=fmax,
            )
        )

    coefficients = np.empty(candidate_array.shape[0], dtype=float)
    for index, candidate in enumerate(candidate_array):
        numerator = float(
            np.dot(candidate[signal_mask], reference[signal_mask])
        )
        denominator = float(
            np.linalg.norm(candidate[signal_mask])
            * np.linalg.norm(reference[signal_mask])
        )
        coefficients[index] = (
            numerator / denominator if denominator > EPSILON else -np.inf
        )

    finite = np.isfinite(coefficients)
    good = np.zeros(coefficients.size, dtype=bool)
    if np.any(finite):
        q25, q75 = np.percentile(coefficients[finite], (25, 75))
        cutoff = max(0.4, float(q25 - 1.25 * (q75 - q25)))
        good[finite] = coefficients[finite] > cutoff
    selected = {
        candidate_stations[index] for index in np.flatnonzero(good)
    }

    if np.count_nonzero(good) < 2:
        output = target
    else:
        output = _peak_normalize(np.mean(candidate_array[good], axis=0))
        if taper_output:
            output = _peak_normalize(
                _taper_and_filter(
                    output,
                    context.times,
                    delta=context.wavefield.delta,
                    tmin=tmin,
                    tmax=tmax,
                    fmin=fmin,
                    fmax=fmax,
                )
            )
    return _PairDetails(
        output=output,
        target_index=target_index,
        target_pair=station_pair,
        tmin=tmin,
        tmax=tmax,
        rows=rows,
        first_inputs=first_inputs,
        second_inputs=second_inputs,
        mechanisms=mechanisms,
        selected=selected,
    )


def _wavefield_from_data(
    wavefield: Wavefield,
    data: np.ndarray,
) -> Wavefield:
    """Construct an independent wavefield from validated output samples.

    ``data`` must have shape ``(wavefield.n_pairs, wavefield.n_samples)``;
    copied trace arrays preserve all input metadata and component identity.
    """
    values = np.asarray(data)
    expected_shape = (wavefield.n_pairs, wavefield.n_samples)
    if values.shape != expected_shape:
        raise ValueError(
            f"output data must have shape {expected_shape}; got {values.shape}"
        )
    stream = wavefield.stream()
    for trace, row in zip(stream, values):
        trace.data = np.asarray(row, dtype=float).copy()
    return Wavefield(
        stream,
        component=wavefield.component,
        copy=False,
        check_distance_order=wavefield.check_distance_order,
    )


def _denoise_pair_output(
    station_pair: tuple[str, str],
    *,
    context: _Context,
    data: np.ndarray,
    include_convolution: bool,
    sqrt_spectrum: bool,
    taper_output: bool,
    fmin: float,
    fmax: float,
    distance_threshold: float,
    signal_vmin: float,
    signal_vmax: float,
    window_padding: float,
) -> tuple[int, np.ndarray]:
    """Return the output row for one target pair, suitable for Pool.map."""
    details = _denoise_pair(
        context,
        data,
        station_pair,
        include_convolution=include_convolution,
        sqrt_spectrum=sqrt_spectrum,
        taper_output=taper_output,
        fmin=fmin,
        fmax=fmax,
        distance_threshold=distance_threshold,
        signal_vmin=signal_vmin,
        signal_vmax=signal_vmax,
        window_padding=window_padding,
    )
    return details.target_index, details.output


def denoise_wavefield_iteratively(
    wavefield: Wavefield,
    example_pair,
    *,
    threshold: float,
    first_iteration_convolution: bool,
    max_iterations: int = 6,
    sqrt_spectrum: bool = True,
    taper_output: bool = False,
    fmin: float = 0.5,
    fmax: float = 5.0,
    distance_threshold: float = 0.0,
    signal_vmin: float = DEFAULT_SIGNAL_VMIN,
    signal_vmax: float = DEFAULT_SIGNAL_VMAX,
    window_padding: float = DEFAULT_WINDOW_PADDING,
    n_jobs: int = 1,
) -> DenoisingResult:
    """Iteratively denoise every trace in a one-component wavefield.

    The first iteration optionally includes inner-station convolutions through
    the required ``first_iteration_convolution`` flag. Every later iteration
    uses both outer-station cross correlations and inner-station convolutions.
    Iteration stops when the whole-wavefield relative L2 change is no larger
    than ``threshold`` or after ``max_iterations`` completed iterations.
    """
    context = _build_context(wavefield)
    pair = _coerce_station_pair(example_pair, context)
    first_iteration_convolution = _require_bool(
        first_iteration_convolution,
        "first_iteration_convolution",
    )
    sqrt_spectrum = _require_bool(sqrt_spectrum, "sqrt_spectrum")
    taper_output = _require_bool(taper_output, "taper_output")
    n_jobs = validate_n_jobs(n_jobs)
    threshold = float(threshold)
    if not np.isfinite(threshold) or threshold < 0:
        raise ValueError("threshold must be finite and non-negative")
    if isinstance(max_iterations, bool) or not isinstance(
        max_iterations,
        (int, np.integer),
    ):
        raise TypeError("max_iterations must be an integer")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1")
    (
        distance_threshold,
        signal_vmin,
        signal_vmax,
        window_padding,
    ) = _validate_common_parameters(
        distance_threshold=distance_threshold,
        signal_vmin=signal_vmin,
        signal_vmax=signal_vmax,
        window_padding=window_padding,
    )
    if taper_output:
        fmin, fmax = validate_frequency_band(
            fmin,
            fmax,
            wavefield.sampling_rate,
        )

    current = np.asarray(
        [_peak_normalize(row) for row in wavefield.data()]
    )
    example_index = context.pair_indices[pair]
    history: list[np.ndarray] = []
    changes: list[float] = []
    converged = False
    target_pairs = tuple(context.pair_indices)
    pool = Pool(processes=n_jobs) if n_jobs > 1 else None
    try:
        for iteration in range(int(max_iterations)):
            include_convolution = (
                first_iteration_convolution if iteration == 0 else True
            )
            compute_output = partial(
                _denoise_pair_output,
                context=context,
                data=current,
                include_convolution=include_convolution,
                sqrt_spectrum=sqrt_spectrum,
                taper_output=taper_output,
                fmin=fmin,
                fmax=fmax,
                distance_threshold=distance_threshold,
                signal_vmin=signal_vmin,
                signal_vmax=signal_vmax,
                window_padding=window_padding,
            )
            outputs = (
                [compute_output(target_pair) for target_pair in target_pairs]
                if pool is None
                else pool.map(compute_output, target_pairs)
            )
            next_data = np.empty_like(current)
            for target_index, output in outputs:
                next_data[target_index] = output
            change = float(
                np.linalg.norm(next_data - current)
                / max(float(np.linalg.norm(current)), EPSILON)
            )
            history.append(next_data[example_index].copy())
            changes.append(change)
            current = next_data
            if change <= threshold:
                converged = True
                break
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    return DenoisingResult(
        final_wavefield=_wavefield_from_data(wavefield, current),
        example_pair=pair,
        example_history=np.asarray(history),
        relative_changes=np.asarray(changes),
        converged=converged,
        stop_reason="threshold" if converged else "max_iterations",
    )


__all__ = [
    "DenoisingResult",
    "denoise_wavefield_iteratively",
]
