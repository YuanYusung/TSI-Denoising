# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Internal data models and signal-processing primitives for TSI denoising.

The helpers build validated station-triplet contexts, enforce common numeric
options, and apply the taper/filter operations used consistently by the core
algorithm and diagnostic visualizations.

Notes
-----
This module is private.  Its arrays are sampled correlation functions with
time in s, distance in km, velocity in km/s, and frequency in Hz.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from obspy import Trace
from scipy.fft import fft, ifft, next_fast_len

from .._validation import station_number
from ..wavefield import Wavefield


EPSILON = 1e-12
DEFAULT_WINDOW_PADDING = 0.2
DEFAULT_SIGNAL_VMIN = 0.2
DEFAULT_SIGNAL_VMAX = 2.0
DEFAULT_GAUSSIAN_ALPHA = 20.0


@dataclass
class _Context:
    """Validated immutable-by-convention geometry used by a TSI operation.

    Attributes hold the input wavefield, station rank, pair lookup, pair
    distances in km, and symmetric correlation times in s.
    """
    wavefield: Wavefield
    station_order: tuple[str, ...]
    rank: dict[str, int]
    pair_indices: dict[tuple[str, str], int]
    distances: np.ndarray
    times: np.ndarray


@dataclass
class _PairDetails:
    """Intermediate candidate products retained for one target station pair.

    ``rows`` has one normalized candidate waveform per station rank; ``tmin``
    and ``tmax`` delimit the velocity-derived signal window in seconds.
    """
    output: np.ndarray
    target_index: int
    target_pair: tuple[str, str]
    tmin: float
    tmax: float
    rows: np.ndarray
    first_inputs: dict[int, np.ndarray]
    second_inputs: dict[int, np.ndarray]
    mechanisms: dict[int, str]
    selected: set[int]


def _require_bool(value, name: str) -> bool:
    """Return ``value`` as bool or raise a named type error.

    NumPy boolean scalars are accepted; integer truth values are intentionally
    rejected so algorithm switches cannot be set accidentally.
    """
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a boolean")
    return bool(value)


def _peak_normalize(data: np.ndarray) -> np.ndarray:
    """Scale a finite waveform by its absolute peak without mutating it.

    Zero-amplitude input returns an equally shaped zero array, avoiding a
    divide-by-zero during candidate construction.
    """
    values = np.asarray(data, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("waveform contains non-finite samples")
    peak = float(np.max(np.abs(values)))
    if peak <= EPSILON:
        return np.zeros_like(values)
    return values / peak


def _validate_common_parameters(
    *,
    distance_threshold: float,
    signal_vmin: float,
    signal_vmax: float,
    window_padding: float,
) -> tuple[float, float, float, float]:
    """Validate shared TSI window options and return finite float values.

    ``distance_threshold`` and ``window_padding`` use km and s; the ordered
    velocity bounds use km/s.
    """
    distance_threshold = float(distance_threshold)
    signal_vmin = float(signal_vmin)
    signal_vmax = float(signal_vmax)
    window_padding = float(window_padding)
    if not np.isfinite(distance_threshold) or distance_threshold < 0:
        raise ValueError("distance_threshold must be finite and non-negative")
    if not np.isfinite(signal_vmin) or not np.isfinite(signal_vmax):
        raise ValueError("signal velocity bounds must be finite")
    if not 0 < signal_vmin < signal_vmax:
        raise ValueError("signal velocity bounds must satisfy 0 < vmin < vmax")
    if not np.isfinite(window_padding) or window_padding < 0:
        raise ValueError("window_padding must be finite and non-negative")
    return distance_threshold, signal_vmin, signal_vmax, window_padding


def _canonical_pair(
    first: str,
    second: str,
    rank: dict[str, int],
) -> tuple[str, str]:
    """Order two known distinct stations according to the context rank."""
    if first not in rank or second not in rank:
        raise ValueError(f"unknown station pair ({first!r}, {second!r})")
    if first == second:
        raise ValueError("a station pair must contain two different stations")
    return (first, second) if rank[first] < rank[second] else (second, first)


def _coerce_station_pair(
    station_pair,
    context: _Context,
) -> tuple[str, str]:
    """Validate, trim, canonicalize, and locate a user station-pair request."""
    if isinstance(station_pair, (str, bytes)):
        raise TypeError(
            "station_pair must be a two-element sequence, not a string"
        )
    try:
        values = tuple(station_pair)
    except TypeError as exc:
        raise TypeError(
            "station_pair must be an iterable containing two station names"
        ) from exc
    if len(values) != 2:
        raise ValueError("station_pair must contain exactly two station names")
    if not all(isinstance(value, str) for value in values):
        raise TypeError("station_pair entries must be strings")
    first, second = (value.strip() for value in values)
    if not first or not second:
        raise ValueError("station_pair must contain non-empty station names")
    pair = _canonical_pair(first, second, context.rank)
    if pair not in context.pair_indices:
        raise ValueError(f"wavefield does not contain station pair {pair}")
    return pair


def _build_context(wavefield: Wavefield) -> _Context:
    """Build TSI pair lookups after validating zero-centered time geometry.

    The method requires at least three stations, uniformly sampled traces, and
    a correlation-time axis symmetric around zero lag.
    """
    if not isinstance(wavefield, Wavefield):
        raise TypeError("wavefield must be a Wavefield")

    pairs = wavefield.pairs
    stations = tuple(
        sorted(
            {station for pair in pairs for station in pair},
            key=station_number,
        )
    )
    if len(stations) < 3:
        raise ValueError(
            "three-station interferometry requires at least three stations"
        )
    times = np.asarray(wavefield.time, dtype=float)
    tolerance = max(1e-9, wavefield.delta * 1e-6)
    if not np.allclose(
        np.diff(times),
        wavefield.delta,
        rtol=0.0,
        atol=tolerance,
    ):
        raise ValueError(
            "TSI requires a uniformly sampled correlation-time axis"
        )
    zero_index = int(np.argmin(np.abs(times)))
    if not np.isclose(
        times[zero_index],
        0.0,
        rtol=0.0,
        atol=tolerance,
    ):
        raise ValueError("TSI correlation-time axis must contain zero lag")
    if not np.allclose(
        times,
        -times[::-1],
        rtol=0.0,
        atol=tolerance,
    ):
        raise ValueError(
            "TSI correlation-time axis must be symmetric about zero lag"
        )

    rank = {station: index for index, station in enumerate(stations)}
    pair_indices: dict[tuple[str, str], int] = {}
    for index, (first, second) in enumerate(pairs):
        pair = _canonical_pair(first, second, rank)
        if pair in pair_indices:
            raise ValueError(
                "wavefield contains duplicate station pairs after station-order "
                f"canonicalization: {pair}"
            )
        pair_indices[pair] = index
    return _Context(
        wavefield=wavefield,
        station_order=stations,
        rank=rank,
        pair_indices=pair_indices,
        distances=wavefield.distances,
        times=times,
    )


def _sqrt_amplitude_spectrum(data: np.ndarray) -> np.ndarray:
    """Apply square-root amplitude spectral whitening while retaining phase.

    A fast doubled FFT reduces circular effects; the returned real waveform
    has the same one-dimensional length as ``data``.
    """
    spectrum = fft(data, n=next_fast_len(data.size * 2))
    root_spectrum = np.sqrt(np.abs(spectrum)) * np.exp(
        1j * np.angle(spectrum)
    )
    return np.real(ifft(root_spectrum))[: data.size]


def _taper_and_filter(
    data: np.ndarray,
    times: np.ndarray,
    *,
    delta: float,
    tmin: float,
    tmax: float,
    fmin: float,
    fmax: float,
) -> np.ndarray:
    """Apply the optional legacy signal-window taper and bandpass filter."""
    keep = (times >= tmin) & (times <= tmax)
    if np.count_nonzero(keep) < 2:
        return np.zeros_like(data)
    output = np.zeros_like(data)
    segment = Trace(np.asarray(data[keep], dtype=float).copy())
    segment.stats.delta = delta
    segment.taper(max_percentage=0.05, type="cosine", side="both")
    output[keep] = segment.data
    filtered = Trace(output)
    filtered.stats.delta = delta
    filtered.filter(
        "bandpass",
        freqmin=fmin,
        freqmax=fmax,
        zerophase=True,
    )
    return np.asarray(filtered.data, dtype=float)
