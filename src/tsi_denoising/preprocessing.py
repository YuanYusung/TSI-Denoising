# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Preprocessing for causal surface-wave cross-correlations.

The public routine symmetrizes each correlation, applies a distance-dependent
causal velocity window with a cosine taper, and band-pass filters the result.
It returns a new wavefield; :meth:`Wavefield.preprocess` is the explicitly
in-place convenience method.

Notes
-----
Distance is km, velocity is km/s, frequencies are Hz, and correlation time is
s.  The causal signal window spans ``distance / vmax`` through
``distance / vmin``.
"""

from __future__ import annotations

import numpy as np
from obspy import Stream

from ._validation import (
    trace_distance,
    trace_pair,
    trace_start_time,
    validate_frequency_band,
    validate_stream,
)
from .wavefield import Wavefield


def _velocity_window(trace, vmin: float, vmax: float) -> np.ndarray:
    """Return the causal sample mask for one trace."""
    time = trace_start_time(trace) + np.arange(trace.stats.npts) * trace.stats.delta
    distance = trace_distance(trace)
    return (time >= distance / vmax) & (time <= distance / vmin)


def _window_trace(trace, *, vmin: float, vmax: float, taper_fraction: float):
    """Symmetrize, window, and return data with the original length."""
    symmetrized = 0.5 * (trace.data + trace.data[::-1])
    keep = _velocity_window(trace, vmin, vmax)
    if keep.sum() < 2:
        raise ValueError(
            f"velocity window for trace {trace_pair(trace)} "
            "contains fewer than two samples"
        )

    windowed = trace.copy()
    windowed.data = np.asarray(symmetrized[keep]).copy()
    if taper_fraction > 0:
        windowed.taper(
            max_percentage=taper_fraction,
            type="cosine",
            side="both",
        )

    data = np.zeros_like(trace.data)
    data[keep] = windowed.data
    return data


def _validate_taper_parameters(vmin: float, vmax: float, taper_fraction: float):
    """Validate finite km/s velocity bounds and a cosine-taper fraction."""
    vmin = float(vmin)
    vmax = float(vmax)
    taper_fraction = float(taper_fraction)
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("velocity bounds must be finite")
    if not 0 < vmin < vmax:
        raise ValueError("velocity bounds must satisfy 0 < vmin < vmax")
    if not np.isfinite(taper_fraction) or not 0 <= taper_fraction <= 0.5:
        raise ValueError("taper_fraction must be between 0 and 0.5")
    return vmin, vmax, taper_fraction


def _velocity_taper_stream(
    stream: Stream,
    *,
    vmin: float,
    vmax: float,
    taper_fraction: float,
) -> Stream:
    """Return a copied stream after distance-dependent velocity tapering.

    Velocity bounds are km/s and taper fraction is applied to each retained
    causal window.  No frequency filtering is performed here.
    """
    """Apply the preprocessing velocity-window taper without filtering."""
    validate_stream(stream)
    output = stream.copy()
    for input_trace, output_trace in zip(stream, output):
        output_trace.data = _window_trace(
            input_trace,
            vmin=vmin,
            vmax=vmax,
            taper_fraction=taper_fraction,
        )
    return output


def _apply_velocity_taper(
    wavefield: Wavefield,
    *,
    vmin: float = 0.1,
    vmax: float = 2.5,
    taper_fraction: float = 0.05,
) -> Wavefield:
    """Return a wavefield after the standard distance-dependent taper."""
    if not isinstance(wavefield, Wavefield):
        raise TypeError("wavefield must be a Wavefield")
    vmin, vmax, taper_fraction = _validate_taper_parameters(
        vmin, vmax, taper_fraction
    )
    stream = _velocity_taper_stream(
        wavefield.stream(),
        vmin=vmin,
        vmax=vmax,
        taper_fraction=taper_fraction,
    )
    return Wavefield(
        stream,
        component=wavefield.component,
        copy=False,
        check_distance_order=wavefield.check_distance_order,
    )


def _preprocess_stream(
    stream: Stream,
    *,
    fmin: float,
    fmax: float,
    vmin: float,
    vmax: float,
    taper_fraction: float,
) -> Stream:
    """Return a copied stream after velocity tapering and zero-phase filtering.

    Frequency limits are Hz; velocity bounds are km/s.  The input stream is
    validated and never modified.
    """
    output = _velocity_taper_stream(
        stream,
        vmin=vmin,
        vmax=vmax,
        taper_fraction=taper_fraction,
    )

    output.filter(
        "bandpass",
        freqmin=fmin,
        freqmax=fmax,
        corners=4,
        zerophase=True,
    )
    return output


def preprocess_stream(
    wavefield: Wavefield,
    fmin: float = 0.5,
    fmax: float = 5.0,
    vmin: float = 0.1,
    vmax: float = 2.5,
    taper_fraction: float = 0.05,
) -> Wavefield:
    """Preprocess one Wavefield.

    The causal window spans ``distance / vmax`` to ``distance / vmin``.
    Distances and velocities are expressed in km and km/s, respectively.
    Traces are first symmetrized with their time-reversed copies. A new
    Wavefield is returned and the input is not modified.
    """
    if not isinstance(wavefield, Wavefield):
        raise TypeError("wavefield must be a Wavefield")

    fmin, fmax = validate_frequency_band(fmin, fmax, wavefield.sampling_rate)
    vmin, vmax, taper_fraction = _validate_taper_parameters(
        vmin, vmax, taper_fraction
    )

    stream = _preprocess_stream(
        wavefield.stream(),
        fmin=fmin,
        fmax=fmax,
        vmin=vmin,
        vmax=vmax,
        taper_fraction=taper_fraction,
    )
    return Wavefield(
        stream,
        component=wavefield.component,
        copy=False,
        check_distance_order=wavefield.check_distance_order,
    )


__all__ = ["preprocess_stream"]
