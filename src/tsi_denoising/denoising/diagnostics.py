# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Diagnostic visualizations for three-station-interferometry denoising.

The module constructs reproducible, Figure-4-style views of input gathers,
TSI candidates, narrow-band products, and iterative denoising results without
mutating the supplied :class:`~tsi_denoising.Wavefield`.

Notes
-----
Plot time is seconds and station-pair distance is km.  Matplotlib figures are
returned to the caller and are not saved implicitly.
"""

from __future__ import annotations

import numpy as np
from scipy.fft import fft, ifft, next_fast_len
from scipy.signal import hilbert

from .._validation import validate_frequency_band
from ..wavefield import Wavefield, _local_distance_spacing, _stable_jitter
from ._common import (
    DEFAULT_GAUSSIAN_ALPHA,
    DEFAULT_SIGNAL_VMIN,
    EPSILON,
    _build_context,
    _canonical_pair,
    _coerce_station_pair,
    _peak_normalize,
    _require_bool,
    _validate_common_parameters,
)
from .three_station import DenoisingResult, _denoise_pair


_MECHANISM_COLORS = {
    "cross_correlation": "k",
    "convolution": "b",
    "target": "r",
    "original": "gray",
    "output": "k",
}
_FONT_SIZES = {
    "title": 10,
    "axis": 12,
    "tick": 12,
    "panel": 14,
    "annotation": 12,
}
_TICK_RANGE_FRACTION = 0.2
_SPECTRUM_INSET_BOUNDS = (0.58, 0.055, 0.38, 0.32)
_SPECTRUM_INSET_RELATIVE_HALF_WIDTH = 0.35
_SPECTRUM_INSET_MIN_BINS = 4
_DENOISED_RESULT_FIGURE_WIDTH = 12.0
_DENOISED_RESULT_ROW_HEIGHT = 4.5


def _coerce_periods(periods) -> np.ndarray:
    """Return a finite one-dimensional positive period array in seconds."""
    values = np.asarray(periods, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(
            "periods must be a non-empty one-dimensional array"
        )
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("periods must contain finite positive values")
    return values


def _gaussian_narrowband(
    data: np.ndarray,
    *,
    delta: float,
    periods: np.ndarray,
) -> np.ndarray:
    """Return the RR-workflow Gaussian narrow-band output for every trace."""
    values = np.asarray(data, dtype=float)
    delta = float(delta)
    if values.ndim != 2:
        raise ValueError(
            "narrow-band data must have shape (n_traces, n_samples)"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("narrow-band data contain non-finite samples")
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError("delta must be finite and positive")
    periods = _coerce_periods(periods)
    center_frequencies = 1.0 / periods
    nyquist = 0.5 / delta
    if np.any(center_frequencies >= nyquist):
        raise ValueError(
            "narrow-band center frequencies must be below "
            f"Nyquist ({nyquist:g} Hz)"
        )

    npts = values.shape[-1]
    nfft = next_fast_len(npts * 2)
    spectrum = fft(values, n=nfft, axis=-1)
    angular_step = 2.0 * np.pi / (nfft * delta)
    angular_frequency = angular_step * np.arange(nfft)
    output = np.empty(
        (values.shape[0], npts, periods.size),
        dtype=float,
    )
    for index, period in enumerate(periods):
        center = 2.0 * np.pi / period
        window = np.exp(
            -((angular_frequency - center) / center) ** 2
            * DEFAULT_GAUSSIAN_ALPHA
        )
        filtered = spectrum * window
        filtered[:, nfft // 2 + 1 :] = 0.0
        filtered[:, 0] /= 2.0
        filtered[:, nfft // 2] = np.real(filtered[:, nfft // 2])
        output[:, :, index] = np.real(ifft(filtered, axis=-1))[:, :npts]
    return output


def _gaussian_narrowband_allowing_missing_rows(
    data: np.ndarray,
    *,
    delta: float,
    periods: np.ndarray,
) -> np.ndarray:
    """Filter finite rows while preserving all-NaN denoising placeholders."""
    values = np.asarray(data, dtype=float)
    if values.ndim != 2:
        raise ValueError(
            "narrow-band data must have shape (n_traces, n_samples)"
        )
    finite_rows = np.all(np.isfinite(values), axis=1)
    missing_rows = np.all(np.isnan(values), axis=1)
    if not np.all(finite_rows | missing_rows):
        raise ValueError(
            "narrow-band rows must be finite or entirely NaN"
        )
    filled = values.copy()
    filled[missing_rows] = 0.0
    output = _gaussian_narrowband(
        filled,
        delta=delta,
        periods=periods,
    )
    output[missing_rows] = np.nan
    return output


def _station_title(station: str) -> str:
    """Return the station identifier exactly as supplied by the wavefield."""
    return station


def _panel_label(index: int) -> str:
    """Return spreadsheet-style lower-case panel labels: (a), ..., (aa)."""
    letters = ""
    value = index
    while True:
        value, remainder = divmod(value, 26)
        letters = chr(ord("a") + remainder) + letters
        if value == 0:
            return f"({letters})"
        value -= 1


def _proportional_ticks(start: float, stop: float) -> np.ndarray:
    """Return labels separated by 20% of an axis's displayed data range."""
    if not np.isfinite(start) or not np.isfinite(stop) or stop < start:
        raise ValueError("tick range must be finite and increasing")
    if np.isclose(start, stop):
        return np.asarray([start], dtype=float)
    step = (stop - start) * _TICK_RANGE_FRACTION
    count = int(np.ceil((stop - start) / step)) + 1
    return np.linspace(start, stop, count)


def _demo_station_geometry(context):
    """Return station coordinates, scale, limits, and real station ticks.

    Tick targets are first spaced proportionally across the displayed station
    range, then snapped to the nearest actual station.  This avoids labels
    such as ``10.2`` when a sparse or non-consecutive station subset is used.
    """
    try:
        positions = np.asarray(context.station_order, dtype=float)
    except (TypeError, ValueError):
        positions = np.arange(1, len(context.station_order) + 1, dtype=float)
    steps = np.diff(positions)
    scale = 0.9 * (float(np.median(np.abs(steps))) if steps.size else 1.0)
    ymin = float(np.min(positions) - 2.0 * scale)
    ymax = float(np.max(positions) + 2.0 * scale)
    target_ticks = _proportional_ticks(
        float(np.min(positions)),
        float(np.max(positions)),
    )
    tick_indices = np.unique(
        [int(np.argmin(np.abs(positions - target))) for target in target_ticks]
    )
    ticks = positions[tick_indices]
    labels = tuple(context.station_order[index] for index in tick_indices)
    return positions, steps, scale, ymin, ymax, ticks, labels


def _diagnostic_distance_ticks(context, station_positions: np.ndarray):
    """Return the right-hand distance-axis label positions for the figure."""
    reference_station = context.station_order[0]
    distances = np.zeros(len(context.station_order), dtype=float)
    known_positions = [0]
    for rank, station in enumerate(context.station_order[1:], start=1):
        pair = _canonical_pair(reference_station, station, context.rank)
        trace_index = context.pair_indices.get(pair)
        if trace_index is None:
            distances[rank] = np.nan
        else:
            distances[rank] = context.distances[trace_index]
            known_positions.append(rank)
    if len(known_positions) > 1:
        known = np.asarray(known_positions)
        distances = np.interp(np.arange(len(distances)), known, distances[known])
    else:
        distances = np.arange(len(distances), dtype=float)
    if np.any(np.diff(distances) <= 0):
        raise ValueError("diagnostic distance axis must increase with station number")

    tick_values = _proportional_ticks(0.0, float(distances[-1]))
    tick_positions = np.interp(tick_values, distances, station_positions)
    return reference_station, tick_positions, tick_values


def _filter_input_gather(
    inputs: dict[int, np.ndarray],
    *,
    nstations: int,
    npts: int,
    delta: float,
    periods: np.ndarray,
) -> np.ndarray:
    """Narrow-band filter every input gather once for all requested periods."""
    rows = np.zeros((nstations, npts), dtype=float)
    present = np.zeros(nstations, dtype=bool)
    for rank, values in inputs.items():
        rows[rank] = values
        present[rank] = True
    filtered = _gaussian_narrowband(rows, delta=delta, periods=periods)
    filtered[~present] = np.nan
    return filtered


def _narrowband_diagnostic_products(
    details,
    data: np.ndarray,
    context,
    pair: tuple[str, str],
    *,
    delta: float,
    periods: np.ndarray,
):
    """Precompute all narrow-band data used by the diagnostic figure."""
    first, second = pair
    nstations, npts = details.rows.shape
    rows = np.nan_to_num(details.rows, nan=0.0)
    rows[context.rank[first]] = data[details.target_index]
    rows[context.rank[second]] = data[details.target_index]
    filtered_rows = _gaussian_narrowband(rows, delta=delta, periods=periods)
    missing = np.isnan(details.rows).all(axis=1)
    missing[context.rank[first]] = False
    missing[context.rank[second]] = False
    filtered_rows[missing] = np.nan
    return (
        _filter_input_gather(
            details.first_inputs,
            nstations=nstations,
            npts=npts,
            delta=delta,
            periods=periods,
        ),
        _filter_input_gather(
            details.second_inputs,
            nstations=nstations,
            npts=npts,
            delta=delta,
            periods=periods,
        ),
        filtered_rows,
        _gaussian_narrowband(
            data[details.target_index][None, :], delta=delta, periods=periods
        )[0],
        _gaussian_narrowband(
            details.output[None, :], delta=delta, periods=periods
        )[0],
    )


def _plot_input_gather(
    axis,
    filtered: np.ndarray,
    *,
    input_ranks: dict[int, np.ndarray],
    mechanisms: dict[int, str],
    times: np.ndarray,
    station_positions: np.ndarray,
    scale: float,
    period_index: int,
):
    """Draw one virtual-shot input gather from precomputed narrow-band data."""
    for rank in input_ranks:
        color = _MECHANISM_COLORS[mechanisms[rank]]
        axis.plot(
            times,
            scale * _peak_normalize(filtered[rank, :, period_index])
            + station_positions[rank],
            color=color,
            linewidth=0.6,
        )


def _normalized_interferograms(interferograms: np.ndarray) -> np.ndarray:
    """Peak-normalize each candidate row while preserving missing rows as NaN."""
    finite = np.where(np.isfinite(interferograms), interferograms, 0.0)
    peaks = np.max(np.abs(finite), axis=1, keepdims=True)
    return np.divide(
        interferograms,
        peaks,
        out=np.full_like(interferograms, np.nan),
        where=peaks > EPSILON,
    )


def _signal_window(output: np.ndarray, details, times: np.ndarray, period: float):
    """Return the four-period display window centred on the output envelope."""
    if np.max(np.abs(output)) <= EPSILON:
        peak_time = float(np.clip(0.5 * (details.tmin + details.tmax), times[0], times[-1]))
    else:
        peak_time = times[int(np.argmax(np.abs(hilbert(output))))]
    return float(peak_time - 2.0 * period), float(peak_time + 2.0 * period)


def _jittered_distances(
    wavefield: Wavefield,
    *,
    enabled: bool,
) -> np.ndarray:
    """Return distances with the Wavefield plot's stable duplicate jitter."""
    distances = wavefield.distances
    if not enabled:
        return distances
    plotted = distances.copy()
    for distance in np.unique(distances):
        indices = np.flatnonzero(distances == distance)
        if len(indices) < 2:
            continue
        spacing = _local_distance_spacing(distances, float(distance))
        if spacing is None:
            continue
        limit = 0.2 * spacing
        for index in indices:
            plotted[index] += _stable_jitter(
                wavefield.pairs[index],
                float(distance),
                limit,
            )
    return plotted


def _normalized_rows(data: np.ndarray) -> np.ndarray:
    """Peak-normalize every waveform while retaining finite zero rows."""
    values = np.asarray(data, dtype=float)
    peak = np.max(np.abs(values), axis=1, keepdims=True)
    return np.divide(values, peak, out=np.zeros_like(values), where=peak > EPSILON)


def _mean_normalized_spectrum(data: np.ndarray, delta: float):
    """Return positive-frequency mean amplitude spectrum and its peak index."""
    values = np.asarray(data, dtype=float)
    spectrum = np.mean(np.abs(np.fft.rfft(values, axis=-1)), axis=0)
    frequency = np.fft.rfftfreq(values.shape[-1], d=delta)
    peak = float(np.max(spectrum))
    if peak <= EPSILON:
        return frequency, np.zeros_like(spectrum), None
    normalized = spectrum / peak
    positive = frequency > 0
    if not np.any(positive):
        return frequency, normalized, None
    candidates = np.flatnonzero(positive)
    return frequency, normalized, int(candidates[np.argmax(normalized[candidates])])


def _centered_spectrum_limits(
    frequency: np.ndarray,
    peak_frequency: float,
) -> tuple[float, float] | None:
    """Return a compact symmetric frequency window centered on a spectral peak.

    The target half-width is 35% of the peak frequency, with a four-bin
    minimum so a low-frequency peak remains legible. It is reduced only when
    available positive-frequency samples would clip a side of the red peak
    marker, so the marker remains at the horizontal midpoint of the inset.
    """
    values = np.asarray(frequency, dtype=float)
    if values.ndim != 1 or values.size < 2:
        return None
    lower_bound = float(values[1])
    upper_bound = float(values[-1])
    frequency_step = float(np.min(np.diff(values)))
    target_half_width = max(
        _SPECTRUM_INSET_RELATIVE_HALF_WIDTH * peak_frequency,
        _SPECTRUM_INSET_MIN_BINS * frequency_step,
    )
    half_width = min(
        target_half_width,
        peak_frequency - lower_bound,
        upper_bound - peak_frequency,
    )
    if not np.isfinite(half_width) or half_width <= EPSILON:
        return None
    return peak_frequency - half_width, peak_frequency + half_width


def _annotate_tsi_panel(
    axis,
    *,
    display_tmax: float,
    target_first: float,
    target_second: float,
    station_positions: np.ndarray,
    station_steps: np.ndarray,
    ymin: float,
    scale: float,
    signal_tmin: float,
    signal_tmax: float,
):
    """Add mechanism and signal-window annotations to the TSI panel."""
    offset = 0.25 * (
        float(np.median(np.abs(station_steps))) if station_steps.size else 1.0
    )
    annotation_x = display_tmax - 0.5
    for first_point, second_point, color in (
        (
            target_first - offset,
            target_second + offset,
            _MECHANISM_COLORS["convolution"],
        ),
        (
            np.min(station_positions),
            target_first + offset,
            _MECHANISM_COLORS["cross_correlation"],
        ),
        (
            target_second - offset,
            np.max(station_positions),
            _MECHANISM_COLORS["cross_correlation"],
        ),
    ):
        axis.annotate(
            "",
            xy=(annotation_x, first_point),
            xytext=(annotation_x, second_point),
            arrowprops=dict(
                color=color,
                arrowstyle="<->",
                mutation_scale=16,
                lw=1.4,
            ),
        )
    text_x = annotation_x - 0.1
    axis.text(
        text_x,
        (target_first + target_second) / 2.0,
        "Convolution",
        rotation=90,
        ha="right",
        va="center",
        fontsize=_FONT_SIZES["annotation"],
        color=_MECHANISM_COLORS["convolution"],
    )
    axis.text(
        text_x,
        (target_first + np.min(station_positions)) / 2.0,
        "CC",
        rotation=90,
        ha="right",
        va="center",
        fontsize=_FONT_SIZES["annotation"],
    )
    axis.text(
        text_x,
        (target_second + np.max(station_positions)) / 2.0,
        "CC",
        rotation=90,
        ha="right",
        va="center",
        fontsize=_FONT_SIZES["annotation"],
    )
    signal_y = ymin + 0.9 * scale
    axis.annotate(
        "",
        xy=(signal_tmin, signal_y),
        xytext=(signal_tmax, signal_y),
        arrowprops=dict(
            color=_MECHANISM_COLORS["cross_correlation"],
            arrowstyle="<->",
            mutation_scale=14,
            lw=1.5,
        ),
    )


def denoise_station_pair_demo(
    wavefield: Wavefield,
    station_pair,
    *,
    include_convolution: bool = True,
    sqrt_spectrum: bool = True,
    taper_output: bool = False,
    fmin: float = 0.5,
    fmax: float = 4.0,
    window_padding: float = 0.1,
    periods=(0.8, 0.3),
    time_limits: tuple[float, float] = (-2.0, 8.0),
):
    """Plot one three-station-interferometry denoising diagnostic.

    The Figure-4-style layout has one three-panel row per Gaussian narrow-band
    period. ``taper_output`` is disabled by default: it controls only the TSI
    candidate-stack taper/bandpass, not preprocessing of ``wavefield``.
    """
    import matplotlib.pyplot as plt

    context = _build_context(wavefield)
    pair = _coerce_station_pair(station_pair, context)
    include_convolution = _require_bool(
        include_convolution,
        "include_convolution",
    )
    sqrt_spectrum = _require_bool(sqrt_spectrum, "sqrt_spectrum")
    taper_output = _require_bool(taper_output, "taper_output")
    periods = _coerce_periods(periods)
    try:
        display_tmin, display_tmax = (
            float(value) for value in time_limits
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("time_limits must contain two finite values") from exc
    if (
        not np.isfinite(display_tmin)
        or not np.isfinite(display_tmax)
        or display_tmin >= display_tmax
    ):
        raise ValueError("time_limits must contain two increasing finite values")
    _, internal_vmin, internal_vmax, window_padding = _validate_common_parameters(
        distance_threshold=0.0,
        signal_vmin=DEFAULT_SIGNAL_VMIN,
        signal_vmax=1.5,
        window_padding=window_padding,
    )
    if taper_output:
        fmin, fmax = validate_frequency_band(
            fmin,
            fmax,
            wavefield.sampling_rate,
        )

    data = np.asarray([_peak_normalize(row) for row in wavefield.data()])
    details = _denoise_pair(
        context,
        data,
        pair,
        include_convolution=include_convolution,
        sqrt_spectrum=sqrt_spectrum,
        taper_output=taper_output,
        fmin=fmin,
        fmax=fmax,
        distance_threshold=0.0,
        signal_vmin=internal_vmin,
        signal_vmax=internal_vmax,
        window_padding=window_padding,
    )

    first, second = pair
    (
        station_positions,
        station_steps,
        scale,
        ymin,
        ymax,
        station_ticks,
        station_tick_labels,
    ) = _demo_station_geometry(context)
    target_first = station_positions[context.rank[first]]
    target_second = station_positions[context.rank[second]]
    (
        filtered_first,
        filtered_second,
        filtered_rows,
        filtered_target,
        filtered_output,
    ) = _narrowband_diagnostic_products(
        details,
        data,
        context,
        pair,
        delta=wavefield.delta,
        periods=periods,
    )
    (
        reference_station,
        distance_tick_positions,
        distance_tick_values,
    ) = _diagnostic_distance_ticks(
        context,
        station_positions,
    )
    time_grid, station_grid = np.meshgrid(context.times, station_positions)

    fig, axes = plt.subplots(
        periods.size,
        3,
        figsize=(10, 4.25 * periods.size),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for period_index, (period, row_axes) in enumerate(zip(periods, axes)):
        target = _peak_normalize(filtered_target[:, period_index])
        output = _peak_normalize(filtered_output[:, period_index])
        _plot_input_gather(
            row_axes[0],
            filtered_first,
            input_ranks=details.first_inputs,
            mechanisms=details.mechanisms,
            times=context.times,
            station_positions=station_positions,
            scale=scale,
            period_index=period_index,
        )
        _plot_input_gather(
            row_axes[1],
            filtered_second,
            input_ranks=details.second_inputs,
            mechanisms=details.mechanisms,
            times=context.times,
            station_positions=station_positions,
            scale=scale,
            period_index=period_index,
        )
        row_axes[0].plot(
            context.times,
            scale * target + target_second,
            f"{_MECHANISM_COLORS['target']}-",
            linewidth=0.8,
        )
        row_axes[1].plot(
            context.times,
            scale * target + target_first,
            f"{_MECHANISM_COLORS['target']}-",
            linewidth=0.8,
        )
        for axis, star_position in (
            (row_axes[0], target_first),
            (row_axes[1], target_second),
        ):
            axis.plot(
                0.0,
                star_position,
                f"{_MECHANISM_COLORS['target']}*",
                markersize=13,
                markeredgecolor="k",
                clip_on=False,
            )

        interferograms = _normalized_interferograms(filtered_rows[:, :, period_index])
        row_axes[2].pcolormesh(
            time_grid,
            station_grid,
            interferograms,
            cmap="coolwarm",
            rasterized=True,
            shading="gouraud",
            vmin=-1,
            vmax=1,
        )
        overlay_scale = (ymax - ymin) / 8.0
        row_axes[2].plot(
            context.times,
            target * overlay_scale + target_first,
            color=_MECHANISM_COLORS["original"],
            linewidth=0.8,
        )
        row_axes[2].plot(
            context.times,
            output * overlay_scale + target_first,
            f"{_MECHANISM_COLORS['output']}-",
            linewidth=0.8,
        )
        signal_tmin, signal_tmax = _signal_window(
            output,
            details,
            context.times,
            period,
        )
        _annotate_tsi_panel(
            row_axes[2],
            display_tmax=display_tmax,
            target_first=target_first,
            target_second=target_second,
            station_positions=station_positions,
            station_steps=station_steps,
            ymin=ymin,
            scale=scale,
            signal_tmin=signal_tmin,
            signal_tmax=signal_tmax,
        )

        for axis_index, axis in enumerate(row_axes):
            if axis_index == 2:
                axis.axvline(
                    signal_tmin,
                    color="k",
                    linestyle="--",
                    linewidth=1.5,
                )
                axis.axvline(
                    signal_tmax,
                    color="k",
                    linestyle="--",
                    linewidth=1.5,
                )
            axis.set_xlim(display_tmin, display_tmax)
            axis.set_ylim(ymin, ymax)
            axis.tick_params(labelsize=_FONT_SIZES["tick"])
            axis.text(
                display_tmin + 0.1,
                ymin + 0.1 * scale,
                _panel_label(3 * period_index + axis_index),
                fontsize=_FONT_SIZES["panel"],
                ha="left",
                va="bottom",
                bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
            )
        row_axes[0].set_ylabel("Station number", fontsize=_FONT_SIZES["axis"])
        row_axes[0].set_yticks(station_ticks)
        row_axes[0].set_yticklabels(station_tick_labels)
        row_axes[0].text(
            display_tmax - 0.1,
            ymin + 0.1 * scale,
            f"{period:.1f}s",
            fontsize=_FONT_SIZES["panel"],
            ha="right",
            va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
        )
        if period_index == 0:
            row_axes[0].set_title(
                f"Virtual shot gather at station {_station_title(first)}",
                fontsize=_FONT_SIZES["title"],
            )
            row_axes[1].set_title(
                f"Virtual shot gather at station {_station_title(second)}",
                fontsize=_FONT_SIZES["title"],
            )
            row_axes[2].set_title(
                "Three-station interferometry",
                fontsize=_FONT_SIZES["title"],
            )
        if period_index == periods.size - 1:
            for axis in row_axes:
                axis.set_xlabel("Correlation time (s)", fontsize=_FONT_SIZES["axis"])
        distance_axis = row_axes[2].twinx()
        distance_axis.set_ylim(ymin, ymax)
        distance_axis.set_yticks(distance_tick_positions)
        distance_axis.set_yticklabels(
            [f"{distance:.1f}" for distance in distance_tick_values]
        )
        distance_axis.tick_params(labelsize=_FONT_SIZES["tick"])
        distance_axis.set_ylabel(
            f"Distance to station {_station_title(reference_station)} (km)",
            fontsize=_FONT_SIZES["axis"],
        )

    fig.subplots_adjust(
        left=0.08,
        right=0.88,
        bottom=0.08,
        top=0.94,
        wspace=0.08,
        hspace=0.03,
    )
    return fig, axes


def plot_denoised_result(
    wavefield: Wavefield,
    result: DenoisingResult,
    *,
    periods=(0.8, 0.3),
    time_limits: tuple[float, float] = (-2.0, 8.0),
    jitter_duplicate_distances: bool = True,
):
    """Plot the narrow-band wavefield before and after iterative denoising.

    One row is drawn per ``period``: the example pair and iteration history
    appear at left, with the original and final denoised wavefields in the
    centre and right panels.  The right panel includes an inset comparison of
    average spectra.  Duplicate distances receive the same stable, small
    vertical jitter used by :meth:`Wavefield.plot` by default.
    """
    import matplotlib.pyplot as plt

    if not isinstance(wavefield, Wavefield):
        raise TypeError("wavefield must be a Wavefield")
    if not isinstance(result, DenoisingResult):
        raise TypeError("result must be a DenoisingResult")
    if not isinstance(jitter_duplicate_distances, bool):
        raise TypeError("jitter_duplicate_distances must be a boolean")
    if result.final_wavefield.pairs != wavefield.pairs:
        raise ValueError("result final_wavefield pairs do not match wavefield")
    try:
        display_tmin, display_tmax = (float(value) for value in time_limits)
    except (TypeError, ValueError) as exc:
        raise ValueError("time_limits must contain two finite values") from exc
    if (
        not np.isfinite(display_tmin)
        or not np.isfinite(display_tmax)
        or display_tmin >= display_tmax
    ):
        raise ValueError("time_limits must contain two increasing finite values")
    periods = _coerce_periods(periods)
    try:
        example_index = wavefield.pairs.index(result.example_pair)
    except ValueError as exc:
        raise ValueError(
            f"wavefield does not contain result example pair {result.example_pair}"
        ) from exc

    raw = wavefield.data()
    denoised = result.final_wavefield.data()
    history = np.asarray(result.example_history, dtype=float)
    if history.shape != (result.iterations, wavefield.n_samples):
        raise ValueError("result example_history has inconsistent shape")
    raw_filtered = _gaussian_narrowband(raw, delta=wavefield.delta, periods=periods)
    denoised_filtered = _gaussian_narrowband_allowing_missing_rows(
        denoised,
        delta=wavefield.delta,
        periods=periods,
    )
    history_filtered = _gaussian_narrowband_allowing_missing_rows(
        history,
        delta=wavefield.delta,
        periods=periods,
    )
    plotted_distances = _jittered_distances(
        wavefield,
        enabled=jitter_duplicate_distances,
    )
    order = np.argsort(plotted_distances, kind="stable")
    time = wavefield.time
    time_grid, distance_grid = np.meshgrid(time, plotted_distances[order])
    fig, axes = plt.subplots(
        periods.size,
        3,
        figsize=(
            _DENOISED_RESULT_FIGURE_WIDTH,
            _DENOISED_RESULT_ROW_HEIGHT * periods.size,
        ),
        sharex=True,
        squeeze=False,
    )

    for row, period in enumerate(periods):
        pair_axis, raw_axis, denoised_axis = axes[row]
        raw_section = _normalized_rows(raw_filtered[:, :, row])
        denoised_section = _normalized_rows(denoised_filtered[:, :, row])
        history_section = _normalized_rows(history_filtered[:, :, row])
        distance = plotted_distances[example_index]
        pair_raw = raw_section[example_index]
        pair_final = denoised_section[example_index]

        for axis, section in ((raw_axis, raw_section), (denoised_axis, denoised_section)):
            axis.pcolormesh(
                time_grid,
                distance_grid,
                section[order],
                cmap="coolwarm",
                shading="auto",
                rasterized=True,
                vmin=-0.8,
                vmax=0.8,
            )
            axis.plot(time, section[example_index] * 0.1 + distance, color="k", linewidth=1.1)

        pair_axis.plot(time, pair_raw - 2.0, color="k", linewidth=1.1)
        pair_axis.text(0.0, -2.0, "ANC", ha="right", va="bottom", fontsize=12)
        for iteration, waveform in enumerate(history_section):
            offset = 2.0 * iteration
            pair_axis.plot(time, waveform + offset, color="r", linewidth=1.0)
            pair_axis.text(
                0.0,
                offset,
                rf"$C^{{{iteration + 3}}}$",
                ha="right",
                va="bottom",
                fontsize=12,
            )

        envelope = np.abs(hilbert(pair_final))
        peak_time = (
            float(time[int(np.argmax(envelope))])
            if np.max(envelope) > EPSILON
            else float(0.5 * (display_tmin + display_tmax))
        )
        for boundary in (peak_time - 2.0 * period, peak_time + 2.0 * period):
            pair_axis.axvline(boundary, color="k", linestyle="--", linewidth=1.5)

        raw_frequency, raw_spectrum, _ = _mean_normalized_spectrum(raw_section, wavefield.delta)
        frequency, denoised_spectrum, peak_index = _mean_normalized_spectrum(
            denoised_section,
            wavefield.delta,
        )
        inset = denoised_axis.inset_axes(_SPECTRUM_INSET_BOUNDS)
        inset.plot(raw_frequency[1:], raw_spectrum[1:], color="k", linewidth=1.0)
        inset.plot(frequency[1:], denoised_spectrum[1:], color="r", linewidth=1.0)
        if peak_index is not None:
            peak_frequency = float(frequency[peak_index])
            inset.axvline(peak_frequency, color="r", linestyle="--", linewidth=1.0)
            spectrum_limits = _centered_spectrum_limits(frequency, peak_frequency)
            if spectrum_limits is not None:
                inset.set_xlim(*spectrum_limits)
            denoised_axis.text(
                0.98,
                0.98,
                f"{peak_frequency:.2f} Hz",
                transform=denoised_axis.transAxes,
                ha="right",
                va="top",
                color="white",
                fontsize=11,
                bbox=dict(facecolor="black", edgecolor="none", alpha=0.55, pad=1.5),
            )
        inset.set(xlabel="Frequency (Hz)", ylabel="Normalized amplitude")
        inset.xaxis.set_label_position("top")
        inset.xaxis.tick_top()
        inset.tick_params(
            axis="x",
            top=True,
            labeltop=True,
            bottom=False,
            labelbottom=False,
            labelsize=8,
        )
        inset.tick_params(axis="y", labelsize=8)

        for column, axis in enumerate(axes[row]):
            axis.set_xlim(display_tmin, display_tmax)
            axis.tick_params(labelsize=12)
            axis.text(
                display_tmin + 0.1,
                0.98,
                _panel_label(3 * row + column),
                transform=axis.get_xaxis_transform(),
                ha="left",
                va="top",
                fontsize=14,
            )
        raw_axis.set_yticks([])
        denoised_axis.yaxis.tick_right()
        denoised_axis.yaxis.set_label_position("right")
        denoised_axis.set_ylabel("Interstation distance (km)", fontsize=14)
        pair_axis.set_yticks([])
        first, second = result.example_pair
        pair_axis.set_title(f"Pair {first}-{second}", fontsize=14)
        raw_axis.set_title(f"ANC at {period:.2f} s", fontsize=14)
        denoised_axis.set_title(
            rf"$C^{{{result.iterations + 2}}}$ at {period:.2f} s",
            fontsize=14,
        )

    for axis in axes[-1]:
        axis.set_xlabel("Correlation time (s)", fontsize=14)
    fig.subplots_adjust(
        left=0.08,
        right=0.92,
        bottom=0.08,
        top=0.94,
        wspace=0.08,
        hspace=0.15,
    )
    return fig, axes


__all__ = ["denoise_station_pair_demo", "plot_denoised_result"]
