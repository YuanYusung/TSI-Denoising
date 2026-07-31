# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Phase-matched filtering for dispersive surface-wave mode separation.

An input correlation is phase-shifted using a frequency--phase-velocity
reference curve, isolated around zero correlation time with a Gaussian window,
and shifted back.  A curve may be supplied directly or picked interactively
from a computed or cached MASW image.

Notes
-----
Reference curves have columns ``(frequency_Hz, phase_velocity_km_s)``.  The
function returns new wavefield objects and does not mutate its input.  The GUI
picker requires an interactive Matplotlib backend.

References
----------
Levshin, A. and Ritzwoller, M. (2001). Automated detection, extraction, and
measurement of regional surface waves. *Pure and Applied Geophysics*, 158,
1531--1545. https://doi.org/10.1007/PL00001233
"""

from __future__ import annotations

import re
import warnings
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from .._validation import validate_frequency_band
from ..masw import MASW
from ..preprocessing import _apply_velocity_taper
from ..wavefield import Wavefield


DEFAULT_T_WINDOW = 0.2
PEAK_EPSILON = 1e-12
_CACHE_COMPONENT_PATTERN = re.compile(r"[^a-z0-9._-]+")


def _coerce_reference_curve(reference_curve) -> np.ndarray:
    """Validate and return a sorted ``(frequency, velocity)`` curve."""
    if isinstance(reference_curve, tuple) and len(reference_curve) == 2:
        frequencies = np.asarray(reference_curve[0], dtype=float)
        velocities = np.asarray(reference_curve[1], dtype=float)
        if frequencies.ndim != 1 or velocities.ndim != 1:
            raise ValueError(
                "reference curve frequency and velocity arrays must be one-dimensional"
            )
        if frequencies.size != velocities.size:
            raise ValueError(
                "reference curve frequency and velocity arrays must have equal length"
            )
        curve = np.column_stack((frequencies, velocities))
    else:
        try:
            curve = np.asarray(reference_curve, dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "reference_curve must be an Nx2 array or a (frequency, velocity) tuple"
            ) from exc

        if curve.ndim != 2 or curve.shape[1] != 2:
            raise ValueError(
                "reference_curve must be an Nx2 array of frequency and velocity"
            )

    if curve.shape[0] < 2:
        raise ValueError("reference curve must contain at least two points")
    if not np.all(np.isfinite(curve)):
        raise ValueError("reference curve must contain only finite values")
    if np.any(curve[:, 0] <= 0) or np.any(curve[:, 1] <= 0):
        raise ValueError("reference curve frequencies and velocities must be positive")

    curve = curve[np.argsort(curve[:, 0], kind="mergesort")]
    if np.any(np.diff(curve[:, 0]) <= 0):
        raise ValueError("reference curve frequencies must be strictly increasing")
    return curve


def _peak_normalize(data: np.ndarray, *, epsilon: float = PEAK_EPSILON) -> np.ndarray:
    """Return a finite waveform normalized by peak amplitude or all zeros."""
    values = np.asarray(data, dtype=float)
    peak = float(np.max(np.abs(values)))
    if not np.isfinite(peak):
        raise ValueError("waveform amplitude must be finite")
    if peak <= epsilon:
        return np.zeros_like(values)
    return values / peak


def _phase_shift(
    data: np.ndarray,
    *,
    delta: float,
    distance: float,
    slowness: Callable[[np.ndarray], np.ndarray],
    fmin: float,
    fmax: float,
    reverse: bool = False,
) -> np.ndarray:
    """Apply one direction of a reference-curve phase shift in the FFT domain.

    ``distance`` is km, ``delta`` is s, and the callable returns slowness in
    s/km for the supplied Hz frequencies.
    """
    normalized = _peak_normalize(data)
    spectrum = np.fft.rfft(normalized)
    frequencies = np.fft.rfftfreq(normalized.size, d=delta)
    in_band = (frequencies >= fmin) & (frequencies <= fmax)

    phase = np.ones(spectrum.shape, dtype=complex)
    if np.any(in_band):
        phase_slowness = np.asarray(slowness(frequencies[in_band]), dtype=float)
        if not np.all(np.isfinite(phase_slowness)) or np.any(phase_slowness <= 0):
            raise ValueError("reference curve gives invalid phase slowness in the frequency band")
        delays = distance * phase_slowness
        sign = -1.0 if reverse else 1.0
        phase[in_band] = np.exp(
            sign * 1j * 2.0 * np.pi * frequencies[in_band] * delays
        )

    shifted = np.fft.irfft(spectrum * phase, n=normalized.size)
    return _peak_normalize(shifted)


def _gaussian_taper(
    times: np.ndarray,
    data: np.ndarray,
    *,
    t_window: float,
    keep_positive: bool = True,
) -> np.ndarray:
    """Taper the unwanted correlation-time branch with a Gaussian in seconds."""
    sigma = t_window / 3.0
    window = np.exp(-0.5 * (times / sigma) ** 2)
    if keep_positive:
        window[times > 0] = 1.0
    else:
        window[times < 0] = 1.0
    return _peak_normalize(np.asarray(data, dtype=float) * window)


def _pick_dispersion_curve(masw: MASW) -> np.ndarray:
    """Display a computed MASW image and collect manual frequency-velocity picks."""
    import matplotlib.pyplot as plt

    try:
        fig, ax = masw.plot()
        ax.set_title(
            "Pick dispersion curve: left-click points, press Enter to finish"
        )
        points = fig.ginput(n=-1, timeout=-1, show_clicks=True)
    except Exception as exc:
        raise RuntimeError(
            "unable to open the dispersion-picking GUI; provide reference_curve explicitly"
        ) from exc
    finally:
        if "fig" in locals():
            try:
                plt.close(fig)
            except TypeError:
                # Keep cleanup tolerant of lightweight Figure doubles used by
                # callers' tests; real Matplotlib figures close normally.
                pass

    if points is None or len(points) < 2:
        raise ValueError(
            "dispersion picking was cancelled or produced fewer than two points"
        )
    curve = np.asarray(points, dtype=float)
    curve = curve[
        np.isfinite(curve).all(axis=1)
        & (curve[:, 0] > 0)
        & (curve[:, 1] > 0)
    ]
    if curve.shape[0] < 2:
        raise ValueError("dispersion picking requires at least two positive points")
    return _coerce_reference_curve(curve)


def _component_cache_filename(component: str | None) -> str | None:
    """Return the safe automatic MASW cache name for a component."""
    if component is None:
        return None
    normalized = _CACHE_COMPONENT_PATTERN.sub(
        "_",
        component.strip().lower(),
    ).strip("._")
    return f"{normalized}_masw.npz" if normalized else None


def _require_computed_masw(masw: MASW, path: Path) -> MASW:
    """Require a loaded cache to contain a computed MASW image."""
    if any(
        value is None
        for value in (masw.velocity, masw.frequency, masw.amplitude)
    ):
        raise ValueError(
            f"MASW cache does not contain computed results: {path}"
        )
    return masw


def _load_explicit_masw_cache(path) -> MASW:
    """Load a required, explicitly named MASW cache."""
    try:
        cache_path = Path(path)
    except TypeError as exc:
        raise TypeError("masw_cache_path must be a path-like value") from exc
    if cache_path.suffix.lower() != ".npz":
        raise ValueError("masw_cache_path must identify an NPZ file")
    return _require_computed_masw(MASW.load(cache_path), cache_path)


def _masw_for_picking(
    wavefield: Wavefield,
    *,
    fmin: float,
    fmax: float,
    masw_cache_path=None,
) -> MASW:
    """Load a usable picking background or compute one from *wavefield*."""
    if masw_cache_path is not None:
        return _load_explicit_masw_cache(masw_cache_path)

    filename = _component_cache_filename(wavefield.component)
    if filename is not None:
        working_directory = Path.cwd()
        candidates = (
            working_directory / filename,
            working_directory / "processed" / filename,
            working_directory / "cache" / filename,
        )
        for candidate in candidates:
            if not candidate.is_file():
                continue
            try:
                return _require_computed_masw(
                    MASW.load(candidate),
                    candidate,
                )
            except (OSError, ValueError) as exc:
                warnings.warn(
                    f"ignoring unusable automatic MASW cache "
                    f"{candidate}: {exc}",
                    UserWarning,
                    stacklevel=2,
                )

    return MASW(wavefield, fmin=fmin, fmax=fmax).compute()


def phase_match_separate(
    wavefield: Wavefield,
    reference_curve=None,
    *,
    fmin: float = 0.5,
    fmax: float = 5.0,
    t_window: float = DEFAULT_T_WINDOW,
    keep_positive: bool = True,
    return_reference: bool = False,
    masw_cache_path: str | Path | None = None,
    vmin: float = 0.1,
    vmax: float = 2.5,
    taper_fraction: float = 0.05,
) -> Wavefield | tuple[Wavefield, np.ndarray]:
    """Separate a mode using phase-matched filtering.

    The reference curve is an ``(N, 2)`` array containing frequency in Hz and
    phase velocity in km/s.  A two-array ``(frequency, velocity)`` tuple is
    also accepted. If no curve is supplied, a computed MASW cache is loaded
    when available and displayed for manual picking; otherwise MASW is
    computed from ``wavefield``. Phase matching is followed by the same
    distance-dependent velocity-window taper used by
    :func:`preprocess_stream`; apply band-pass filtering separately when it
    is required.

    Parameters
    ----------
    wavefield:
        One input wavefield.
    reference_curve:
        Optional reference dispersion curve.
    fmin, fmax:
        Frequency band in which the reference phase is applied.
    t_window:
        Gaussian-window width in seconds around zero correlation time.
    keep_positive:
        If true, retain the positive-time half-axis and taper the negative
        half-axis.  If false, retain the negative-time half-axis instead.
    return_reference:
        If true, return ``(separated_wavefield, used_curve)``.
    masw_cache_path:
        Optional path to a required computed MASW NPZ cache used as the
        picking background. Without an explicit path, caches named
        ``<component>_masw.npz`` are checked in the current directory,
        ``processed/``, and ``cache/`` before MASW is recomputed.
    vmin, vmax, taper_fraction:
        Parameters for the same distance-dependent velocity-window taper used
        by :func:`preprocess_stream`.  The taper is applied after separation.
    """
    if not isinstance(wavefield, Wavefield):
        raise TypeError("wavefield must be a Wavefield")
    if reference_curve is not None and masw_cache_path is not None:
        raise ValueError(
            "masw_cache_path cannot be used when reference_curve is provided"
        )
    fmin, fmax = validate_frequency_band(fmin, fmax, wavefield.sampling_rate)
    t_window = float(t_window)
    if not np.isfinite(t_window) or t_window <= 0:
        raise ValueError("t_window must be finite and positive")
    if not isinstance(keep_positive, (bool, np.bool_)):
        raise TypeError("keep_positive must be a boolean")

    if reference_curve is None:
        masw = _masw_for_picking(
            wavefield,
            fmin=fmin,
            fmax=fmax,
            masw_cache_path=masw_cache_path,
        )
        curve = _pick_dispersion_curve(masw)
    else:
        curve = _coerce_reference_curve(reference_curve)

    frequencies = curve[:, 0]
    velocities = curve[:, 1]
    velocity_interp = interp1d(
        frequencies,
        velocities,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
    )

    def slowness(query_frequencies):
        """Interpolate phase velocity (km/s) then convert it to slowness (s/km)."""
        # Interpolate phase velocity first, then convert to phase slowness.
        # This avoids the nonphysical negative values that can result from
        # linearly extrapolating 1 / velocity directly.
        phase_velocity = np.asarray(
            velocity_interp(query_frequencies), dtype=float
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            return 1.0 / phase_velocity

    times = wavefield.time
    output_stream = wavefield.stream().copy()
    for output_trace, input_trace, distance in zip(
        output_stream, wavefield.stream(), wavefield.distances
    ):
        zero_phase = _phase_shift(
            input_trace.data,
            delta=wavefield.delta,
            distance=float(distance),
            slowness=slowness,
            fmin=fmin,
            fmax=fmax,
        )
        windowed = _gaussian_taper(
            times,
            zero_phase,
            t_window=t_window,
            keep_positive=bool(keep_positive),
        )
        output_trace.data = _phase_shift(
            windowed,
            delta=wavefield.delta,
            distance=float(distance),
            slowness=slowness,
            fmin=fmin,
            fmax=fmax,
            reverse=True,
        )

    separated = Wavefield(
        output_stream,
        component="phase_matched",
        copy=False,
        check_distance_order=wavefield.check_distance_order,
    )
    separated = _apply_velocity_taper(
        separated,
        vmin=vmin,
        vmax=vmax,
        taper_fraction=taper_fraction,
    )
    if return_reference:
        return separated, curve.copy()
    return separated


__all__ = ["phase_match_separate"]
