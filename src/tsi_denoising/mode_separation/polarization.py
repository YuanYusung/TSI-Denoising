# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Polarization separation for four-component surface-wave wavefields.

Aligned ``ZZ``, ``RZ``, ``ZR``, and ``RR`` cross-correlations are combined
with Hilbert-transformed cross components to emphasize retrograde or prograde
Rayleigh-wave particle motion.  The input wavefields are never modified.

Notes
-----
All inputs must share station pairs, samples, sampling interval, time axis,
and pair distances.  Velocity-taper parameters use km/s and output values are
new :class:`~tsi_denoising.Wavefield` instances.

References
----------
Gribler, G. and Mikesell, T. D. (2019). Methods to isolate retrograde and
prograde Rayleigh-wave signals. *Geophysical Journal International*, 219,
975--994. https://doi.org/10.1093/gji/ggz341
"""

from __future__ import annotations

import numpy as np
from obspy import Stream
from scipy.signal import hilbert

from ..preprocessing import _apply_velocity_taper
from ..wavefield import Wavefield


RMS_EPSILON = 1e-12


def _rms_normalize(data: np.ndarray, *, epsilon: float = RMS_EPSILON) -> np.ndarray:
    """Return a copy of *data* normalized to unit RMS.

    A zero-amplitude input is represented by zeros rather than by NaNs.  This
    is useful for incomplete or deliberately masked component wavefields.
    """
    values = np.asarray(data, dtype=float)
    rms = float(np.sqrt(np.mean(values**2)))
    if not np.isfinite(rms):
        raise ValueError("waveform RMS must be finite")
    if rms <= epsilon:
        return np.zeros_like(values)
    return values / rms


def _validate_component_wavefields(
    zz: Wavefield,
    rz: Wavefield,
    zr: Wavefield,
    rr: Wavefield,
) -> dict[str, tuple[int, ...]]:
    """Validate four aligned components and map each one to ``ZZ`` pair order.

    Matching pair sets, sampling, correlation times, and distances are
    required before their arrays can be combined safely.
    """
    components = {"ZZ": zz, "RZ": rz, "ZR": zr, "RR": rr}
    for name, wavefield in components.items():
        if not isinstance(wavefield, Wavefield):
            raise TypeError(f"{name} must be a Wavefield")

    reference_pairs = zz.pairs
    reference_pair_set = set(reference_pairs)
    if len(reference_pair_set) != len(reference_pairs):
        raise ValueError("ZZ contains duplicate station pairs")

    indices: dict[str, tuple[int, ...]] = {"ZZ": tuple(range(zz.n_pairs))}
    for name, wavefield in components.items():
        if name == "ZZ":
            continue
        pairs = wavefield.pairs
        if set(pairs) != reference_pair_set:
            raise ValueError(
                f"{name} station pairs do not match ZZ; "
                "all four components must contain the same pairs"
            )
        pair_to_index = {pair: index for index, pair in enumerate(pairs)}
        indices[name] = tuple(pair_to_index[pair] for pair in reference_pairs)

    for name, wavefield in components.items():
        if wavefield.n_samples != zz.n_samples:
            raise ValueError(
                f"{name} has {wavefield.n_samples} samples; "
                f"expected {zz.n_samples}"
            )
        if not np.isclose(wavefield.delta, zz.delta, rtol=0.0, atol=1e-12):
            raise ValueError(
                f"{name} sampling interval {wavefield.delta} does not match ZZ"
            )
        if not np.isclose(
            wavefield.time[0], zz.time[0], rtol=0.0, atol=1e-9
        ):
            raise ValueError(f"{name} start time does not match ZZ")

        distances = wavefield.distances[np.asarray(indices[name])]
        if not np.allclose(distances, zz.distances, rtol=0.0, atol=1e-9):
            raise ValueError(f"{name} distances do not match ZZ")

    return indices


def polarization_separate(
    *,
    zz: Wavefield,
    rz: Wavefield,
    zr: Wavefield,
    rr: Wavefield,
    swap_polarization: bool = False,
    use_rr: bool = True,
    vmin: float = 0.1,
    vmax: float = 2.5,
    taper_fraction: float = 0.05,
) -> dict[str, Wavefield]:
    """Separate four-component wavefields into retrograde and prograde parts.

    Each station-pair component is independently RMS-normalized.  The radial-
    vertical cross-components (``RZ`` and ``ZR``) use the imaginary Hilbert
    component before normalization, following the Step 1 formulation used by
    the RR/SJFZ workflows.  The returned mapping contains new ``Wavefield``
    objects and the input wavefields are not modified.

    Parameters
    ----------
    zz, rz, zr, rr:
        Four component wavefields.  The arguments are keyword-only to make
        component order explicit at call sites.
    swap_polarization:
        If true, exchange the retrograde and prograde output definitions.
    use_rr:
        If false, omit the RMS-normalized ``RR`` component from the prograde
        combination.  The remaining three terms are averaged by three.
    vmin, vmax, taper_fraction:
        Parameters for the same distance-dependent velocity-window taper used
        by :func:`preprocess_stream`.  The taper is applied after separation.
    """
    if not isinstance(swap_polarization, (bool, np.bool_)):
        raise TypeError("swap_polarization must be a boolean")
    if not isinstance(use_rr, (bool, np.bool_)):
        raise TypeError("use_rr must be a boolean")

    indices = _validate_component_wavefields(zz, rz, zr, rr)
    component_wavefields = {"ZZ": zz, "RZ": rz, "ZR": zr, "RR": rr}
    normalized: dict[str, list[np.ndarray]] = {
        name: [] for name in component_wavefields
    }

    for pair_index in range(zz.n_pairs):
        for name, wavefield in component_wavefields.items():
            source_index = indices[name][pair_index]
            data = np.asarray(wavefield[source_index].data, dtype=float)
            if name in {"RZ", "ZR"}:
                data = np.imag(hilbert(data))
            normalized[name].append(_rms_normalize(data))

    if use_rr:
        prograde_data = [
            (
                normalized["ZZ"][i]
                + normalized["ZR"][i]
                - normalized["RZ"][i]
                + normalized["RR"][i]
            )
            / 4.0
            for i in range(zz.n_pairs)
        ]
        retrograde_data = [
            (
                normalized["ZZ"][i]
                - normalized["ZR"][i]
                + normalized["RZ"][i]
                + normalized["RR"][i]
            )
            / 4.0
            for i in range(zz.n_pairs)
        ]
    else:
        prograde_data = [
            (
                2 * normalized["ZZ"][i]
                + normalized["ZR"][i]
                - normalized["RZ"][i]
            )
            / 4.0
            for i in range(zz.n_pairs)
        ]
        retrograde_data = [
            (
                2 * normalized["ZZ"][i]
                - normalized["ZR"][i]
                + normalized["RZ"][i]
            )
            / 4.0
            for i in range(zz.n_pairs)
        ]

    if swap_polarization:
        retrograde_data, prograde_data = prograde_data, retrograde_data

    retrograde_stream = Stream()
    prograde_stream = Stream()
    for trace, retrograde, prograde in zip(
        zz.stream(), retrograde_data, prograde_data
    ):
        retrograde_trace = trace.copy()
        retrograde_trace.data = np.asarray(retrograde, dtype=float)
        retrograde_stream.append(retrograde_trace)

        prograde_trace = trace.copy()
        prograde_trace.data = np.asarray(prograde, dtype=float)
        prograde_stream.append(prograde_trace)

    separated = {
        "retrograde": Wavefield(
            retrograde_stream,
            component="retrograde",
            copy=False,
            check_distance_order=zz.check_distance_order,
        ),
        "prograde": Wavefield(
            prograde_stream,
            component="prograde",
            copy=False,
            check_distance_order=zz.check_distance_order,
        ),
    }
    return {
        name: _apply_velocity_taper(
            wavefield,
            vmin=vmin,
            vmax=vmax,
            taper_fraction=taper_fraction,
        )
        for name, wavefield in separated.items()
    }


__all__ = ["polarization_separate"]
