# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Phase-shift multi-channel analysis of surface waves (MASW).

The :class:`MASW` state object computes, plots, and persists a normalized
frequency--phase-velocity image from a single surface-wave wavefield.  The
functional wrapper is retained for callers that only need the three arrays.

Notes
-----
Frequency is Hz, phase velocity is km/s, distance is km, and image amplitude
has shape ``(n_frequencies, n_velocities)``.  Computation may use processes
when ``n_jobs`` exceeds one; the input wavefield remains unchanged.

References
----------
Park, C. B., Miller, R. D., and Xia, J. (1998). Imaging dispersion curves of
surface waves on multi-channel record. *SEG Expanded Abstracts*, 1377--1380.
https://doi.org/10.1190/1.1820161
"""

from __future__ import annotations

import os
import tempfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from obspy import Stream, Trace
from scipy.fft import rfft, rfftfreq

from ._validation import (
    trace_start_time,
    validate_frequency_band,
    validate_n_jobs,
)
from .wavefield import Wavefield


DEFAULT_VELOCITIES = np.linspace(0.2, 2.5, 231)
SPECTRAL_EPSILON = np.finfo(float).eps
_MASW_FORMAT = "surface-wave-masw"
_MASW_FORMAT_VERSION = 1


def _validate_velocities(velocities) -> np.ndarray:
    """Return an independent strictly increasing positive velocity grid in km/s."""
    values = np.asarray(velocities, dtype=float)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("velocities must be a non-empty one-dimensional array")
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("velocities must contain finite positive values")
    if np.any(np.diff(values) <= 0):
        raise ValueError("velocities must be strictly increasing")
    return values.copy()


def _compute_frequency_row(
    task: tuple[float, np.ndarray],
    *,
    distance_over_velocity: np.ndarray,
) -> np.ndarray:
    """Compute one normalized phase-shift MASW row."""
    current_frequency, phase_spectrum = task
    steering = np.exp(
        1j * 2.0 * np.pi * current_frequency * distance_over_velocity
    )
    stack = np.abs(np.sum(steering * phase_spectrum, axis=0))
    peak = float(np.max(stack))
    return stack / peak if peak > SPECTRAL_EPSILON else np.zeros_like(stack)


class MASW:
    """Phase-shift multi-channel analysis of one wavefield component."""

    def __init__(
        self,
        wavefield: Wavefield,
        velocities=None,
        fmin: float = 0.5,
        fmax: float = 5.0,
        padding_factor: int = 5,
        dist_threshold: float = 0.2,
    ) -> None:
        """Configure a deferred MASW computation.

        Parameters are phase-velocity samples in km/s, frequency limits in Hz,
        FFT zero-padding factor, and the minimum pair distance in km.  Results
        remain ``None`` until :meth:`compute` is called.
        """
        if not isinstance(wavefield, Wavefield):
            raise TypeError("wavefield must be a Wavefield")
        self.wavefield = wavefield
        self.velocities = _validate_velocities(
            DEFAULT_VELOCITIES if velocities is None else velocities
        )
        self.fmin = fmin
        self.fmax = fmax
        self.padding_factor = padding_factor
        self.dist_threshold = dist_threshold
        self.velocity: np.ndarray | None = None
        self.frequency: np.ndarray | None = None
        self.amplitude: np.ndarray | None = None

    def compute(self, n_jobs: int = 1) -> "MASW":
        """Compute and store the MASW image, then return ``self``."""
        n_jobs = validate_n_jobs(n_jobs)
        stream = self.wavefield.stream()
        fmin, fmax = validate_frequency_band(
            self.fmin,
            self.fmax,
            self.wavefield.sampling_rate,
        )

        if isinstance(self.padding_factor, bool) or not isinstance(
            self.padding_factor,
            (int, np.integer),
        ):
            raise TypeError("padding_factor must be an integer")
        if self.padding_factor < 1:
            raise ValueError("padding_factor must be at least 1")

        data = np.stack([np.asarray(trace.data, dtype=float) for trace in stream])
        distance = self.wavefield.distances
        mask = distance >= self.dist_threshold

        distance = distance[mask]

        data = data[mask]

        nfft = int(self.padding_factor) * data.shape[1]
        all_frequency = rfftfreq(nfft, d=stream[0].stats.delta)
        in_band = (all_frequency >= fmin) & (all_frequency <= fmax)
        frequency = all_frequency[in_band]
        if frequency.size == 0:
            raise ValueError("frequency band contains no FFT bins")

        spectrum = rfft(data, n=nfft, axis=-1)[:, in_band]
        magnitude = np.abs(spectrum)
        phase_spectrum = np.zeros_like(spectrum)
        np.divide(
            spectrum,
            magnitude,
            out=phase_spectrum,
            where=magnitude > SPECTRAL_EPSILON,
        )

        distance_over_velocity = distance[:, None] / self.velocities[None, :]
        compute_row = partial(
            _compute_frequency_row,
            distance_over_velocity=distance_over_velocity,
        )
        tasks = [
            (float(current_frequency), phase_spectrum[:, index, None])
            for index, current_frequency in enumerate(frequency)
        ]
        if n_jobs == 1:
            rows = [compute_row(task) for task in tasks]
        else:
            with Pool(processes=n_jobs) as pool:
                rows = pool.map(compute_row, tasks)
        amplitude = np.asarray(rows, dtype=float)

        self.velocity = self.velocities.copy()
        self.frequency = frequency
        self.amplitude = amplitude
        return self

    def print(self, label: str = "MASW") -> None:
        """Print the computed grid and its strongest normalized sample."""
        label = str(label).strip()
        if not label:
            raise ValueError("label must be a non-empty string")
        if self.frequency is None or self.velocity is None or self.amplitude is None:
            raise RuntimeError("MASW must be computed before printing its summary")

        peak_index = np.unravel_index(
            np.nanargmax(self.amplitude),
            self.amplitude.shape,
        )
        print(
            f"[{label}] grid={self.amplitude.shape}; "
            f"frequency={self.frequency[0]:.3f}.."
            f"{self.frequency[-1]:.3f} Hz; "
            f"velocity={self.velocity[0]:.3f}.."
            f"{self.velocity[-1]:.3f} km/s"
        )
        print(
            "  strongest normalized energy: "
            f"f={self.frequency[peak_index[0]]:.3f} Hz, "
            f"v={self.velocity[peak_index[1]]:.3f} km/s"
        )

    def save(self, path: str | Path, *, overwrite: bool = False) -> Path:
        """Save the MASW configuration, wavefield, and optional results to NPZ."""
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean")
        output_path = Path(path)
        if not output_path.parent.is_dir():
            raise FileNotFoundError(
                f"parent directory does not exist: {output_path.parent}"
            )
        if output_path.exists() and not overwrite:
            raise FileExistsError(f"MASW file already exists: {output_path}")

        computed = all(
            value is not None
            for value in (self.velocity, self.frequency, self.amplitude)
        )
        payload = {
            "format_name": np.asarray(_MASW_FORMAT),
            "format_version": np.asarray(_MASW_FORMAT_VERSION, dtype=np.int64),
            "velocities": self.velocities,
            "fmin": np.asarray(self.fmin, dtype=float),
            "fmax": np.asarray(self.fmax, dtype=float),
            "padding_factor": np.asarray(self.padding_factor, dtype=np.int64),
            "dist_threshold": np.asarray(self.dist_threshold, dtype=float),
            "computed": np.asarray(computed, dtype=np.bool_),
            "frequency": (
                self.frequency if computed else np.empty(0, dtype=float)
            ),
            "amplitude": (
                self.amplitude if computed else np.empty((0, 0), dtype=float)
            ),
            "wavefield_data": self.wavefield.data(),
            "wavefield_sources": np.asarray(
                self.wavefield.sources,
                dtype=np.str_,
            ),
            "wavefield_receivers": np.asarray(
                self.wavefield.receivers,
                dtype=np.str_,
            ),
            "wavefield_distances": self.wavefield.distances,
            "wavefield_delta": np.asarray(self.wavefield.delta, dtype=float),
            "wavefield_start_time": np.asarray(
                trace_start_time(self.wavefield.stream(copy=False)[0]),
                dtype=float,
            ),
            "wavefield_component_present": np.asarray(
                self.wavefield.component is not None,
                dtype=np.bool_,
            ),
            "wavefield_component": np.asarray(
                self.wavefield.component or "",
                dtype=np.str_,
            ),
            "wavefield_check_distance_order": np.asarray(
                self.wavefield.check_distance_order,
                dtype=np.bool_,
            ),
        }

        descriptor, temporary_name = tempfile.mkstemp(
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
        )
        os.close(descriptor)
        temporary_path = Path(temporary_name)
        try:
            with temporary_path.open("wb") as file:
                np.savez_compressed(file, **payload)
            if overwrite:
                os.replace(temporary_path, output_path)
            else:
                os.link(temporary_path, output_path)
                temporary_path.unlink()
        finally:
            if temporary_path.exists():
                temporary_path.unlink()
        return output_path

    @classmethod
    def load(cls, path: str | Path) -> "MASW":
        """Load a MASW object saved by :meth:`save`."""
        input_path = Path(path)
        if not input_path.exists():
            raise FileNotFoundError(f"MASW file does not exist: {input_path}")

        required = {
            "format_name",
            "format_version",
            "velocities",
            "fmin",
            "fmax",
            "padding_factor",
            "dist_threshold",
            "computed",
            "frequency",
            "amplitude",
            "wavefield_data",
            "wavefield_sources",
            "wavefield_receivers",
            "wavefield_distances",
            "wavefield_delta",
            "wavefield_start_time",
            "wavefield_component_present",
            "wavefield_component",
            "wavefield_check_distance_order",
        }
        try:
            with np.load(input_path, allow_pickle=False) as archive:
                missing = sorted(required - set(archive.files))
                if missing:
                    raise ValueError(f"missing required fields: {missing}")

                def scalar(name):
                    """Extract one schema scalar from the currently open NPZ archive."""
                    value = np.asarray(archive[name])
                    if value.shape != ():
                        raise ValueError(f"{name} must be a scalar")
                    return value.item()

                if scalar("format_name") != _MASW_FORMAT:
                    raise ValueError("unsupported MASW format")
                version = scalar("format_version")
                if version != _MASW_FORMAT_VERSION:
                    raise ValueError(f"unsupported MASW format version {version}")

                data = np.asarray(archive["wavefield_data"]).copy()
                sources = np.asarray(archive["wavefield_sources"]).astype(str)
                receivers = np.asarray(archive["wavefield_receivers"]).astype(str)
                distances = np.asarray(
                    archive["wavefield_distances"],
                    dtype=float,
                )
                if data.ndim != 2:
                    raise ValueError("wavefield_data must be two-dimensional")
                n_pairs = data.shape[0]
                if not (
                    sources.shape == receivers.shape == distances.shape == (n_pairs,)
                ):
                    raise ValueError("wavefield metadata does not match its data")

                delta = float(scalar("wavefield_delta"))
                start_time = float(scalar("wavefield_start_time"))
                stream = Stream()
                for index in range(n_pairs):
                    trace = Trace(data=data[index].copy())
                    trace.stats.delta = delta
                    trace.stats.source = str(sources[index])
                    trace.stats.station = str(receivers[index])
                    trace.stats.sac = {
                        "b": start_time,
                        "dist": float(distances[index]),
                        "kevnm": str(sources[index]),
                        "kstnm": str(receivers[index]),
                    }
                    stream.append(trace)

                component_present = scalar("wavefield_component_present")
                component = (
                    str(scalar("wavefield_component"))
                    if component_present
                    else None
                )
                wavefield = Wavefield(
                    stream,
                    component=component,
                    copy=False,
                    check_distance_order=bool(
                        scalar("wavefield_check_distance_order")
                    ),
                )
                masw = cls(
                    wavefield,
                    velocities=np.asarray(archive["velocities"], dtype=float),
                    fmin=float(scalar("fmin")),
                    fmax=float(scalar("fmax")),
                    padding_factor=int(scalar("padding_factor")),
                    dist_threshold=float(scalar("dist_threshold")),
                )
                if bool(scalar("computed")):
                    frequency = np.asarray(archive["frequency"], dtype=float)
                    amplitude = np.asarray(archive["amplitude"], dtype=float)
                    if frequency.ndim != 1 or amplitude.shape != (
                        frequency.size,
                        masw.velocities.size,
                    ):
                        raise ValueError("MASW result arrays have inconsistent shapes")
                    masw.velocity = masw.velocities.copy()
                    masw.frequency = frequency.copy()
                    masw.amplitude = amplitude.copy()
                return masw
        except FileNotFoundError:
            raise
        except Exception as exc:
            raise ValueError(f"invalid MASW file {input_path}: {exc}") from exc

    def plot(self, ax=None, *, cmap: str = "coolwarm"):
        """Plot the computed frequency-phase-velocity image.

        Returns ``(fig, ax)`` and does not call ``show``. Call ``compute``
        before plotting.
        """
        import matplotlib.pyplot as plt

        if self.velocity is None or self.frequency is None or self.amplitude is None:
            raise ValueError("MASW results are not available; call compute() first")

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        frequency_grid, velocity_grid = np.meshgrid(
            self.frequency,
            self.velocity,
        )
        image = ax.pcolormesh(
            frequency_grid,
            velocity_grid,
            self.amplitude.T,
            shading="auto",
            cmap=cmap,
            vmin=0,
            vmax=1,
            rasterized=True,
        )
        ax.set(
            xlabel="Frequency (Hz)",
            ylabel="Phase velocity (km/s)",
            title=(
                f"MASW ({self.wavefield.component})"
                if self.wavefield.component
                else "MASW"
            ),
        )
        #fig.colorbar(image, ax=ax, label="Normalized amplitude")
        return fig, ax


def compute_masw(
    wavefield: Wavefield,
    velomin: float = 0.2,
    velomax: float = 2.5,
    fmin: float = 0.5,
    fmax: float = 5.0,
    padding_factor: int = 5,
    *,
    velocities=None,
    dist_threshold: float = 0.2,
    n_jobs: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compatibility wrapper around :class:`MASW`.

    ``velomin``/``velomax`` define the default linear velocity grid. A custom
    grid can be supplied with the keyword-only ``velocities`` argument.
    """
    velocity_grid = (
        np.linspace(velomin, velomax, 231)
        if velocities is None
        else velocities
    )
    masw = MASW(
        wavefield,
        velocities=velocity_grid,
        fmin=fmin,
        fmax=fmax,
        padding_factor=padding_factor,
        dist_threshold=dist_threshold,
    ).compute(n_jobs=n_jobs)
    return masw.velocity, masw.frequency, masw.amplitude


__all__ = ["MASW", "compute_masw"]
