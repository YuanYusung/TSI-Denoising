# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Validated data model and persistence utilities for one wavefield component.

:class:`Wavefield` owns a normalized ObsPy stream of pair-aligned correlation
traces together with optional component metadata.  It provides defensive data
access, plotting, NPZ persistence, and an explicitly in-place preprocessing
method.

Notes
-----
Each trace must share sample count, interval, and correlation-time origin.
Pair distance is km, time is s, and waveform matrices have shape
``(n_pairs, n_samples)``.  File operations use a versioned NPZ schema.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path

import numpy as np
from obspy import Stream, Trace

from ._validation import (
    normalize_pair_directions,
    trace_distance,
    trace_pair,
    trace_start_time,
    station_number,
    validate_reference_distance_order,
    validate_station_numbers,
    validate_stream,
)

_WAVEFIELD_FORMAT = "surface-wave-wavefield"
_WAVEFIELD_FORMAT_VERSION = 1
_WAVEFIELD_NPZ_FIELDS = frozenset(
    {
        "format_name",
        "format_version",
        "data",
        "sources",
        "receivers",
        "distances",
        "delta",
        "start_time",
        "component_present",
        "component",
        "check_distance_order",
        "n_pairs",
        "n_samples",
    }
)


def _scalar(array, name: str):
    """Extract an NPZ scalar field and reject non-scalar archive entries."""
    value = np.asarray(array)
    if value.shape != ():
        raise ValueError(f"{name} must be a scalar")
    return value.item()


def _canonical_plot_pair(example_pair) -> tuple[str, str]:
    """Validate and normalize a user-supplied station pair for plotting."""
    if isinstance(example_pair, (str, bytes)):
        raise TypeError(
            "example_pair must be a two-element sequence, not a string"
        )
    try:
        values = tuple(example_pair)
    except TypeError as exc:
        raise TypeError(
            "example_pair must be an iterable containing two station names"
        ) from exc
    if len(values) != 2:
        raise ValueError("example_pair must contain exactly two station names")
    if not all(isinstance(value, str) for value in values):
        raise TypeError("example_pair entries must be strings")

    first, second = (value.strip() for value in values)
    if not first or not second:
        raise ValueError("example_pair must contain non-empty station names")
    if first == second or station_number(first) == station_number(second):
        raise ValueError("example_pair must contain two different stations")
    return (first, second) if station_number(first) < station_number(second) else (second, first)


def _local_distance_spacing(distances: np.ndarray, distance: float) -> float | None:
    """Return the nearest distinct distance spacing, if one exists."""
    distinct_distances = np.unique(distances)
    index = int(np.searchsorted(distinct_distances, distance))
    gaps = []
    if index > 0:
        gaps.append(distance - distinct_distances[index - 1])
    if index + 1 < len(distinct_distances):
        gaps.append(distinct_distances[index + 1] - distance)
    return min(gaps) if gaps else None


def _stable_jitter(pair: tuple[str, str], distance: float, limit: float) -> float:
    """Return a stable pseudo-random offset in ``[-limit, limit]``."""
    token = f"{pair[0]}\0{pair[1]}\0{distance.hex()}".encode()
    value = int.from_bytes(hashlib.sha256(token).digest()[:8], "big")
    fraction = value / ((1 << 64) - 1)
    return (2.0 * fraction - 1.0) * limit


class Wavefield:
    """A validated collection of pair-aligned surface-wave traces.

    One ``Wavefield`` contains one ObsPy ``Stream``. The optional ``component``
    name is metadata only; different components should be stored in separate
    ``Wavefield`` objects.
    """

    def __init__(
        self,
        stream: Stream,
        *,
        component: str | None = None,
        copy: bool = True,
        check_distance_order: bool = True,
    ) -> None:
        """Validate and store a pair-aligned ObsPy stream.

        Parameters
        ----------
        stream : obspy.Stream
            Non-empty correlations with compatible sample geometry.
        component : str or None, optional
            Metadata label such as ``ZZ``; it does not alter processing.
        copy : bool, default=True
            Defensively copy input traces before direction normalization.
        check_distance_order : bool, default=True
            Enforce monotonic reference-station distance ordering.

        Notes
        -----
        Station-pair directions are normalized using numeric station suffixes.
        Input metadata must provide pair names, distance in km, and time origin.
        """
        if not isinstance(stream, Stream):
            raise TypeError("stream must be an obspy.Stream")
        if len(stream) == 0:
            raise ValueError("stream is empty")
        if not isinstance(check_distance_order, bool):
            raise TypeError("check_distance_order must be a boolean")
        if component is not None:
            component = str(component).strip()
            if not component:
                raise ValueError("component must be a non-empty name")

        working_stream = stream.copy() if copy else stream
        raw_pairs = tuple(trace_pair(trace) for trace in working_stream)
        validate_station_numbers(raw_pairs)
        working_stream = normalize_pair_directions(working_stream)
        validate_stream(working_stream)
        if check_distance_order:
            validate_reference_distance_order(working_stream)

        self._stream = working_stream
        self._component = component
        self._check_distance_order = check_distance_order
        self._pairs = tuple(trace_pair(trace) for trace in self._stream)

    @property
    def component(self) -> str | None:
        """Optional component name, such as ``ZZ`` or ``RR``."""
        return self._component

    @property
    def check_distance_order(self) -> bool:
        """Whether linear-array distance ordering is checked on construction."""
        return self._check_distance_order

    @property
    def pairs(self) -> tuple[tuple[str, str], ...]:
        """Source-receiver pairs in stream order."""
        return self._pairs

    @property
    def n_pairs(self) -> int:
        """Number of source-receiver correlations stored in this wavefield."""
        return len(self._pairs)

    @property
    def n_samples(self) -> int:
        """Samples per correlation trace."""
        return int(self._stream[0].stats.npts)

    @property
    def delta(self) -> float:
        """Shared correlation sampling interval in seconds."""
        return float(self._stream[0].stats.delta)

    @property
    def sampling_rate(self) -> float:
        """Shared sampling rate in Hz (the reciprocal of :attr:`delta`)."""
        return float(self._stream[0].stats.sampling_rate)

    @property
    def time(self) -> np.ndarray:
        """Return an independent one-dimensional correlation-time axis in seconds."""
        start = trace_start_time(self._stream[0])
        return start + np.arange(self.n_samples) * self.delta

    @property
    def distances(self) -> np.ndarray:
        """Return source-receiver distances in km, aligned with :attr:`pairs`."""
        return np.asarray(
            [trace_distance(trace) for trace in self._stream],
            dtype=float,
        )

    @property
    def sources(self) -> tuple[str, ...]:
        """Return source station names in stream order."""
        return tuple(pair[0] for pair in self._pairs)

    @property
    def receivers(self) -> tuple[str, ...]:
        """Return receiver station names in stream order."""
        return tuple(pair[1] for pair in self._pairs)

    def stream(self, *, copy: bool = True) -> Stream:
        """Return a copy of the ObsPy ``Stream`` unless ``copy=False``."""
        return self._stream.copy() if copy else self._stream

    def data(self) -> np.ndarray:
        """Return a new ``(n_pairs, n_samples)`` data array."""
        return np.stack([trace.data for trace in self._stream])

    def save(
        self,
        path: str | Path,
        *,
        overwrite: bool = False,
    ) -> Path:
        """Save this wavefield to a versioned, compressed NPZ file."""
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean")

        output_path = Path(path)
        if not output_path.parent.is_dir():
            raise FileNotFoundError(
                f"parent directory does not exist: {output_path.parent}"
            )
        if output_path.exists() and not overwrite:
            raise FileExistsError(f"Wavefield file already exists: {output_path}")

        payload = {
            "format_name": np.asarray(_WAVEFIELD_FORMAT),
            "format_version": np.asarray(
                _WAVEFIELD_FORMAT_VERSION,
                dtype=np.int64,
            ),
            "data": self.data(),
            "sources": np.asarray(self.sources, dtype=np.str_),
            "receivers": np.asarray(self.receivers, dtype=np.str_),
            "distances": self.distances,
            "delta": np.asarray(self.delta, dtype=np.float64),
            "start_time": np.asarray(
                trace_start_time(self._stream[0]),
                dtype=np.float64,
            ),
            "component_present": np.asarray(
                self.component is not None,
                dtype=np.bool_,
            ),
            "component": np.asarray(self.component or "", dtype=np.str_),
            "check_distance_order": np.asarray(
                self.check_distance_order,
                dtype=np.bool_,
            ),
            "n_pairs": np.asarray(self.n_pairs, dtype=np.int64),
            "n_samples": np.asarray(self.n_samples, dtype=np.int64),
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
                file.flush()
                os.fsync(file.fileno())

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
    def load(cls, path: str | Path) -> "Wavefield":
        """Load and validate a wavefield saved by :meth:`save`."""
        input_path = Path(path)
        if not input_path.exists():
            raise FileNotFoundError(f"Wavefield file does not exist: {input_path}")

        try:
            with np.load(input_path, allow_pickle=False) as archive:
                missing = sorted(_WAVEFIELD_NPZ_FIELDS - set(archive.files))
                if missing:
                    raise ValueError(f"missing required fields: {missing}")

                if _scalar(archive["format_name"], "format_name") != _WAVEFIELD_FORMAT:
                    raise ValueError("unsupported Wavefield format")
                version = _scalar(archive["format_version"], "format_version")
                if version != _WAVEFIELD_FORMAT_VERSION:
                    raise ValueError(
                        f"unsupported Wavefield format version {version}"
                    )

                n_pairs = int(_scalar(archive["n_pairs"], "n_pairs"))
                n_samples = int(_scalar(archive["n_samples"], "n_samples"))
                data = np.asarray(archive["data"]).copy()
                sources = np.asarray(archive["sources"]).astype(str)
                receivers = np.asarray(archive["receivers"]).astype(str)
                distances = np.asarray(archive["distances"], dtype=float)
                if data.shape != (n_pairs, n_samples):
                    raise ValueError(
                        "data shape does not match n_pairs and n_samples"
                    )
                if not (
                    sources.shape == receivers.shape == distances.shape == (n_pairs,)
                ):
                    raise ValueError("wavefield metadata does not match its data")

                delta = float(_scalar(archive["delta"], "delta"))
                start_time = float(_scalar(archive["start_time"], "start_time"))
                component = (
                    str(_scalar(archive["component"], "component"))
                    if bool(
                        _scalar(
                            archive["component_present"],
                            "component_present",
                        )
                    )
                    else None
                )
                check_distance_order = bool(
                    _scalar(
                        archive["check_distance_order"],
                        "check_distance_order",
                    )
                )

            stream = Stream()
            for index, (source, receiver, distance) in enumerate(
                zip(sources, receivers, distances)
            ):
                trace = Trace(data=data[index].copy())
                trace.stats.delta = delta
                trace.stats.source = str(source)
                trace.stats.station = str(receiver)
                trace.stats.sac = {
                    "b": start_time,
                    "dist": float(distance),
                    "kevnm": str(source),
                    "kstnm": str(receiver),
                }
                stream.append(trace)

            return cls(
                stream,
                component=component,
                copy=False,
                check_distance_order=check_distance_order,
            )
        except FileNotFoundError:
            raise
        except Exception as exc:
            raise ValueError(
                f"invalid Wavefield file {input_path}: {exc}"
            ) from exc

    def plot(
        self,
        ax=None,
        *,
        tshow_min: float | None = -1.0,
        tshow_max: float | None = 5.0,
        normalize: bool = True,
        cmap: str = "coolwarm",
        jitter_duplicate_distances: bool = False,
        example_pair: tuple[str, str] | None = None,
    ):
        """Plot the wavefield as a distance-time image.

        Traces are sorted by distance and optionally normalized individually.
        Set ``jitter_duplicate_distances`` to separate pairs at identical
        distances with stable, small vertical offsets.  ``example_pair`` may
        name either direction of a stored pair; its normalized trace is
        overlaid as a black wiggle.  Returns ``(fig, ax)`` without calling
        ``show``.
        """
        import matplotlib.pyplot as plt

        if not isinstance(jitter_duplicate_distances, bool):
            raise TypeError("jitter_duplicate_distances must be a boolean")

        data = np.asarray(self.data(), dtype=float)
        if normalize:
            scale = np.max(np.abs(data), axis=1, keepdims=True)
            data = np.divide(
                data,
                scale,
                out=np.zeros_like(data),
                where=scale > 0,
            )

        distances = self.distances
        plotted_distances = distances.copy()
        if jitter_duplicate_distances:
            for distance in np.unique(distances):
                indices = np.flatnonzero(distances == distance)
                if len(indices) < 2:
                    continue
                spacing = _local_distance_spacing(distances, distance)
                if spacing is None:
                    continue
                limit = 0.2 * spacing
                for index in indices:
                    plotted_distances[index] += _stable_jitter(
                        self.pairs[index],
                        float(distance),
                        limit,
                    )

        order = np.argsort(plotted_distances, kind="stable")
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        time_grid, distance_grid = np.meshgrid(
            self.time,
            plotted_distances[order],
        )
        image = ax.pcolormesh(
            time_grid,
            distance_grid,
            data[order],
            shading="auto",
            cmap=cmap,
            vmin=-1 if normalize else None,
            vmax=1 if normalize else None,
            rasterized=True,
        )

        if example_pair is not None:
            pair = _canonical_plot_pair(example_pair)
            try:
                pair_index = self.pairs.index(pair)
            except ValueError as exc:
                raise ValueError(
                    f"wavefield does not contain station pair {pair}"
                ) from exc

            trace = data[pair_index]
            trace_scale = np.max(np.abs(trace))
            normalized_trace = (
                trace / trace_scale if trace_scale > 0 else np.zeros_like(trace)
            )
            ax.plot(
                self.time,
                plotted_distances[pair_index]
                + 0.2 * np.max(distances) * normalized_trace,
                color="black",
                linewidth=1.0,
            )
        ax.set_xlim(tshow_min, tshow_max)
        ax.set(
            xlabel="Correlation time (s)",
            ylabel="Distance (km)",
            title=f"Wavefield ({self.component})" if self.component else "Wavefield",
        )

        return fig, ax

    def copy(self) -> "Wavefield":
        """Return an independent wavefield with copied traces and metadata."""
        return Wavefield(
            self._stream,
            component=self.component,
            check_distance_order=self.check_distance_order,
        )

    def preprocess(
        self,
        fmin: float = 0.5,
        fmax: float = 5.0,
        vmin: float = 0.1,
        vmax: float = 2.5,
        taper_fraction: float = 0.05,
    ) -> "Wavefield":
        """Preprocess and replace this wavefield in place.

        This is the mutating counterpart of :func:`preprocess_stream`.  The
        processing is completed on an independent copy first, so a validation
        or filtering error leaves the current wavefield unchanged.  The
        method returns ``self`` for optional method chaining.
        """
        from .preprocessing import preprocess_stream

        processed = preprocess_stream(
            self,
            fmin=fmin,
            fmax=fmax,
            vmin=vmin,
            vmax=vmax,
            taper_fraction=taper_fraction,
        )
        self._stream = processed.stream()
        self._component = processed.component
        self._check_distance_order = processed.check_distance_order
        self._pairs = processed.pairs
        return self

    def __getitem__(self, index):
        """Return the underlying trace or trace slice for sequence-style access."""
        return self._stream[index]

    def __iter__(self):
        """Iterate over the stored ObsPy traces in canonical pair order."""
        return iter(self._stream)

    def __len__(self) -> int:
        """Return :attr:`n_pairs` for standard collection semantics."""
        return self.n_pairs

    def __repr__(self) -> str:
        """Return a compact representation with geometry and component metadata."""
        name = f", component={self.component!r}" if self.component else ""
        return (
            f"Wavefield(n_pairs={self.n_pairs}, n_samples={self.n_samples}, "
            f"sampling_rate={self.sampling_rate:g}{name})"
        )


__all__ = ["Wavefield"]
