# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Validation and metadata-normalization helpers for surface-wave streams.

These internal routines establish the package-wide contract for ObsPy trace
metadata, station-pair ordering, finite waveform samples, and common sampling
geometry.  They never mutate an input stream except where an explicitly named
normalization helper returns a copied and normalized stream.

Notes
-----
Station order is inferred from the final contiguous integer in a station name.
Distances are km; sampling interval and correlation-time origin are seconds.
"""

from __future__ import annotations

import re
import warnings
from collections import Counter

import numpy as np
from obspy import Stream


def _nested_attribute(obj, path: str):
    """Read a dotted attribute path, returning None if any part is missing."""
    value = obj
    for name in path.split("."):
        value = getattr(value, name, None)
        if value is None:
            return None
    return value


def _station_name_from_fields(
    trace,
    fields: tuple[str, ...],
    *,
    role: str,
) -> str:
    """Return the first non-empty station name found in candidate fields."""
    attempted: list[str] = []

    for field in fields:
        attempted.append(field)
        value = _nested_attribute(trace, field)

        if value is None:
            continue

        name = str(value).strip()
        if name:
            return name

    field_list = ", ".join(attempted)
    raise ValueError(
        f"trace is missing a valid {role} station name; "
        f"checked fields: {field_list}"
    )


SOURCE_FIELDS = (
    "stats.sac.kevnm",
    "stats.source",
)

RECEIVER_FIELDS = (
    "stats.sac.kstnm",
    "stats.station",
)


def trace_pair(trace) -> tuple[str, str]:
    """Return the source-receiver pair stored in trace metadata."""
    source = _station_name_from_fields(
        trace,
        SOURCE_FIELDS,
        role="source",
    )
    receiver = _station_name_from_fields(
        trace,
        RECEIVER_FIELDS,
        role="receiver",
    )
    return source, receiver


_STATION_NUMBER_PATTERN = re.compile(r"\d+")


def station_number(name: str) -> int:
    """Extract the station number from a station name.

    The last contiguous digit group is used, so names such as
    ``LINE1_STA023`` are interpreted as station 23.
    """
    name = str(name).strip()
    matches = _STATION_NUMBER_PATTERN.findall(name)

    if not matches:
        raise ValueError(
            f"station name {name!r} does not contain a numeric identifier"
        )

    return int(matches[-1])


def validate_station_numbers(
    pairs: tuple[tuple[str, str], ...],
) -> dict[str, int]:
    """Validate station names and return their numeric identifiers."""
    names = {
        name
        for pair in pairs
        for name in pair
    }

    numbers = {
        name: station_number(name)
        for name in names
    }

    number_to_names: dict[int, list[str]] = {}
    for name, number in numbers.items():
        number_to_names.setdefault(number, []).append(name)

    conflicts = {
        number: sorted(names)
        for number, names in number_to_names.items()
        if len(names) > 1
    }

    if conflicts:
        raise ValueError(
            "multiple station names map to the same numeric identifier: "
            f"{conflicts}"
        )

    return numbers


def _reverse_pair_trace(trace, source: str, receiver: str):
    """Return a reversed copy with canonical source-receiver metadata."""
    reversed_trace = trace.copy()
    sac = getattr(reversed_trace.stats, "sac", None)
    if sac is None:
        raise ValueError(
            f"cannot reverse station pair {(source, receiver)} without SAC headers"
        )

    had_stats_source = hasattr(reversed_trace.stats, "source")
    had_stats_station = hasattr(reversed_trace.stats, "station")
    reversed_trace.data = np.asarray(reversed_trace.data)[::-1].copy()
    sac.kevnm = receiver
    sac.kstnm = source
    if had_stats_source:
        reversed_trace.stats.source = receiver
    if had_stats_station:
        reversed_trace.stats.station = source
    return reversed_trace


def normalize_pair_directions(stream: Stream) -> Stream:
    """Remove self pairs and normalize cross-correlation directions.

    Reverse duplicates retain the trace whose source has the smaller station
    number. A lone trace whose source number is larger than its receiver
    number is reversed, including its sample order, while ``sac.b`` and
    ``sac.dist`` remain unchanged.
    """
    if not isinstance(stream, Stream):
        raise TypeError("stream must be an obspy.Stream")
    if len(stream) == 0:
        raise ValueError("stream is empty")

    pairs = tuple(trace_pair(trace) for trace in stream)
    numbers = validate_station_numbers(pairs)
    pair_counts = Counter(pairs)
    duplicates = sorted(pair for pair, count in pair_counts.items() if count > 1)
    if duplicates:
        raise ValueError(f"stream contains duplicate station pairs: {duplicates}")

    grouped: dict[
        tuple[int, int],
        list[tuple[int, object, str, str, int, int]],
    ] = {}
    for index, (trace, (source, receiver)) in enumerate(zip(stream, pairs)):
        source_number = numbers[source]
        receiver_number = numbers[receiver]
        if source == receiver:
            warnings.warn(
                f"dropping self-pair trace {source!r} -> {receiver!r}",
                UserWarning,
                stacklevel=2,
            )
            continue

        key = tuple(sorted((source_number, receiver_number)))
        grouped.setdefault(key, []).append(
            (
                index,
                trace,
                source,
                receiver,
                source_number,
                receiver_number,
            )
        )

    retained: list[tuple[int, object]] = []
    for records in grouped.values():
        if len(records) == 2:
            preferred = min(records, key=lambda record: record[4])
            discarded = max(records, key=lambda record: record[4])
            warnings.warn(
                "dropping reverse duplicate "
                f"{discarded[2]!r} -> {discarded[3]!r}; "
                f"keeping {preferred[2]!r} -> {preferred[3]!r} because its "
                "source station number is smaller",
                UserWarning,
                stacklevel=2,
            )
            retained.append((preferred[0], preferred[1]))
            continue

        if len(records) != 1:
            raise ValueError("station pair direction normalization is ambiguous")

        index, trace, source, receiver, source_number, receiver_number = records[0]
        if source_number < receiver_number:
            retained.append((index, trace))
            continue

        warnings.warn(
            f"reversing station pair {source!r} -> {receiver!r} to "
            f"{receiver!r} -> {source!r}",
            UserWarning,
            stacklevel=2,
        )
        retained.append(
            (index, _reverse_pair_trace(trace, source, receiver))
        )

    if not retained:
        raise ValueError("no usable interstation traces remain")

    retained.sort(key=lambda item: item[0])
    return Stream(traces=[trace for _, trace in retained])


def trace_distance(trace) -> float:
    """Return a finite, positive interstation distance in kilometres."""
    sac = getattr(trace.stats, "sac", None)
    if sac is None or not hasattr(sac, "dist"):
        raise ValueError("trace is missing required SAC header sac.dist")
    distance = float(sac.dist)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("SAC distance must be finite and positive")
    return distance


def trace_start_time(trace) -> float:
    """Return the SAC correlation-time origin in seconds."""
    sac = getattr(trace.stats, "sac", None)
    if sac is None or not hasattr(sac, "b"):
        raise ValueError("trace is missing required SAC header sac.b")
    start = float(sac.b)
    if not np.isfinite(start):
        raise ValueError("SAC start time must be finite")
    return start


def validate_n_jobs(n_jobs: int) -> int:
    """Return a validated positive process count."""
    if isinstance(n_jobs, bool) or not isinstance(n_jobs, (int, np.integer)):
        raise TypeError("n_jobs must be an integer")
    if n_jobs < 1:
        raise ValueError("n_jobs must be at least 1")
    return int(n_jobs)


def validate_stream(stream: Stream, *, require_unique_pairs: bool = True) -> None:
    """Validate metadata and common sampling for a non-empty ObsPy Stream."""
    if not isinstance(stream, Stream):
        raise TypeError("stream must be an obspy.Stream")
    if len(stream) == 0:
        raise ValueError("stream is empty")

    reference = stream[0]
    reference_npts = int(reference.stats.npts)
    reference_delta = float(reference.stats.delta)
    reference_start = trace_start_time(reference)
    if reference_npts <= 0:
        raise ValueError("traces must contain at least one sample")
    if not np.isfinite(reference_delta) or reference_delta <= 0:
        raise ValueError("sampling interval must be finite and positive")

    pairs: list[tuple[str, str]] = []
    for trace in stream:
        pair = trace_pair(trace)
        trace_distance(trace)
        start = trace_start_time(trace)
        delta = float(trace.stats.delta)
        data = np.asarray(trace.data)

        if int(trace.stats.npts) != reference_npts:
            raise ValueError(
                f"trace {pair} has {trace.stats.npts} samples; "
                f"expected {reference_npts}"
            )
        if not np.isclose(delta, reference_delta, rtol=0.0, atol=1e-12):
            raise ValueError(
                f"trace {pair} has sampling interval {delta}; "
                f"expected {reference_delta}"
            )
        if not np.isclose(start, reference_start, rtol=0.0, atol=1e-9):
            raise ValueError(
                f"trace {pair} starts at {start}; expected {reference_start}"
            )
        if data.ndim != 1 or data.size != reference_npts:
            raise ValueError(f"trace {pair} data must be one-dimensional")
        if not np.all(np.isfinite(data)):
            raise ValueError(f"trace {pair} contains non-finite samples")
        pairs.append(pair)

    if require_unique_pairs and len(set(pairs)) != len(pairs):
        pair_counts = Counter(pairs)
        duplicates = sorted(
            pair for pair, count in pair_counts.items() if count > 1
        )
        raise ValueError(f"stream contains duplicate station pairs: {duplicates}")


def validate_reference_distance_order(stream: Stream) -> None:
    """Require distances from the minimum station to increase by station number."""
    if not isinstance(stream, Stream):
        raise TypeError("stream must be an obspy.Stream")
    if len(stream) == 0:
        raise ValueError("stream is empty")

    pairs = tuple(trace_pair(trace) for trace in stream)
    numbers = validate_station_numbers(pairs)
    reference_name = min(numbers, key=numbers.__getitem__)
    reference_number = numbers[reference_name]
    records: list[tuple[int, str, float]] = []

    for trace, (source, receiver) in zip(stream, pairs):
        if source == reference_name:
            other_name = receiver
        elif receiver == reference_name:
            other_name = source
        else:
            continue
        records.append(
            (numbers[other_name], other_name, trace_distance(trace))
        )

    records.sort(key=lambda record: record[0])
    for previous, current in zip(records, records[1:]):
        previous_number, previous_name, previous_distance = previous
        current_number, current_name, current_distance = current
        if current_distance <= previous_distance:
            raise ValueError(
                "distance ordering is inconsistent with station numbering: "
                f"reference station {reference_name!r} "
                f"(number {reference_number}), distance to "
                f"{previous_name!r} (number {previous_number}) is "
                f"{previous_distance:g} km, but distance to "
                f"{current_name!r} (number {current_number}) is "
                f"{current_distance:g} km"
            )


def ordered_pairs(stream: Stream) -> tuple[tuple[str, str], ...]:
    """Return station-pair keys in stream order after validation."""
    validate_stream(stream)
    return tuple(trace_pair(trace) for trace in stream)


def validate_frequency_band(
    fmin: float, fmax: float, sampling_rate: float
) -> tuple[float, float]:
    """Validate a positive frequency band below Nyquist."""
    fmin = float(fmin)
    fmax = float(fmax)
    sampling_rate = float(sampling_rate)
    nyquist = sampling_rate / 2.0
    if not np.isfinite(fmin) or not np.isfinite(fmax):
        raise ValueError("frequency bounds must be finite")
    if not 0 < fmin < fmax:
        raise ValueError("frequency bounds must satisfy 0 < fmin < fmax")
    if fmax >= nyquist:
        raise ValueError(
            f"fmax must be below the Nyquist frequency ({nyquist:g} Hz)"
        )
    return fmin, fmax


__all__ = [
    "normalize_pair_directions",
    "ordered_pairs",
    "station_number",
    "trace_distance",
    "trace_pair",
    "trace_start_time",
    "validate_frequency_band",
    "validate_reference_distance_order",
    "validate_station_numbers",
    "validate_stream",
]
