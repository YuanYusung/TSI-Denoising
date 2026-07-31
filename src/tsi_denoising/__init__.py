# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Public API for dense-linear-array surface-wave processing.

The package reads SAC cross-correlations into :class:`Wavefield` objects and
provides preprocessing, phase-shift MASW imaging, modal separation, and
three-station-interferometry (TSI) denoising.  Public symbols are re-exported
here so applications should import from :mod:`tsi_denoising` rather than from
implementation modules.

Notes
-----
Distances are expressed in km, velocities in km/s, correlation times in s,
and frequencies in Hz unless an individual API documents otherwise.
"""

from .io import read_sac_directory
from .masw import MASW, compute_masw
from .mode_separation import phase_match_separate, polarization_separate
from .preprocessing import preprocess_stream
from .denoising import (
    DenoisingResult,
    denoise_station_pair_demo,
    denoise_wavefield_iteratively,
    plot_denoised_result,
)
from .wavefield import Wavefield

__all__ = [
    "Wavefield",
    "MASW",
    "compute_masw",
    "polarization_separate",
    "phase_match_separate",
    "preprocess_stream",
    "read_sac_directory",
    "DenoisingResult",
    "denoise_station_pair_demo",
    "denoise_wavefield_iteratively",
    "plot_denoised_result",
]
