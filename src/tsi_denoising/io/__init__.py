# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Public waveform-input helpers.

Currently the subpackage provides recursive SAC ingestion into a validated
:class:`~tsi_denoising.Wavefield` with normalized station-pair directions.
"""

from .sac import read_sac_directory

__all__ = ["read_sac_directory"]
