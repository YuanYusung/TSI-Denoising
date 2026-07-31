# Copyright (c) 2026 Yusong Yuan and Hongrui Qiu
#
# This file is part of TSI-Denoising.

"""Public three-station-interferometry denoising API.

The subpackage exposes iterative wavefield denoising, publication-oriented
diagnostic plots, and the serializable :class:`DenoisingResult` container.

References
----------
Qiu, H., Niu, F., and Qin, L. (2021). Denoising surface waves extracted from
ambient noise recorded by 1-D linear array using three-station
interferometry of direct waves. *JGR: Solid Earth*, 126, e2021JB021712.
https://doi.org/10.1029/2021JB021712
"""

from .diagnostics import denoise_station_pair_demo, plot_denoised_result
from .three_station import DenoisingResult, denoise_wavefield_iteratively

__all__ = [
    "DenoisingResult",
    "denoise_station_pair_demo",
    "plot_denoised_result",
    "denoise_wavefield_iteratively",
]
