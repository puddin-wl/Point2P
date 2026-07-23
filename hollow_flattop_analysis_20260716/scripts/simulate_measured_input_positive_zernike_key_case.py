"""Run the measured-Gaussian key case with both Zernike signs reversed."""

from __future__ import annotations

import simulate_measured_input_zernike_key_case as simulation


simulation.Z40_RMS_WAVES = +0.10625
simulation.Z20_RMS_WAVES = +0.25000
simulation.OUTPUT_DIR = (
    simulation.ANALYSIS_ROOT
    / "results"
    / "13_measured_input_positive_zernike_key_case"
)


if __name__ == "__main__":
    simulation.main()
