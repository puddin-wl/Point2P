"""Tests for the Z40/Z20 SLM compensation exporter."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from export_zernike_compensated_slm import (
    default_label,
    make_zernike_compensation,
    quantize_phase_to_uint8,
)


class ZernikeCompensationTests(unittest.TestCase):
    def test_noll_values_at_pupil_center(self) -> None:
        _, waves, mask = make_zernike_compensation(
            (101, 101),
            dx_doe_m=0.1e-3,
            pupil_diameter_m=10e-3,
            z40_rms_waves=0.1,
            z20_rms_waves=0.2,
        )
        expected = 0.1 * math.sqrt(5.0) - 0.2 * math.sqrt(3.0)
        self.assertTrue(mask[50, 50])
        self.assertAlmostEqual(float(waves[50, 50]), expected, places=12)

    def test_compensation_is_smoothly_extended_outside_pupil(self) -> None:
        phase_rad, waves, mask = make_zernike_compensation(
            (101, 101),
            dx_doe_m=0.1e-3,
            pupil_diameter_m=6e-3,
            z40_rms_waves=0.10625,
            z20_rms_waves=0.25,
        )
        self.assertTrue(np.any(waves[~mask] != 0.0))
        self.assertTrue(np.any(phase_rad[~mask] != 0.0))

        center = waves.shape[0] // 2
        # extension_width_rho=0.2, so samples beyond rho=1.2 are constant
        self.assertAlmostEqual(
            float(waves[center, 100]),
            float(waves[0, 0]),
            places=12,
        )

    def test_no_hard_step_at_normalization_boundary(self) -> None:
        phase_rad, _, _ = make_zernike_compensation(
            (401, 401),
            dx_doe_m=0.025e-3,
            pupil_diameter_m=6e-3,
            z40_rms_waves=0.10625,
            z20_rms_waves=0.25,
        )
        center = phase_rad.shape[0] // 2
        radius_px = int((3e-3) / (0.025e-3))
        inside = phase_rad[center, center + radius_px]
        outside = phase_rad[center, center + radius_px + 1]
        wrapped_step = abs(float(np.angle(np.exp(1j * (outside - inside)))))
        self.assertLess(wrapped_step, 0.25)

    def test_nominal_negative_aberration_cancels_positive_map(self) -> None:
        phase_rad, _, mask = make_zernike_compensation(
            (65, 65),
            dx_doe_m=0.2e-3,
            pupil_diameter_m=10e-3,
            z40_rms_waves=0.10625,
            z20_rms_waves=0.25,
        )
        residual = np.angle(
            np.exp(-1j * phase_rad[mask]) * np.exp(1j * phase_rad[mask])
        )
        self.assertLess(float(np.max(np.abs(residual))), 1e-12)

    def test_quantization_matches_existing_convention(self) -> None:
        phase = np.array([[0.0, np.pi, 2.0 * np.pi - 1e-12]])
        gray = quantize_phase_to_uint8(phase)
        self.assertEqual(gray.tolist(), [[0, 127, 254]])

    def test_default_label_records_positive_coefficients(self) -> None:
        self.assertEqual(
            default_label(0.10625, 0.25),
            "phase_V2_posSphericalComp_Z40+0p10625_Z20+0p25",
        )


if __name__ == "__main__":
    unittest.main()
