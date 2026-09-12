"""Tests for the compensated whole-phase Y sweep exporter."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from export_zernike_compensated_y_sweep import (
    output_filename,
    parse_shift_values,
)


class CompensatedYSweepTests(unittest.TestCase):
    def test_parse_shift_values_preserves_loading_order(self) -> None:
        self.assertEqual(
            parse_shift_values("-10,-5,0,5,10,15,20"),
            [-10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0],
        )

    def test_duplicate_shift_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            parse_shift_values("0,5,5")

    def test_output_filename_is_explicit(self) -> None:
        self.assertEqual(
            output_filename(5.0, -10.0),
            "phase_V2_posZernikeComp_shiftX+5_Y-10.bmp",
        )


if __name__ == "__main__":
    unittest.main()
