"""Default configuration for real-world DOE tolerance simulation.

This project intentionally treats the previous refinement output as read-only
input data. It does not import or modify the refinement program.
"""

from __future__ import annotations

from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
BASELINE_ARTIFACT_DIR = Path(
    r"E:\program\Point2P\rtad_mraf_gs_python\artifacts\fixed_baseline_bg0p9_initial_compare_20260429-175148"
)


CONFIG = {
    "paths": {
        "baseline_artifact_dir": BASELINE_ARTIFACT_DIR,
        "phase_refined": BASELINE_ARTIFACT_DIR / "phase_refined.npy",
        "target_npz": BASELINE_ARTIFACT_DIR / "target.npz",
        "baseline_config": BASELINE_ARTIFACT_DIR / "config_used.json",
        "output_root": PROJECT_DIR / "artifacts",
    },
    "physical": {
        "wavelength_m": 532e-9,
        "focal_length_m": 429e-3,
        "input_gaussian_1e2_diameter_m": 5e-3,
        "clear_aperture_m": 15e-3,
    },
    "grid": {
        "N": 2048,
        "focal_dx_um": 2.5,
        "focal_dy_um": 2.5,
        "dx_doe_m": 4.457578125e-05,
    },
    "target": {
        "W50_um": 330.0,
        "H50_um": 120.0,
        "delta_x_um": 15.0,
        "delta_y_um": 8.0,
        "guard_x_um": 20.0,
        "guard_y_um": 12.0,
        "constraint_mode": "truncated_rtad",
        "release_level": 0.1353352832366127,
        "target_mode": "separable",
    },
    "nominal": {
        "defocus_mm": 0.0,
        "offset_x_mm": 0.0,
        "offset_y_mm": 0.0,
        "diameter_1e2_x_mm": 5.0,
        "diameter_1e2_y_mm": 5.0,
        "divergence_edge_mrad": 0.0,
        "pointing_shift_x_um": 0.0,
        "pointing_shift_y_um": 0.0,
        "clear_aperture_mm": 15.0,
    },
    "runtime": {
        "figure_dpi": 150,
        "smoke_size": 512,
    },
}


STRESS_SWEEPS = {
    "nominal": [{"parameter": "nominal", "value": 0.0, "updates": {}}],
    "defocus": [
        {"parameter": "defocus_mm", "value": v, "updates": {"defocus_mm": v}}
        for v in [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
    ],
    "beam_offset_x": [
        {"parameter": "offset_x_mm", "value": v, "updates": {"offset_x_mm": v}}
        for v in [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
    ],
    "beam_offset_y": [
        {"parameter": "offset_y_mm", "value": v, "updates": {"offset_y_mm": v}}
        for v in [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
    ],
    "beam_size": [
        {
            "parameter": "diameter_1e2_mm",
            "value": v,
            "updates": {"diameter_1e2_x_mm": v, "diameter_1e2_y_mm": v},
        }
        for v in [3.5, 4, 4.5, 5, 5.5, 6, 6.5]
    ],
    "divergence": [
        {"parameter": "divergence_edge_mrad", "value": v, "updates": {"divergence_edge_mrad": v}}
        for v in [-0.5, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.5]
    ],
    "pointing_x": [
        {
            "parameter": "pointing_shift_x_um",
            "value": v,
            "updates": {"pointing_shift_x_um": v},
        }
        for v in [-100, -50, -25, 0, 25, 50, 100]
    ],
    "pointing_y": [
        {
            "parameter": "pointing_shift_y_um",
            "value": v,
            "updates": {"pointing_shift_y_um": v},
        }
        for v in [-100, -50, -25, 0, 25, 50, 100]
    ],
    "aperture": [
        {"parameter": "clear_aperture_mm", "value": v, "updates": {"clear_aperture_mm": v}}
        for v in [10, 12, 13, 14, 15]
    ],
    "ellipticity": [
        {
            "parameter": "diameter_1e2_xy_mm",
            "value": [dx, dy],
            "updates": {"diameter_1e2_x_mm": dx, "diameter_1e2_y_mm": dy},
        }
        for dx, dy in [
            (5.0, 5.0),
            (4.5, 5.0),
            (5.5, 5.0),
            (5.0, 4.5),
            (5.0, 5.5),
            (4.5, 5.5),
            (5.5, 4.5),
        ]
    ],
}


MILD_SWEEPS = {
    "nominal": [{"parameter": "nominal", "value": 0.0, "updates": {}}],
    "defocus": [
        {"parameter": "defocus_mm", "value": v, "updates": {"defocus_mm": v}}
        for v in [-0.5, -0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.25, 0.5]
    ],
    "beam_offset_x": [
        {"parameter": "offset_x_mm", "value": v, "updates": {"offset_x_mm": v}}
        for v in [-0.5, -0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.25, 0.5]
    ],
    "beam_offset_y": [
        {"parameter": "offset_y_mm", "value": v, "updates": {"offset_y_mm": v}}
        for v in [-0.5, -0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.25, 0.5]
    ],
    "beam_size": [
        {
            "parameter": "diameter_1e2_mm",
            "value": v,
            "updates": {"diameter_1e2_x_mm": v, "diameter_1e2_y_mm": v},
        }
        for v in [4.5, 4.75, 5, 5.25, 5.5]
    ],
    "divergence": [
        {"parameter": "divergence_edge_mrad", "value": v, "updates": {"divergence_edge_mrad": v}}
        for v in [-0.1, -0.05, -0.02, -0.01, 0, 0.01, 0.02, 0.05, 0.1]
    ],
    "pointing_x": [
        {
            "parameter": "pointing_shift_x_um",
            "value": v,
            "updates": {"pointing_shift_x_um": v},
        }
        for v in [-25, -10, -5, 0, 5, 10, 25]
    ],
    "pointing_y": [
        {
            "parameter": "pointing_shift_y_um",
            "value": v,
            "updates": {"pointing_shift_y_um": v},
        }
        for v in [-25, -10, -5, 0, 5, 10, 25]
    ],
    "aperture": [
        {"parameter": "clear_aperture_mm", "value": v, "updates": {"clear_aperture_mm": v}}
        for v in [13, 14, 14.5, 15]
    ],
    "ellipticity": [
        {
            "parameter": "diameter_1e2_xy_mm",
            "value": [dx, dy],
            "updates": {"diameter_1e2_x_mm": dx, "diameter_1e2_y_mm": dy},
        }
        for dx, dy in [
            (5.0, 5.0),
            (4.75, 5.0),
            (5.25, 5.0),
            (5.0, 4.75),
            (5.0, 5.25),
            (4.75, 5.25),
            (5.25, 4.75),
        ]
    ],
}


SWEEP_PROFILES = {
    "mild": MILD_SWEEPS,
    "stress": STRESS_SWEEPS,
}


# Keep the original broad ranges available for compatibility.
SWEEPS = STRESS_SWEEPS


ALL_SWEEPS = [
    "defocus",
    "beam_offset_x",
    "beam_offset_y",
    "beam_size",
    "divergence",
    "pointing_x",
    "pointing_y",
    "aperture",
    "ellipticity",
]


SUMMARY_FIELDS = [
    "sweep_name",
    "sweep_parameter",
    "sweep_value",
    "size50_x_um",
    "size50_y_um",
    "size13p5_x_um",
    "size13p5_y_um",
    "transition_13p5_90_x_um",
    "transition_13p5_90_y_um",
    "rms_nonuniformity_percent",
    "efficiency_e2_percent",
    "aperture_throughput_percent",
    "center_offset_x_um",
    "center_offset_y_um",
]
