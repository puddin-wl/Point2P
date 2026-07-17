"""Default configuration for RTAD MRAF/GS refinement — formal pipeline (f=429mm).

Target SLM: 1024×1024 @ 17 μm pixel pitch.
Refinement defaults reflect the verified best-params baseline
(artifacts/fixed_baseline_bg0p9_initial_compare_20260429-175148).
"""

from __future__ import annotations


CONFIG = {
    "physical": {
        "wavelength_m": 532e-9,
        "focal_length_m": 429e-3,
        "input_gaussian_1e2_diameter_m": 6e-3,  # 6mm baseline (2026-06-02), 原5mm
        "clear_aperture_m": 15e-3,
    },
    "grid": {
        "N": 2048,
        "focal_dx_um": 2.5,
        "focal_dy_um": 2.5,
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
    "refinement": {
        "method": "wgs",
        "num_iters": 200,
        "mraf_iters": 150,
        "wgs_iters": 50,
        "mraf_factor": 0.8,
        "wgs_after_iters": 0,
        "feedback_exponent": 0.7,
        "wgs_update_mask": "flat",
        "wgs_feedback": "amplitude",
        "wgs_strategy": "flat_local",
        "wgs_xy_iters": 20,
        "wgs_xonly_iters": 30,
        "wgs_feedback_exponent": 0.8,
        "wgs_xy_feedback_exponent": 0.3,
        "wgs_x_feedback_exponent": 0.45,
        "wgs_update_every": 5,
        "wgs_xy_update_every": 5,
        "wgs_x_update_every": 5,
        "wgs_weight_min": 0.5,
        "wgs_weight_max": 1.5,
        "wgs_xy_weight_min": 0.5,
        "wgs_xy_weight_max": 2.0,
        "wgs_x_weight_min": 0.5,
        "wgs_x_weight_max": 2.5,
        "wgs_normalize_weights": True,
        "wgs_x_normalize": True,
        "wgs_update_region": "flat",
        "bg_mode": "attenuate",
        "bg_factor": 0.9,
        "wgs_clip_min": 0.5,
        "wgs_clip_max": 2.0,
        "metrics_interval": 10,
    },
    "paths": {
        "phase_mat": None,
        "phase_var": "phase0",
        "out_root": "artifacts",
        "transpose_h5": False,
        "swap_phase_xy": True,
    },
    "runtime": {
        "use_cupy": True,
        "device_id": 0,
        "random_seed": 1,
        "smoke_shape": 512,
        "smoke_phase": "random",
        "figure_dpi": 150,
    },
    "slm": {
        "width": 1024,
        "height": 1024,
        "pitch_um": 17.0,
    },
}
