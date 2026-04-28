"""Default configuration for one RTAD MRAF/GS refinement run."""

from __future__ import annotations


CONFIG = {
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
    },
    "target": {
        "W50_um": 330.0,
        "H50_um": 120.0,
        "delta_x_um": 15.0,
        "delta_y_um": 8.0,
        "guard_x_um": 20.0,
        "guard_y_um": 12.0,
        "mode": "separable",
    },
    "refinement": {
        "method": "mraf",
        "num_iters": 200,
        "mraf_factor": 0.4,
        "wgs_after_iters": 0,
        "feedback_exponent": 0.7,
        "wgs_update_region": "flat",
        "bg_mode": "attenuate",
        "bg_factor": 0.25,
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
}
