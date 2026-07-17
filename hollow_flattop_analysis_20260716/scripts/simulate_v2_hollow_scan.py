"""Reproduce the 2026-07-16 flat-top center depression with the fixed V2 phase.

The script deliberately reuses the production forward propagator from
``rtad_mraf_gs_python_test_20260605/src/propagation.py``.  It does not modify
the V2 artifact or any real-test data.

The experimental target is the average of all frames in 20260716-1.bmData.
Simulation and experiment are compared in physical units.  In particular,
the experimental 5 x 5 center statistic covers 36.9 x 36.9 um, so the
simulation uses the same physical window rather than five 2.5-um pixels.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
POINT2P_ROOT = ANALYSIS_ROOT.parent
PROJECT_ROOT = POINT2P_ROOT / "rtad_mraf_gs_python_test_20260605"
REAL_TEST_ROOT = POINT2P_ROOT / "real_test"
DEFAULT_CASE_DIR = (
    PROJECT_ROOT
    / "artifacts"
    / "run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm"
)
DEFAULT_BMDATA = REAL_TEST_ROOT / "20260716" / "20260716-1.bmData"
DEFAULT_PERCENT_SUMMARY = (
    REAL_TEST_ROOT
    / "20260716"
    / "analysis_20260716_1_main_spot_only"
    / "20260716-1_percent_energy_summary.json"
)
DEFAULT_HOLLOW_SUMMARY = (
    REAL_TEST_ROOT
    / "20260716"
    / "analysis_20260716_1_main_spot_only_sequence"
    / "20260716-1_hollow_sequence_summary.json"
)

for module_root in (PROJECT_ROOT, REAL_TEST_ROOT):
    value = str(module_root)
    if value not in sys.path:
        sys.path.insert(0, value)

from analyze_rect_flattop_size import load_spiricon_frame, robust_background  # noqa: E402
from src.diagnostics import _load_phase_from_case, load_case_data  # noqa: E402
from src.propagation import forward_fft, intensity, make_input_gaussian  # noqa: E402


def normalize_power(field: np.ndarray) -> np.ndarray:
    """Return a complex64 or float32 field with unit discrete power."""
    power = float(np.sum(np.abs(field) ** 2, dtype=np.float64))
    if not np.isfinite(power) or power <= 0:
        raise ValueError(f"Input field power must be positive; got {power}")
    return (field / math.sqrt(power)).astype(field.dtype, copy=False)


def load_experimental_target(
    bmdata: Path,
    percent_summary_path: Path,
    hollow_summary_path: Path,
) -> dict[str, Any]:
    """Load and average every BeamGage frame, preserving physical coordinates."""
    percent_summary = json.loads(percent_summary_path.read_text(encoding="utf-8"))
    hollow_summary = json.loads(hollow_summary_path.read_text(encoding="utf-8"))
    with h5py.File(bmdata, "r") as h5:
        frame_keys = sorted(h5["BG_DATA"].keys(), key=int)

    signals: list[np.ndarray] = []
    frame_maxima: list[float] = []
    for key in frame_keys:
        image, _ = load_spiricon_frame(bmdata, frame_key=key)
        background, _ = robust_background(image, corner_px=30)
        signal = np.clip(image - background, 0.0, None)
        signals.append(signal.astype(np.float32))
        frame_maxima.append(float(np.max(image)))

    average = np.mean(np.stack(signals, axis=0), axis=0, dtype=np.float64)
    sx = float(percent_summary["metadata"]["pixel_scale_x_um"])
    sy = float(percent_summary["metadata"]["pixel_scale_y_um"])
    cx, cy = (float(v) for v in percent_summary["center_px"])
    x0, y0, width, height = (int(v) for v in hollow_summary["size50_analysis_box_px"])
    x1, y1 = x0 + width, y0 + height
    core = average[y0:y1, x0:x1]
    core_mean = float(np.mean(core))
    normalized = average / core_mean
    x_rel_um = (np.arange(average.shape[1], dtype=np.float64) - cx) * sx
    y_rel_um = (np.arange(average.shape[0], dtype=np.float64) - cy) * sy
    center_width_x_um = 5.0 * sx
    center_width_y_um = 5.0 * sy

    metrics = morphology_metrics(
        normalized=normalized,
        x_um=x_rel_um,
        y_um=y_rel_um,
        core_width_x_um=width * sx,
        core_width_y_um=height * sy,
        center_width_x_um=center_width_x_um,
        center_width_y_um=center_width_y_um,
    )
    metrics["aggregate_from_existing_analysis"] = hollow_summary["aggregate"]
    return {
        "average_signal": average,
        "normalized": normalized,
        "x_um": x_rel_um,
        "y_um": y_rel_um,
        "frame_count": len(frame_keys),
        "frame_maxima": frame_maxima,
        "pixel_scale_um": [sx, sy],
        "center_px": [cx, cy],
        "core_box_px": [x0, y0, width, height],
        "core_width_um": [width * sx, height * sy],
        "center_window_um": [center_width_x_um, center_width_y_um],
        "metrics": metrics,
    }


def band_profile(
    normalized: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
    axis: str,
    band_width_um: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return a center-band mean profile in physical coordinates."""
    if axis == "x":
        rows = np.abs(y_um) <= band_width_um / 2.0
        return x_um, np.mean(normalized[rows, :], axis=0)
    if axis == "y":
        columns = np.abs(x_um) <= band_width_um / 2.0
        return y_um, np.mean(normalized[:, columns], axis=1)
    raise ValueError(f"Unsupported profile axis: {axis}")


def profile_on_unit_core(
    normalized: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
    core_width_x_um: float,
    core_width_y_um: float,
    center_width_x_um: float,
    center_width_y_um: float,
    samples: int = 81,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate x/y center-band profiles onto the common unit core."""
    unit = np.linspace(-1.0, 1.0, samples)
    xp_x, profile_x = band_profile(normalized, x_um, y_um, "x", center_width_y_um)
    xp_y, profile_y = band_profile(normalized, x_um, y_um, "y", center_width_x_um)
    common_x = unit * core_width_x_um / 2.0
    common_y = unit * core_width_y_um / 2.0
    return (
        unit,
        np.interp(common_x, xp_x, profile_x),
        np.interp(common_y, xp_y, profile_y),
    )


def morphology_metrics(
    normalized: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
    core_width_x_um: float,
    core_width_y_um: float,
    center_width_x_um: float,
    center_width_y_um: float,
) -> dict[str, Any]:
    """Measure the same physical center and core regions in both datasets."""
    x_core = np.abs(x_um) <= core_width_x_um / 2.0
    y_core = np.abs(y_um) <= core_width_y_um / 2.0
    core = normalized[np.ix_(y_core, x_core)]
    x_center = np.abs(x_um) <= center_width_x_um / 2.0
    y_center = np.abs(y_um) <= center_width_y_um / 2.0
    center = normalized[np.ix_(y_center, x_center)]

    core_x = x_um[x_core]
    middle_columns = np.abs(core_x) <= core_width_x_um / 6.0
    side_columns = ~middle_columns
    middle = core[:, middle_columns]
    sides = core[:, side_columns]
    top = core[: max(1, core.shape[0] // 3), :]
    bottom = core[-max(1, core.shape[0] // 3) :, :]
    left = core[:, : max(1, core.shape[1] // 3)]
    right = core[:, -max(1, core.shape[1] // 3) :]
    unit, profile_x, profile_y = profile_on_unit_core(
        normalized,
        x_um,
        y_um,
        core_width_x_um,
        core_width_y_um,
        center_width_x_um,
        center_width_y_um,
    )
    return {
        "core_mean": float(np.mean(core)),
        "center_window_over_core_mean": float(np.mean(center) / np.mean(core)),
        "middle_third_over_side_thirds": float(np.mean(middle) / np.mean(sides)),
        "core_rms_fraction": float(np.std(core) / np.mean(core)),
        "top_over_bottom": float(np.mean(top) / np.mean(bottom)),
        "left_over_right": float(np.mean(left) / np.mean(right)),
        "unit_profile_axis": unit.tolist(),
        "unit_profile_x": profile_x.tolist(),
        "unit_profile_y": profile_y.tolist(),
    }


def candidate_score(metrics: dict[str, Any], target: dict[str, Any]) -> dict[str, float]:
    """Weighted mismatch score; lower is better."""
    profile_x = np.asarray(metrics["unit_profile_x"], dtype=np.float64)
    profile_y = np.asarray(metrics["unit_profile_y"], dtype=np.float64)
    target_x = np.asarray(target["unit_profile_x"], dtype=np.float64)
    target_y = np.asarray(target["unit_profile_y"], dtype=np.float64)
    central = slice(8, -8)
    rmse_x = float(np.sqrt(np.mean((profile_x[central] - target_x[central]) ** 2)))
    rmse_y = float(np.sqrt(np.mean((profile_y[central] - target_y[central]) ** 2)))
    components = {
        "center": (metrics["center_window_over_core_mean"] - target["center_window_over_core_mean"]) / 0.04,
        "middle": (metrics["middle_third_over_side_thirds"] - target["middle_third_over_side_thirds"]) / 0.04,
        "rms": (metrics["core_rms_fraction"] - target["core_rms_fraction"]) / 0.05,
        "profile_x": rmse_x / 0.15,
        "profile_y": rmse_y / 0.15,
    }
    score = float(math.sqrt(sum(value * value for value in components.values())))
    return {
        "score": score,
        "profile_rmse_x": rmse_x,
        "profile_rmse_y": rmse_y,
        **{f"scaled_{key}": float(value) for key, value in components.items()},
    }


def build_scan_cases() -> list[dict[str, Any]]:
    """Construct a broad, bounded first-pass scan."""
    cases: list[dict[str, Any]] = [{"family": "baseline", "name": "baseline", "params": {}}]

    for family in ("defocus", "astigmatism", "spherical"):
        for value in np.linspace(-1.6, 1.6, 17):
            cases.append(
                {
                    "family": family,
                    "name": f"{family}_{value:+.3f}waves",
                    "params": {f"{family}_waves": float(value)},
                }
            )

    for value in np.linspace(-1.2, 1.2, 13):
        if abs(value) < 1e-12:
            continue
        cases.append(
            {
                "family": "beam_offset_x",
                "name": f"beam_offset_x_{value:+.3f}mm",
                "params": {"beam_offset_x_mm": float(value)},
            }
        )
        cases.append(
            {
                "family": "beam_offset_y",
                "name": f"beam_offset_y_{value:+.3f}mm",
                "params": {"beam_offset_y_mm": float(value)},
            }
        )

    for diameter in np.linspace(5.0, 14.0, 13):
        cases.append(
            {
                "family": "centered_aperture",
                "name": f"centered_aperture_{diameter:.3f}mm",
                "params": {"aperture_diameter_mm": float(diameter)},
            }
        )

    for offset in np.linspace(-2.0, 2.0, 9):
        if abs(offset) < 1e-12:
            continue
        cases.append(
            {
                "family": "aperture_offset_x",
                "name": f"aperture7mm_offset_x_{offset:+.3f}mm",
                "params": {
                    "aperture_diameter_mm": 7.0,
                    "aperture_offset_x_mm": float(offset),
                },
            }
        )
        cases.append(
            {
                "family": "aperture_offset_y",
                "name": f"aperture7mm_offset_y_{offset:+.3f}mm",
                "params": {
                    "aperture_diameter_mm": 7.0,
                    "aperture_offset_y_mm": float(offset),
                },
            }
        )

    depths = np.linspace(0.15, 0.90, 6)
    sigmas = np.linspace(0.35, 2.10, 6)
    for depth in depths:
        for sigma in sigmas:
            for family in ("center_dip", "horizontal_input_dip", "vertical_input_dip"):
                cases.append(
                    {
                        "family": family,
                        "name": f"{family}_d{depth:.3f}_s{sigma:.3f}mm",
                        "params": {
                            "amplitude_dip_family": family,
                            "amplitude_dip_depth": float(depth),
                            "amplitude_dip_sigma_mm": float(sigma),
                        },
                    }
                )

    for coefficient in np.linspace(-1.0, 1.0, 21):
        if abs(coefficient) < 1e-12:
            continue
        cases.append(
            {
                "family": "HG20_mix",
                "name": f"HG20_mix_{coefficient:+.3f}",
                "params": {"hg20_coefficient": float(coefficient)},
            }
        )
        cases.append(
            {
                "family": "HG02_mix",
                "name": f"HG02_mix_{coefficient:+.3f}",
                "params": {"hg02_coefficient": float(coefficient)},
            }
        )

    cases.append(
        {
            "family": "measured_ellipse",
            "name": "measured_Gaussian_X6p40_Y6p31mm",
            "params": {"diameter_x_mm": 6.40, "diameter_y_mm": 6.31},
        }
    )
    return cases


class Simulator:
    """Fixed V2 phase plus reusable input-plane arrays."""

    def __init__(self, case_dir: Path):
        self.case_dir = case_dir
        self.base = load_case_data(case_dir)
        self.phase = _load_phase_from_case(case_dir).astype(np.float32)
        self.config = self.base.config
        self.dx_m = float(self.config["grid"]["dx_doe_m"])
        self.diameter_m = float(self.config["physical"]["input_gaussian_1e2_diameter_m"])
        self.clear_aperture_m = float(self.config["physical"]["clear_aperture_m"])
        ny, nx = self.phase.shape
        x = (np.arange(nx, dtype=np.float32) - nx // 2) * np.float32(self.dx_m)
        y = (np.arange(ny, dtype=np.float32) - ny // 2) * np.float32(self.dx_m)
        self.X, self.Y = np.meshgrid(x, y)
        self.R2 = self.X * self.X + self.Y * self.Y
        self.reference_radius_m = self.diameter_m / 2.0
        self.normalized_x = self.X / np.float32(self.reference_radius_m)
        self.normalized_y = self.Y / np.float32(self.reference_radius_m)
        self.normalized_r2 = (
            self.normalized_x * self.normalized_x
            + self.normalized_y * self.normalized_y
        )
        self.phase_factor = np.exp(1j * self.phase).astype(np.complex64)
        self.base_amplitude = make_input_gaussian(
            shape=self.phase.shape,
            dx_doe_m=self.dx_m,
            gaussian_1e2_diameter_m=self.diameter_m,
            clear_aperture_m=self.clear_aperture_m,
            xp=np,
            dtype=np.float32,
        )
        self.hg20 = self._orthogonal_mode("x")
        self.hg02 = self._orthogonal_mode("y")

    def _orthogonal_mode(self, axis: str) -> np.ndarray:
        coordinate = self.X if axis == "x" else self.Y
        raw = self.base_amplitude * (
            2.0 * (coordinate / np.float32(self.reference_radius_m)) ** 2 - 0.5
        )
        projection = float(np.sum(raw * self.base_amplitude, dtype=np.float64))
        raw = raw - projection * self.base_amplitude
        return normalize_power(raw.astype(np.float32))

    def gaussian_amplitude(
        self,
        diameter_x_mm: float,
        diameter_y_mm: float,
        offset_x_mm: float,
        offset_y_mm: float,
    ) -> np.ndarray:
        wx = np.float32(diameter_x_mm * 0.5e-3)
        wy = np.float32(diameter_y_mm * 0.5e-3)
        ox = np.float32(offset_x_mm * 1e-3)
        oy = np.float32(offset_y_mm * 1e-3)
        amp = np.exp(
            -((self.X - ox) ** 2 / (wx * wx) + (self.Y - oy) ** 2 / (wy * wy))
        ).astype(np.float32)
        amp[self.R2 > np.float32(self.clear_aperture_m / 2.0) ** 2] = 0.0
        return normalize_power(amp)

    def reconstruct(self, params: dict[str, Any]) -> np.ndarray:
        diameter_x_mm = float(params.get("diameter_x_mm", self.diameter_m * 1e3))
        diameter_y_mm = float(params.get("diameter_y_mm", self.diameter_m * 1e3))
        offset_x_mm = float(params.get("beam_offset_x_mm", 0.0))
        offset_y_mm = float(params.get("beam_offset_y_mm", 0.0))
        if (
            diameter_x_mm != self.diameter_m * 1e3
            or diameter_y_mm != self.diameter_m * 1e3
            or offset_x_mm != 0.0
            or offset_y_mm != 0.0
        ):
            amplitude = self.gaussian_amplitude(
                diameter_x_mm, diameter_y_mm, offset_x_mm, offset_y_mm
            )
        else:
            amplitude = self.base_amplitude.copy()

        dip_family = params.get("amplitude_dip_family")
        if dip_family:
            depth = np.float32(params["amplitude_dip_depth"])
            sigma = np.float32(
                params.get(
                    "amplitude_dip_sigma_mm",
                    params.get("amplitude_dip_sigma_x_mm", 1.0),
                )
                * 1e-3
            )
            dip_offset_x = np.float32(
                float(params.get("amplitude_dip_offset_x_mm", 0.0)) * 1e-3
            )
            dip_offset_y = np.float32(
                float(params.get("amplitude_dip_offset_y_mm", 0.0)) * 1e-3
            )
            if dip_family == "center_dip":
                exponent = -0.5 * (
                    (self.X - dip_offset_x) ** 2
                    + (self.Y - dip_offset_y) ** 2
                ) / (sigma * sigma)
            elif dip_family == "elliptical_input_dip":
                sigma_x = np.float32(params["amplitude_dip_sigma_x_mm"] * 1e-3)
                sigma_y = np.float32(params["amplitude_dip_sigma_y_mm"] * 1e-3)
                exponent = -0.5 * (
                    (self.X - dip_offset_x) ** 2 / (sigma_x * sigma_x)
                    + (self.Y - dip_offset_y) ** 2 / (sigma_y * sigma_y)
                )
            elif dip_family == "horizontal_input_dip":
                exponent = (
                    -0.5
                    * (self.Y - dip_offset_y)
                    * (self.Y - dip_offset_y)
                    / (sigma * sigma)
                )
            elif dip_family == "vertical_input_dip":
                exponent = (
                    -0.5
                    * (self.X - dip_offset_x)
                    * (self.X - dip_offset_x)
                    / (sigma * sigma)
                )
            else:
                raise ValueError(f"Unsupported amplitude dip family: {dip_family}")
            amplitude *= 1.0 - depth * np.exp(exponent).astype(np.float32)

        aperture_diameter_mm = params.get("aperture_diameter_mm")
        if aperture_diameter_mm is not None:
            aperture_radius = np.float32(float(aperture_diameter_mm) * 0.5e-3)
            aperture_x = np.float32(float(params.get("aperture_offset_x_mm", 0.0)) * 1e-3)
            aperture_y = np.float32(float(params.get("aperture_offset_y_mm", 0.0)) * 1e-3)
            aperture = (
                (self.X - aperture_x) ** 2 + (self.Y - aperture_y) ** 2
                <= aperture_radius * aperture_radius
            )
            amplitude *= aperture

        hg20 = float(params.get("hg20_coefficient", 0.0))
        hg02 = float(params.get("hg02_coefficient", 0.0))
        if hg20 != 0.0 or hg02 != 0.0:
            amplitude = (
                amplitude.astype(np.float32)
                + np.float32(hg20) * self.hg20
                + np.float32(hg02) * self.hg02
            )

        amplitude = normalize_power(amplitude.astype(np.float32))
        defocus = float(params.get("defocus_waves", 0.0))
        astigmatism = float(params.get("astigmatism_waves", 0.0))
        spherical = float(params.get("spherical_waves", 0.0))
        if defocus != 0.0 or astigmatism != 0.0 or spherical != 0.0:
            phase_waves = (
                defocus * self.normalized_r2
                + astigmatism
                * (
                    self.normalized_x * self.normalized_x
                    - self.normalized_y * self.normalized_y
                )
                + spherical * self.normalized_r2 * self.normalized_r2
            )
            aberration = np.exp(1j * np.float32(2.0 * np.pi) * phase_waves).astype(
                np.complex64
            )
            field = amplitude * self.phase_factor * aberration
        else:
            field = amplitude * self.phase_factor
        return intensity(forward_fft(field.astype(np.complex64), np), np).astype(
            np.float32
        )


def evaluate_image(
    image: np.ndarray,
    simulator: Simulator,
    experiment: dict[str, Any],
) -> tuple[np.ndarray, dict[str, Any], dict[str, float]]:
    """Normalize one simulated focal image and compare it with experiment."""
    x_um = simulator.base.x_um
    y_um = simulator.base.y_um
    core_width_x_um, core_width_y_um = experiment["core_width_um"]
    center_width_x_um, center_width_y_um = experiment["center_window_um"]
    x_core = np.abs(x_um) <= core_width_x_um / 2.0
    y_core = np.abs(y_um) <= core_width_y_um / 2.0
    core_mean = float(np.mean(image[np.ix_(y_core, x_core)]))
    normalized = image.astype(np.float64) / core_mean
    metrics = morphology_metrics(
        normalized,
        x_um,
        y_um,
        core_width_x_um,
        core_width_y_um,
        center_width_x_um,
        center_width_y_um,
    )
    score = candidate_score(metrics, experiment["metrics"])
    return normalized, metrics, score


def plot_candidate_overview(
    output_path: Path,
    experiment: dict[str, Any],
    simulator: Simulator,
    candidates: list[dict[str, Any]],
) -> None:
    """Plot experiment, baseline, and the best simulated candidates."""
    selected = candidates[:5]
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]] = [
        (
            "Experiment: average of 16 frames",
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            experiment["metrics"],
        )
    ]
    for candidate in selected:
        image = simulator.reconstruct(candidate["params"])
        normalized, metrics, _ = evaluate_image(image, simulator, experiment)
        panels.append(
            (
                f"{candidate['name']}\nscore={candidate['score']:.3f}",
                normalized,
                simulator.base.x_um,
                simulator.base.y_um,
                metrics,
            )
        )

    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    for ax, (title, image, x_um, y_um, metrics) in zip(axes.ravel(), panels):
        xmask = np.abs(x_um) <= 220.0
        ymask = np.abs(y_um) <= 100.0
        roi = image[np.ix_(ymask, xmask)]
        xa = x_um[xmask]
        ya = y_um[ymask]
        im = ax.imshow(
            roi,
            origin="upper",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            cmap="turbo",
            vmin=0.35,
            vmax=1.50,
            aspect="equal",
        )
        ax.set_title(
            f"{title}\ncenter={metrics['center_window_over_core_mean']:.3f}, "
            f"middle/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"RMS={100.0 * metrics['core_rms_fraction']:.1f}%"
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("Fixed V2 phase: experimental center depression versus scan candidates")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_profile_comparison(
    output_path: Path,
    experiment: dict[str, Any],
    simulator: Simulator,
    candidates: list[dict[str, Any]],
) -> None:
    """Compare common-core profiles for experiment, baseline, and best cases."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    target = experiment["metrics"]
    unit = np.asarray(target["unit_profile_axis"])
    axes[0].plot(unit, target["unit_profile_x"], lw=2.2, label="experiment")
    axes[1].plot(unit, target["unit_profile_y"], lw=2.2, label="experiment")
    for candidate in candidates[:5]:
        image = simulator.reconstruct(candidate["params"])
        _, metrics, _ = evaluate_image(image, simulator, experiment)
        axes[0].plot(unit, metrics["unit_profile_x"], label=candidate["name"])
        axes[1].plot(unit, metrics["unit_profile_y"], label=candidate["name"])
    for ax, axis_name in zip(axes, ("X", "Y")):
        ax.axhline(1.0, color="gray", ls=":")
        ax.set_title(f"{axis_name} center-band profile")
        ax.set_xlabel("normalized coordinate inside experimental size50 box")
        ax.set_ylabel("I / mean(core)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_scan_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "rank",
        "score",
        "family",
        "name",
        "center_window_over_core_mean",
        "middle_third_over_side_thirds",
        "core_rms_fraction",
        "top_over_bottom",
        "left_over_right",
        "profile_rmse_x",
        "profile_rmse_y",
        "params_json",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            metrics = row["metrics"]
            writer.writerow(
                {
                    "rank": rank,
                    "score": row["score"],
                    "family": row["family"],
                    "name": row["name"],
                    "center_window_over_core_mean": metrics[
                        "center_window_over_core_mean"
                    ],
                    "middle_third_over_side_thirds": metrics[
                        "middle_third_over_side_thirds"
                    ],
                    "core_rms_fraction": metrics["core_rms_fraction"],
                    "top_over_bottom": metrics["top_over_bottom"],
                    "left_over_right": metrics["left_over_right"],
                    "profile_rmse_x": row["score_details"]["profile_rmse_x"],
                    "profile_rmse_y": row["score_details"]["profile_rmse_y"],
                    "params_json": json.dumps(row["params"], ensure_ascii=False),
                }
            )


def write_summary(
    path: Path,
    case_dir: Path,
    experiment: dict[str, Any],
    baseline_validation: dict[str, Any],
    rows: list[dict[str, Any]],
    elapsed_seconds: float,
) -> None:
    best_by_family: dict[str, dict[str, Any]] = {}
    for row in rows:
        best_by_family.setdefault(row["family"], row)
    best = rows[0]
    target = experiment["metrics"]
    lines = [
        "# V2 平顶光中心空洞第一轮正向传播扫描",
        "",
        f"- V2 相位：`{case_dir / 'phase_refined.npy'}`",
        "- 正向传播：复用项目 `src.propagation.forward_fft`（居中、正交 FFT）。",
        f"- 实验目标：`20260716-1.bmData` 的 {experiment['frame_count']} 帧平均。",
        f"- 扫描耗时：{elapsed_seconds:.1f} s。",
        "",
        "## 基准复现",
        "",
        f"- 与保存的 `reconstruction_refined.npy` 相对 L2 误差："
        f"`{baseline_validation['relative_l2_error']:.6g}`。",
        f"- 基准中心窗口/核心均值："
        f"`{baseline_validation['baseline_metrics']['center_window_over_core_mean']:.4f}`。",
        "",
        "## 实验目标",
        "",
        f"- 物理核心框：`{experiment['core_width_um'][0]:.2f} × "
        f"{experiment['core_width_um'][1]:.2f} um`。",
        f"- 中心窗口：`{experiment['center_window_um'][0]:.2f} × "
        f"{experiment['center_window_um'][1]:.2f} um`。",
        f"- 中心窗口/核心均值：`{target['center_window_over_core_mean']:.4f}`。",
        f"- 中间三分之一区域/左右区域：`{target['middle_third_over_side_thirds']:.4f}`。",
        f"- 核心 RMS：`{100.0 * target['core_rms_fraction']:.2f}%`。",
        "",
        "## 当前最佳候选",
        "",
        f"- 名称：`{best['name']}`",
        f"- 参数：`{json.dumps(best['params'], ensure_ascii=False)}`",
        f"- 综合差异分数：`{best['score']:.4f}`（越低越接近实验）。",
        f"- 中心窗口/核心均值："
        f"`{best['metrics']['center_window_over_core_mean']:.4f}`。",
        f"- 中间三分之一区域/左右区域："
        f"`{best['metrics']['middle_third_over_side_thirds']:.4f}`。",
        f"- 核心 RMS：`{100.0 * best['metrics']['core_rms_fraction']:.2f}%`。",
        "",
        "## 每类物理因素的最佳结果",
        "",
        "| 因素 | 最佳候选 | 分数 | 中心比 | 中部/两侧 | RMS |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for family, row in sorted(
        best_by_family.items(), key=lambda item: item[1]["score"]
    ):
        lines.append(
            f"| {family} | `{row['name']}` | {row['score']:.3f} | "
            f"{row['metrics']['center_window_over_core_mean']:.3f} | "
            f"{row['metrics']['middle_third_over_side_thirds']:.3f} | "
            f"{100.0 * row['metrics']['core_rms_fraction']:.1f}% |"
        )
    lines.extend(
        [
            "",
            "## 说明",
            "",
            "- 这一轮固定 V2 相位，只改变入射振幅、孔径或附加波前。",
            "- 闪耀光栅只平移焦斑，不影响中心形态，因此没有加入。",
            "- 安装补偿等价于相位和入射光的相对横向偏移，本轮通过 beam offset 扫描。",
            "- 第一轮结果用于锁定能够产生同类中心暗带的因素，不能单凭拟合分数认定真实物理原因。",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        args.bmdata, args.percent_summary, args.hollow_summary
    )
    simulator = Simulator(args.case_dir)

    stored = np.load(args.case_dir / "reconstruction_refined.npy").astype(np.float64)
    baseline_image = simulator.reconstruct({})
    baseline_normalized, baseline_metrics, baseline_score = evaluate_image(
        baseline_image, simulator, experiment
    )
    relative_l2 = float(
        np.linalg.norm(baseline_image.astype(np.float64) - stored)
        / np.linalg.norm(stored)
    )
    stored_scale = float(np.sum(baseline_image, dtype=np.float64) / np.sum(stored))
    baseline_validation = {
        "relative_l2_error": relative_l2,
        "sum_ratio_recomputed_over_stored": stored_scale,
        "baseline_metrics": baseline_metrics,
        "baseline_score": baseline_score,
    }
    (args.outdir / "baseline_validation.json").write_text(
        json.dumps(baseline_validation, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    cases = build_scan_cases()
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    for index, case in enumerate(cases, start=1):
        image = baseline_image if case["family"] == "baseline" else simulator.reconstruct(case["params"])
        _, metrics, score_details = evaluate_image(image, simulator, experiment)
        rows.append(
            {
                **case,
                "score": score_details["score"],
                "metrics": metrics,
                "score_details": score_details,
            }
        )
        if index == 1 or index % 40 == 0 or index == len(cases):
            print(
                f"[{index:4d}/{len(cases)}] {case['name']}: "
                f"score={score_details['score']:.3f}",
                flush=True,
            )
    elapsed = time.perf_counter() - started
    rows.sort(key=lambda row: row["score"])

    serializable = {
        "case_dir": str(args.case_dir),
        "phase_source": str(args.case_dir / "phase_refined.npy"),
        "bmdata": str(args.bmdata),
        "forward_propagator": str(PROJECT_ROOT / "src" / "propagation.py"),
        "elapsed_seconds": elapsed,
        "candidate_count": len(rows),
        "experiment": {
            key: value
            for key, value in experiment.items()
            if key not in {"average_signal", "normalized", "x_um", "y_um"}
        },
        "baseline_validation": baseline_validation,
        "candidates": rows,
    }
    (args.outdir / "scan_results.json").write_text(
        json.dumps(serializable, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_scan_csv(args.outdir / "scan_results.csv", rows)

    best = rows[0]
    best_image = simulator.reconstruct(best["params"])
    np.save(args.outdir / "best_candidate_intensity.npy", best_image)
    (args.outdir / "best_candidate.json").write_text(
        json.dumps(best, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_candidate_overview(
        args.outdir / "candidate_overview.png", experiment, simulator, rows
    )
    plot_profile_comparison(
        args.outdir / "profile_comparison.png", experiment, simulator, rows
    )
    write_summary(
        args.outdir / "SUMMARY.md",
        args.case_dir,
        experiment,
        baseline_validation,
        rows,
        elapsed,
    )
    print(
        json.dumps(
            {
                "candidate_count": len(rows),
                "elapsed_seconds": elapsed,
                "best": {
                    "name": best["name"],
                    "family": best["family"],
                    "params": best["params"],
                    "score": best["score"],
                    "metrics": {
                        key: value
                        for key, value in best["metrics"].items()
                        if not key.startswith("unit_profile")
                    },
                },
                "outputs": {
                    "summary": str(args.outdir / "SUMMARY.md"),
                    "overview": str(args.outdir / "candidate_overview.png"),
                    "profiles": str(args.outdir / "profile_comparison.png"),
                },
            },
            ensure_ascii=False,
            indent=2,
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--bmdata", type=Path, default=DEFAULT_BMDATA)
    parser.add_argument(
        "--percent-summary", type=Path, default=DEFAULT_PERCENT_SUMMARY
    )
    parser.add_argument("--hollow-summary", type=Path, default=DEFAULT_HOLLOW_SUMMARY)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ANALYSIS_ROOT / "results" / "01_initial_scan",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
