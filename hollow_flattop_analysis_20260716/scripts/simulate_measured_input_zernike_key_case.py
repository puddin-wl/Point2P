"""Apply the key Z40+Z20 case to the measured Gaussian input amplitude."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scan_measured_input_wavefront import (
    INPUT_BGDATA,
    INPUT_SUMMARY,
    measured_amplitude_on_doe,
)
from simulate_v2_hollow_scan import (
    ANALYSIS_ROOT,
    DEFAULT_BMDATA,
    DEFAULT_CASE_DIR,
    DEFAULT_HOLLOW_SUMMARY,
    DEFAULT_PERCENT_SUMMARY,
    Simulator,
    band_profile,
    evaluate_image,
    load_experimental_target,
)
from src.propagation import forward_fft, intensity


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "08_measured_input_zernike_key_case"
Z40_RMS_WAVES = -0.10625
Z20_RMS_WAVES = -0.25000


def compact_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in metrics.items()
        if not key.startswith("unit_profile")
    }


def render_comparison(
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, dict[str, Any], dict[str, float] | None]],
    experiment: dict[str, Any],
    output_path: Path,
    colors: list[str],
    figsize: tuple[float, float],
    suptitle: str,
) -> None:
    """Render image maps and matched band profiles for selected panels."""
    fig, axes = plt.subplots(
        2,
        len(panels),
        figsize=figsize,
        constrained_layout=True,
        squeeze=False,
    )
    for column, (title, image, px_um, py_um, metrics, footprint) in enumerate(
        panels
    ):
        ax = axes[0, column]
        xmask = np.abs(px_um) <= 230.0
        ymask = np.abs(py_um) <= 105.0
        roi = image[np.ix_(ymask, xmask)]
        xa = px_um[xmask]
        ya = py_um[ymask]
        im = ax.imshow(
            roi,
            origin="upper",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            cmap="turbo",
            vmin=0.30,
            vmax=1.55,
            aspect="equal",
        )
        footprint_text = (
            ""
            if footprint is None
            else f"\nIoU={footprint['iou_vs_ideal_v2_baseline']:.3f}"
        )
        ax.set_title(
            f"{title}\ncenter={metrics['center_window_over_core_mean']:.3f}, "
            f"middle/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"\nRMS={100.0 * metrics['core_rms_fraction']:.1f}%"
            f"{footprint_text}",
            fontsize=12,
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

        profile_ax = axes[1, column]
        profile_x_axis, profile_x = band_profile(
            image,
            px_um,
            py_um,
            "x",
            experiment["center_window_um"][1],
        )
        profile_y_axis, profile_y = band_profile(
            image,
            px_um,
            py_um,
            "y",
            experiment["center_window_um"][0],
        )
        profile_ax.plot(
            profile_x_axis,
            profile_x,
            color=colors[column],
            label="x profile",
        )
        profile_ax.plot(
            profile_y_axis,
            profile_y,
            color=colors[column],
            ls="--",
            label="y profile",
        )
        profile_ax.set_xlim(-210.0, 210.0)
        profile_ax.set_ylim(0.0, 1.65)
        profile_ax.set_xlabel("position / um")
        profile_ax.set_ylabel("normalized intensity")
        profile_ax.grid(True, alpha=0.25)
        profile_ax.legend()
    fig.suptitle(suptitle, fontsize=16)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        DEFAULT_BMDATA, DEFAULT_PERCENT_SUMMARY, DEFAULT_HOLLOW_SUMMARY
    )
    simulator = Simulator(DEFAULT_CASE_DIR)
    measured_amplitude = measured_amplitude_on_doe(simulator)

    pupil_radius_m = float(simulator.clear_aperture_m / 2.0)
    rho2 = simulator.R2 / np.float32(pupil_radius_m * pupil_radius_m)
    z20_map = np.float32(math.sqrt(3.0)) * (2.0 * rho2 - 1.0)
    z40_map = np.float32(math.sqrt(5.0)) * (
        6.0 * rho2 * rho2 - 6.0 * rho2 + 1.0
    )
    phase_waves = (
        np.float32(Z40_RMS_WAVES) * z40_map
        + np.float32(Z20_RMS_WAVES) * z20_map
    )
    aberration = np.exp(
        1j * np.float32(2.0 * np.pi) * phase_waves
    ).astype(np.complex64)

    def propagate(amplitude: np.ndarray, use_aberration: bool) -> np.ndarray:
        field = amplitude * simulator.phase_factor
        if use_aberration:
            field = field * aberration
        return intensity(
            forward_fft(field.astype(np.complex64), np), np
        ).astype(np.float32)

    ideal_key_raw = propagate(simulator.base_amplitude, True)
    measured_key_raw = propagate(measured_amplitude, True)
    measured_flat_raw = propagate(measured_amplitude, False)
    baseline_raw = simulator.reconstruct({})

    ideal_key, ideal_key_metrics, _ = evaluate_image(
        ideal_key_raw, simulator, experiment
    )
    measured_key, measured_key_metrics, _ = evaluate_image(
        measured_key_raw, simulator, experiment
    )
    measured_flat, measured_flat_metrics, _ = evaluate_image(
        measured_flat_raw, simulator, experiment
    )
    baseline, _, _ = evaluate_image(baseline_raw, simulator, experiment)

    x_um = simulator.base.x_um
    y_um = simulator.base.y_um
    X, Y = np.meshgrid(x_um, y_um)
    window = (np.abs(X) <= 240.0) & (np.abs(Y) <= 110.0)
    baseline_binary = (baseline >= 0.5) & window

    def footprint_metrics(normalized: np.ndarray) -> dict[str, float]:
        binary = (normalized >= 0.5) & window
        intersection = int(np.count_nonzero(binary & baseline_binary))
        union = int(np.count_nonzero(binary | baseline_binary))
        yy, xx = np.nonzero(binary)
        return {
            "iou_vs_ideal_v2_baseline": (
                float(intersection / union) if union else 0.0
            ),
            "size50_x_um": (
                float(x_um[xx.max()] - x_um[xx.min()])
                if xx.size
                else float("nan")
            ),
            "size50_y_um": (
                float(y_um[yy.max()] - y_um[yy.min()])
                if yy.size
                else float("nan")
            ),
        }

    ideal_footprint = footprint_metrics(ideal_key)
    measured_footprint = footprint_metrics(measured_key)
    result = {
        "phase_source": str(DEFAULT_CASE_DIR / "phase_refined.npy"),
        "measured_input_source": str(INPUT_BGDATA),
        "measured_input_summary": str(INPUT_SUMMARY),
        "zernike": {
            "Z40_primary_spherical_rms_waves": Z40_RMS_WAVES,
            "Z20_extra_defocus_rms_waves": Z20_RMS_WAVES,
            "pupil_diameter_mm": simulator.clear_aperture_m * 1e3,
            "Z20_definition": "sqrt(3)*(2*rho^2-1)",
            "Z40_definition": "sqrt(5)*(6*rho^4-6*rho^2+1)",
        },
        "experiment_metrics": compact_metrics(experiment["metrics"]),
        "ideal_gaussian_key_case": {
            "metrics": compact_metrics(ideal_key_metrics),
            "footprint": ideal_footprint,
        },
        "measured_gaussian_key_case": {
            "metrics": compact_metrics(measured_key_metrics),
            "footprint": measured_footprint,
        },
        "measured_gaussian_flat_wavefront": {
            "metrics": compact_metrics(measured_flat_metrics),
        },
    }
    np.save(OUTPUT_DIR / "measured_gaussian_zernike_focal_intensity.npy", measured_key_raw)
    (OUTPUT_DIR / "measured_gaussian_zernike_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    panels = [
        (
            "Experiment",
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            experiment["metrics"],
            None,
        ),
        (
            "Ideal Gaussian + key Z40/Z20",
            ideal_key,
            x_um,
            y_um,
            ideal_key_metrics,
            ideal_footprint,
        ),
        (
            "Measured Gaussian + same Z40/Z20",
            measured_key,
            x_um,
            y_um,
            measured_key_metrics,
            measured_footprint,
        ),
        (
            "Measured Gaussian + flat wavefront",
            measured_flat,
            x_um,
            y_um,
            measured_flat_metrics,
            footprint_metrics(measured_flat),
        ),
    ]
    parameter_title = (
        f"Z40={Z40_RMS_WAVES:+.5f}, Z20={Z20_RMS_WAVES:+.5f} RMS waves"
    )
    render_comparison(
        panels,
        experiment,
        OUTPUT_DIR / "MEASURED_GAUSSIAN_ZERNIKE_COMPARISON.png",
        ["black", "#1b9e77", "#d95f02", "#7570b3"],
        (22, 9),
        "Key Zernike case with measured 2026-07-16 Gaussian input\n"
        + parameter_title,
    )
    render_comparison(
        [panels[0], panels[2]],
        experiment,
        OUTPUT_DIR / "EXPERIMENT_VS_MEASURED_GAUSSIAN_ZERNIKE.png",
        ["black", "#d95f02"],
        (12, 9),
        "Experiment vs measured Gaussian with spherical aberration + defocus\n"
        + parameter_title,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
