"""Propagate the fixed V2 phase with the measured 2026-07-16 input intensity."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
from scipy.ndimage import map_coordinates

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulate_v2_hollow_scan import (
    ANALYSIS_ROOT,
    DEFAULT_BMDATA,
    DEFAULT_CASE_DIR,
    DEFAULT_HOLLOW_SUMMARY,
    DEFAULT_PERCENT_SUMMARY,
    REAL_TEST_ROOT,
    Simulator,
    evaluate_image,
    load_experimental_target,
    normalize_power,
)
from analyze_rect_flattop_size import load_spiricon_frame, robust_background
from src.propagation import forward_fft, intensity


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "03_measured_input_forward"
INPUT_BGDATA = REAL_TEST_ROOT / "20260716" / "G-光斑-1.bgData"
INPUT_SUMMARY = (
    REAL_TEST_ROOT
    / "20260716"
    / "analysis_G_spot_1"
    / "G-光斑-1_beam_size_summary.json"
)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        DEFAULT_BMDATA, DEFAULT_PERCENT_SUMMARY, DEFAULT_HOLLOW_SUMMARY
    )
    simulator = Simulator(DEFAULT_CASE_DIR)
    input_summary = json.loads(INPUT_SUMMARY.read_text(encoding="utf-8"))
    measured_image, _ = load_spiricon_frame(INPUT_BGDATA)
    background, background_sigma = robust_background(measured_image, corner_px=200)
    measured_signal = np.clip(measured_image - background, 0.0, None).astype(np.float32)

    x0, y0, width, height = (
        int(value) for value in input_summary["beam_region_bbox_px"]
    )
    component_signal = np.zeros_like(measured_signal)
    component_signal[y0 : y0 + height, x0 : x0 + width] = measured_signal[
        y0 : y0 + height, x0 : x0 + width
    ]
    center_x_px, center_y_px = (
        float(value) for value in input_summary["moments"]["center_px_global"]
    )
    source_sx_um = float(input_summary["pixel_scale_x_um"])
    source_sy_um = float(input_summary["pixel_scale_y_um"])
    source_x_px = center_x_px + simulator.X * np.float32(1e6 / source_sx_um)
    source_y_px = center_y_px + simulator.Y * np.float32(1e6 / source_sy_um)
    sampled_intensity = map_coordinates(
        component_signal,
        [source_y_px, source_x_px],
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32)
    sampled_amplitude = normalize_power(np.sqrt(sampled_intensity).astype(np.float32))
    measured_field = sampled_amplitude * simulator.phase_factor
    focal_intensity = intensity(
        forward_fft(measured_field.astype(np.complex64), np), np
    ).astype(np.float32)
    focal_normalized, focal_metrics, focal_score = evaluate_image(
        focal_intensity, simulator, experiment
    )
    baseline_intensity = simulator.reconstruct({})
    baseline_normalized, baseline_metrics, baseline_score = evaluate_image(
        baseline_intensity, simulator, experiment
    )

    result = {
        "phase_source": str(DEFAULT_CASE_DIR / "phase_refined.npy"),
        "measured_input_source": str(INPUT_BGDATA),
        "measured_input_summary": str(INPUT_SUMMARY),
        "assumption": (
            "Use sqrt(background-subtracted measured intensity) as the DOE-plane "
            "amplitude and assume a flat input phase. The selected Gaussian component "
            "is centered on the V2 phase."
        ),
        "source_background": background,
        "source_background_robust_sigma": background_sigma,
        "source_center_px": [center_x_px, center_y_px],
        "source_component_bbox_px": [x0, y0, width, height],
        "source_pixel_scale_um": [source_sx_um, source_sy_um],
        "experiment_metrics": {
            key: value
            for key, value in experiment["metrics"].items()
            if not key.startswith("unit_profile")
        },
        "baseline_metrics": {
            key: value
            for key, value in baseline_metrics.items()
            if not key.startswith("unit_profile")
        },
        "baseline_score": baseline_score,
        "measured_amplitude_forward_metrics": {
            key: value
            for key, value in focal_metrics.items()
            if not key.startswith("unit_profile")
        },
        "measured_amplitude_forward_score": focal_score,
    }
    (OUTPUT_DIR / "measured_input_forward_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    np.save(OUTPUT_DIR / "measured_input_resampled_intensity.npy", sampled_intensity)
    np.save(OUTPUT_DIR / "measured_input_focal_intensity.npy", focal_intensity)

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    source_x_mm = (
        np.arange(measured_signal.shape[1], dtype=np.float64) - center_x_px
    ) * source_sx_um * 1e-3
    source_y_mm = (
        np.arange(measured_signal.shape[0], dtype=np.float64) - center_y_px
    ) * source_sy_um * 1e-3
    sx_mask = np.abs(source_x_mm) <= 7.5
    sy_mask = np.abs(source_y_mm) <= 7.5
    source_roi = component_signal[np.ix_(sy_mask, sx_mask)]
    source_roi /= max(float(np.max(source_roi)), 1e-20)
    im = axes[0, 0].imshow(
        source_roi,
        origin="upper",
        extent=[-7.5, 7.5, 7.5, -7.5],
        cmap="turbo",
        vmin=0,
        vmax=1,
        aspect="equal",
    )
    axes[0, 0].set_title("Measured input component")
    axes[0, 0].set_xlabel("x / mm")
    axes[0, 0].set_ylabel("y / mm")
    fig.colorbar(im, ax=axes[0, 0], label="I / peak")

    doe_x_mm = simulator.X[0, :].astype(np.float64) * 1e3
    doe_y_mm = simulator.Y[:, 0].astype(np.float64) * 1e3
    dx_mask = np.abs(doe_x_mm) <= 7.5
    dy_mask = np.abs(doe_y_mm) <= 7.5
    doe_roi = sampled_intensity[np.ix_(dy_mask, dx_mask)].astype(np.float64)
    doe_roi /= max(float(np.max(doe_roi)), 1e-20)
    im = axes[0, 1].imshow(
        doe_roi,
        origin="upper",
        extent=[-7.5, 7.5, 7.5, -7.5],
        cmap="turbo",
        vmin=0,
        vmax=1,
        aspect="equal",
    )
    axes[0, 1].set_title("Measured intensity resampled on V2 DOE grid")
    axes[0, 1].set_xlabel("x / mm")
    axes[0, 1].set_ylabel("y / mm")
    fig.colorbar(im, ax=axes[0, 1], label="I / peak")

    panels = (
        (
            axes[1, 0],
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            "Experiment: average of 16 frames",
            experiment["metrics"],
        ),
        (
            axes[1, 1],
            focal_normalized,
            simulator.base.x_um,
            simulator.base.y_um,
            "V2 with measured input amplitude, flat wavefront",
            focal_metrics,
        ),
    )
    for ax, image, x_um, y_um, title, metrics in panels:
        xmask = np.abs(x_um) <= 220
        ymask = np.abs(y_um) <= 100
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
        fig.colorbar(im, ax=ax, label="I / mean(core)")
    figure_path = OUTPUT_DIR / "measured_input_forward_diagnostics.png"
    fig.suptitle("Fixed V2 phase driven by the measured 2026-07-16 input amplitude")
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)
    result["output_figure"] = str(figure_path)
    (OUTPUT_DIR / "measured_input_forward_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
