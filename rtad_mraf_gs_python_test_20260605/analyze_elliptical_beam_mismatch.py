"""Compare one refined DOE phase under two elliptical Gaussian inputs."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_simulation_percent_energy import percent_energy_bbox
from src.diagnostics import _load_phase_from_case, compute_diagnostics, load_case_data
from src.propagation import forward_fft, intensity, make_input_gaussian


def reconstruct(
    phase: np.ndarray,
    config: dict[str, Any],
    diameter_x_mm: float,
    diameter_y_mm: float,
) -> np.ndarray:
    """Forward propagate the fixed phase with an elliptical Gaussian input."""
    amplitude = make_input_gaussian(
        shape=phase.shape,
        dx_doe_m=float(config["grid"]["dx_doe_m"]),
        gaussian_1e2_diameter_m=None,
        gaussian_1e2_diameter_x_m=diameter_x_mm * 1e-3,
        gaussian_1e2_diameter_y_m=diameter_y_mm * 1e-3,
        clear_aperture_m=float(config["physical"]["clear_aperture_m"]),
        xp=np,
        dtype=np.float32,
    )
    field = amplitude * np.exp(1j * phase).astype(np.complex64)
    return intensity(forward_fft(field, np), np).astype(np.float64)


def summarize(base: Any, image: np.ndarray, diameter_x_mm: float, diameter_y_mm: float) -> dict[str, Any]:
    """Compute the standard diagnostics and center/flat statistics."""
    source = f"recomputed with elliptical Gaussian X={diameter_x_mm:.3f} mm, Y={diameter_y_mm:.3f} mm"
    data = replace(base, intensity_raw=image, intensity_source=source)
    metrics, _, warnings = compute_diagnostics(data)
    energy = percent_energy_bbox(image, base.x_um, base.y_um)
    normalized = image / float(np.mean(image[base.mask_flat]))
    cx = int(np.argmin(np.abs(base.x_um)))
    cy = int(np.argmin(np.abs(base.y_um)))
    center5 = normalized[cy - 2 : cy + 3, cx - 2 : cx + 3]
    center15 = normalized[cy - 7 : cy + 8, cx - 7 : cx + 8]
    flat_values = normalized[base.mask_flat]
    return {
        "input_gaussian_1e2_diameter_x_mm": diameter_x_mm,
        "input_gaussian_1e2_diameter_y_mm": diameter_y_mm,
        "size50_um": [metrics["size50_x_um"], metrics["size50_y_um"]],
        "size90_um": [metrics["size90_x_um"], metrics["size90_y_um"]],
        "size13p5_um": [metrics["size13p5_x_um"], metrics["size13p5_y_um"]],
        "rms_nonuniformity_percent": metrics["rms_nonuniformity_percent"],
        "uniformity_rms_percent": metrics["uniformity_rms_percent"],
        "efficiency_e2_percent": metrics["efficiency_e2_percent"],
        "percent_energy_86p5_span_um": energy["bbox_center_span_um"],
        "center_5x5_mean_over_flat_mean": float(np.mean(center5)),
        "center_15x15_mean_over_flat_mean": float(np.mean(center15)),
        "flat_min_over_flat_mean": float(np.min(flat_values)),
        "flat_max_over_flat_mean": float(np.max(flat_values)),
        "warnings": warnings,
    }


def analyze(
    case_dir: Path,
    reference_x_mm: float,
    reference_y_mm: float,
    test_x_mm: float,
    test_y_mm: float,
) -> dict[str, Any]:
    """Run the fixed-phase mismatch comparison and save JSON plus a diagnostic figure."""
    case_dir = case_dir.resolve()
    base = load_case_data(case_dir)
    phase = _load_phase_from_case(case_dir)
    reference = reconstruct(phase, base.config, reference_x_mm, reference_y_mm)
    test = reconstruct(phase, base.config, test_x_mm, test_y_mm)
    ref_summary = summarize(base, reference, reference_x_mm, reference_y_mm)
    test_summary = summarize(base, test, test_x_mm, test_y_mm)

    name = (
        f"input_beam_X{reference_x_mm:g}_Y{reference_y_mm:g}_vs_"
        f"X{test_x_mm:g}_Y{test_y_mm:g}_analysis"
    ).replace(".", "p")
    outdir = case_dir / name
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "elliptical_beam_mismatch_summary.json"
    figure_path = outdir / "elliptical_beam_mismatch_diagnostics.png"

    delta = {
        "size50_um": (np.asarray(test_summary["size50_um"]) - np.asarray(ref_summary["size50_um"])).tolist(),
        "size90_um": (np.asarray(test_summary["size90_um"]) - np.asarray(ref_summary["size90_um"])).tolist(),
        "size13p5_um": (np.asarray(test_summary["size13p5_um"]) - np.asarray(ref_summary["size13p5_um"])).tolist(),
        "rms_nonuniformity_percentage_points": test_summary["rms_nonuniformity_percent"] - ref_summary["rms_nonuniformity_percent"],
        "efficiency_e2_percentage_points": test_summary["efficiency_e2_percent"] - ref_summary["efficiency_e2_percent"],
        "center_5x5_mean_over_flat_mean": test_summary["center_5x5_mean_over_flat_mean"] - ref_summary["center_5x5_mean_over_flat_mean"],
        "center_15x15_mean_over_flat_mean": test_summary["center_15x15_mean_over_flat_mean"] - ref_summary["center_15x15_mean_over_flat_mean"],
    }
    result = {
        "case_dir": str(case_dir),
        "phase_source": str(case_dir / "phase_refined.npy"),
        "note": "The 6.5 x 6.5 mm V2 refined phase is fixed. Only the modeled input Gaussian changes; no installation shift or blaze is applied.",
        "reference": ref_summary,
        "test": test_summary,
        "delta_test_minus_reference": delta,
        "outputs": {"summary_json": str(json_path), "diagnostic_png": str(figure_path)},
    }
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    xmask = np.abs(base.x_um) <= 240.0
    ymask = np.abs(base.y_um) <= 110.0
    xa = base.x_um[xmask]
    ya = base.y_um[ymask]
    ref_norm = reference / float(np.mean(reference[base.mask_flat]))
    test_norm = test / float(np.mean(test[base.mask_flat]))
    diff_percent = 100.0 * (test_norm - ref_norm)
    fig, axes = plt.subplots(2, 3, figsize=(18, 9), constrained_layout=True)
    cases = [
        (ref_norm, ref_summary, f"Reference X={reference_x_mm:g}, Y={reference_y_mm:g} mm"),
        (test_norm, test_summary, f"Test X={test_x_mm:g}, Y={test_y_mm:g} mm"),
    ]
    for column, (image, summary, title) in enumerate(cases):
        roi = image[np.ix_(ymask, xmask)]
        im = axes[0, column].imshow(
            roi,
            origin="upper",
            cmap="turbo",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            vmin=0.85,
            vmax=1.15,
            aspect="equal",
        )
        axes[0, column].set_title(
            f"{title}\ncenter={summary['center_5x5_mean_over_flat_mean']:.4f}, "
            f"RMS={summary['rms_nonuniformity_percent']:.3f}%"
        )
        axes[0, column].set_xlabel("x / um")
        axes[0, column].set_ylabel("y / um")
        fig.colorbar(im, ax=axes[0, column], label="I / mean(flat)")

    diff_roi = diff_percent[np.ix_(ymask, xmask)]
    limit = max(1.0, float(np.percentile(np.abs(diff_roi), 99)))
    im = axes[0, 2].imshow(
        diff_roi,
        origin="upper",
        cmap="coolwarm",
        extent=[xa[0], xa[-1], ya[-1], ya[0]],
        vmin=-limit,
        vmax=limit,
        aspect="equal",
    )
    axes[0, 2].set_title("Test - reference\npercentage points of normalized intensity")
    axes[0, 2].set_xlabel("x / um")
    axes[0, 2].set_ylabel("y / um")
    fig.colorbar(im, ax=axes[0, 2], label="delta I / mean(flat), %")

    cx = int(np.argmin(np.abs(base.x_um)))
    cy = int(np.argmin(np.abs(base.y_um)))
    band = 2
    for image, label in ((ref_norm, "6.5 x 6.5"), (test_norm, "6.5 x 6.3")):
        axes[1, 0].plot(base.x_um, image[cy - band : cy + band + 1, :].mean(axis=0), label=label)
        axes[1, 1].plot(base.y_um, image[:, cx - band : cx + band + 1].mean(axis=1), label=label)
    axes[1, 0].set_xlim(-240, 240)
    axes[1, 0].set_title("X center profile")
    axes[1, 0].set_xlabel("x / um")
    axes[1, 0].set_ylabel("I / mean(flat)")
    axes[1, 1].set_xlim(-110, 110)
    axes[1, 1].set_title("Y center profile")
    axes[1, 1].set_xlabel("y / um")
    axes[1, 1].set_ylabel("I / mean(flat)")
    for ax in axes[1, :2]:
        ax.axhline(1.0, color="gray", ls=":")
        ax.grid(True, alpha=0.25)
        ax.legend()

    axes[1, 2].axis("off")
    lines = [
        "Test minus reference",
        f"size50 X/Y: {delta['size50_um'][0]:+.3f} / {delta['size50_um'][1]:+.3f} um",
        f"size90 X/Y: {delta['size90_um'][0]:+.3f} / {delta['size90_um'][1]:+.3f} um",
        f"RMS: {delta['rms_nonuniformity_percentage_points']:+.3f} percentage points",
        f"e^-2 efficiency: {delta['efficiency_e2_percentage_points']:+.3f} percentage points",
        f"center 5x5: {delta['center_5x5_mean_over_flat_mean']:+.4f}",
        f"center 15x15: {delta['center_15x15_mean_over_flat_mean']:+.4f}",
    ]
    axes[1, 2].text(0.02, 0.95, "\n".join(lines), va="top", family="monospace", fontsize=12)
    fig.suptitle(f"Fixed 6.5 x 6.5 V2 phase under elliptical input mismatch\n{case_dir.name}")
    fig.savefig(figure_path, dpi=190)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--reference-x-mm", type=float, default=6.5)
    parser.add_argument("--reference-y-mm", type=float, default=6.5)
    parser.add_argument("--test-x-mm", type=float, default=6.5)
    parser.add_argument("--test-y-mm", type=float, default=6.3)
    args = parser.parse_args()
    result = analyze(
        args.case_dir,
        args.reference_x_mm,
        args.reference_y_mm,
        args.test_x_mm,
        args.test_y_mm,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
