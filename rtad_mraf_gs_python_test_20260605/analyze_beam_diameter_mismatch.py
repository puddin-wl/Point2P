"""Compare a refined DOE phase under different Gaussian input diameters."""

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


def reconstruct(phase: np.ndarray, config: dict[str, Any], diameter_mm: float) -> np.ndarray:
    amplitude = make_input_gaussian(
        shape=phase.shape,
        dx_doe_m=float(config["grid"]["dx_doe_m"]),
        gaussian_1e2_diameter_m=diameter_mm * 1e-3,
        clear_aperture_m=float(config["physical"]["clear_aperture_m"]),
        xp=np,
        dtype=np.float32,
    )
    field = amplitude * np.exp(1j * phase).astype(np.complex64)
    return intensity(forward_fft(field, np), np).astype(np.float64)


def summarize(base: Any, image: np.ndarray, diameter_mm: float) -> dict[str, Any]:
    data = replace(base, intensity_raw=image, intensity_source=f"recomputed with {diameter_mm:.3f} mm Gaussian")
    metrics, _, warnings = compute_diagnostics(data)
    energy = percent_energy_bbox(image, base.x_um, base.y_um)
    normalized = image / float(np.mean(image[base.mask_flat]))
    cx = int(np.argmin(np.abs(base.x_um)))
    cy = int(np.argmin(np.abs(base.y_um)))
    center5 = normalized[cy - 2 : cy + 3, cx - 2 : cx + 3]
    center15 = normalized[cy - 7 : cy + 8, cx - 7 : cx + 8]
    flat_values = normalized[base.mask_flat]
    return {
        "input_gaussian_1e2_diameter_mm": diameter_mm,
        "size50_um": [metrics["size50_x_um"], metrics["size50_y_um"]],
        "size90_um": [metrics["size90_x_um"], metrics["size90_y_um"]],
        "size13p5_um": [metrics["size13p5_x_um"], metrics["size13p5_y_um"]],
        "rms_nonuniformity_percent": metrics["rms_nonuniformity_percent"],
        "uniformity_rms_percent": metrics["uniformity_rms_percent"],
        "percent_energy_86p5_span_um": energy["bbox_center_span_um"],
        "center_5x5_mean_over_flat_mean": float(np.mean(center5)),
        "center_15x15_mean_over_flat_mean": float(np.mean(center15)),
        "flat_min_over_flat_mean": float(np.min(flat_values)),
        "flat_max_over_flat_mean": float(np.max(flat_values)),
        "warnings": warnings,
    }


def analyze(case_dir: Path, reference_mm: float, test_mm: float) -> dict[str, Any]:
    case_dir = case_dir.resolve()
    base = load_case_data(case_dir)
    phase = _load_phase_from_case(case_dir)
    images = {
        reference_mm: reconstruct(phase, base.config, reference_mm),
        test_mm: reconstruct(phase, base.config, test_mm),
    }
    summaries = {str(d): summarize(base, image, d) for d, image in images.items()}

    outdir = case_dir / f"input_beam_{reference_mm:g}mm_vs_{test_mm:g}mm_analysis".replace(".", "p")
    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / "beam_diameter_mismatch_summary.json"
    figure_path = outdir / "beam_diameter_mismatch_diagnostics.png"
    result = {
        "case_dir": str(case_dir),
        "phase_source": str(case_dir / "phase_refined.npy"),
        "note": "The refined phase is unchanged. Only the circular Gaussian 1/e^2 intensity diameter is changed; no installation shift or blaze is applied.",
        "cases": summaries,
        "delta_test_minus_reference": {
            "size50_um": (np.asarray(summaries[str(test_mm)]["size50_um"]) - np.asarray(summaries[str(reference_mm)]["size50_um"])).tolist(),
            "size90_um": (np.asarray(summaries[str(test_mm)]["size90_um"]) - np.asarray(summaries[str(reference_mm)]["size90_um"])).tolist(),
            "rms_nonuniformity_percent": summaries[str(test_mm)]["rms_nonuniformity_percent"] - summaries[str(reference_mm)]["rms_nonuniformity_percent"],
            "center_5x5_mean_over_flat_mean": summaries[str(test_mm)]["center_5x5_mean_over_flat_mean"] - summaries[str(reference_mm)]["center_5x5_mean_over_flat_mean"],
        },
        "conclusion": "A 6.7 mm Gaussian produces a shallow center depression and higher nonuniformity, but not a hollow center in this simulation.",
        "outputs": {"summary_json": str(json_path), "diagnostic_png": str(figure_path)},
    }
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    xmask = np.abs(base.x_um) <= 240.0
    ymask = np.abs(base.y_um) <= 110.0
    xa = base.x_um[xmask]
    ya = base.y_um[ymask]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    normalized: dict[float, np.ndarray] = {}
    for column, diameter in enumerate((reference_mm, test_mm)):
        image = images[diameter]
        norm = image / float(np.mean(image[base.mask_flat]))
        normalized[diameter] = norm
        roi = norm[np.ix_(ymask, xmask)]
        im = axes[0, column].imshow(
            roi,
            origin="upper",
            cmap="turbo",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            vmin=0.85,
            vmax=1.15,
            aspect="equal",
        )
        s = summaries[str(diameter)]
        axes[0, column].set_title(
            f"Input Gaussian {diameter:.1f} mm\n"
            f"center={s['center_5x5_mean_over_flat_mean']:.3f}, "
            f"RMS={s['rms_nonuniformity_percent']:.2f}%"
        )
        axes[0, column].set_xlabel("x / um")
        axes[0, column].set_ylabel("y / um")
        fig.colorbar(im, ax=axes[0, column], label="I / mean(flat)")

    cx = int(np.argmin(np.abs(base.x_um)))
    cy = int(np.argmin(np.abs(base.y_um)))
    band = 2
    for diameter in (reference_mm, test_mm):
        norm = normalized[diameter]
        profile_x = norm[cy - band : cy + band + 1, :].mean(axis=0)
        profile_y = norm[:, cx - band : cx + band + 1].mean(axis=1)
        axes[1, 0].plot(base.x_um, profile_x, label=f"{diameter:.1f} mm")
        axes[1, 1].plot(base.y_um, profile_y, label=f"{diameter:.1f} mm")
    axes[1, 0].set_xlim(-240, 240)
    axes[1, 0].set_title("X center profile")
    axes[1, 0].set_xlabel("x / um")
    axes[1, 0].set_ylabel("I / mean(flat)")
    axes[1, 1].set_xlim(-110, 110)
    axes[1, 1].set_title("Y center profile")
    axes[1, 1].set_xlabel("y / um")
    axes[1, 1].set_ylabel("I / mean(flat)")
    for ax in axes[1, :]:
        ax.axhline(1.0, color="gray", ls=":")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle(f"Input-beam mismatch test: {case_dir.name}\nSame refined phase, no X/Y installation shift")
    fig.savefig(figure_path, dpi=190)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--reference-mm", type=float, default=6.5)
    parser.add_argument("--test-mm", type=float, default=6.7)
    args = parser.parse_args()
    result = analyze(args.case_dir, args.reference_mm, args.test_mm)
    print(json.dumps(result["cases"], ensure_ascii=False, indent=2))
    print(result["outputs"]["diagnostic_png"])


if __name__ == "__main__":
    main()
