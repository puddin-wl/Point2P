"""Render the key local Z40+Z20 cases against experiment and ideal V2."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

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


RESULT_DIR = (
    ANALYSIS_ROOT / "results" / "07_zernike_spherical_defocus_local_refine"
)


def read_rows(path: Path) -> list[dict[str, Any]]:
    numeric = {
        "z40_rms_waves",
        "z20_rms_waves",
        "envelope_score",
        "footprint_iou",
        "outer_band_rmse",
        "inside_fill_fraction",
        "outside_fraction",
        "size50_x_um",
        "size50_y_um",
        "center_ratio",
        "middle_over_sides",
        "core_rms_fraction",
        "hole_target_score",
    }
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = {
                key: (float(value) if key in numeric else value)
                for key, value in raw.items()
            }
            row["rectangle_accepted"] = (
                str(raw["rectangle_accepted"]).lower() == "true"
            )
            rows.append(row)
    return rows


def main() -> None:
    experiment = load_experimental_target(
        DEFAULT_BMDATA, DEFAULT_PERCENT_SUMMARY, DEFAULT_HOLLOW_SUMMARY
    )
    simulator = Simulator(DEFAULT_CASE_DIR)
    target = experiment["metrics"]
    selected = [
        row
        for row in read_rows(RESULT_DIR / "best_defocus_for_each_spherical.csv")
        if row["rectangle_accepted"]
    ]
    all_accepted = [
        row
        for row in read_rows(RESULT_DIR / "all_cases.csv")
        if row["rectangle_accepted"]
    ]
    center_match = min(
        selected,
        key=lambda row: abs(
            row["center_ratio"] - target["center_window_over_core_mean"]
        ),
    )
    morphology_match = min(
        all_accepted, key=lambda row: row["hole_target_score"]
    )

    pupil_radius_m = float(simulator.clear_aperture_m / 2.0)
    rho2 = simulator.R2 / np.float32(pupil_radius_m * pupil_radius_m)
    z20_map = np.float32(math.sqrt(3.0)) * (2.0 * rho2 - 1.0)
    z40_map = np.float32(math.sqrt(5.0)) * (
        6.0 * rho2 * rho2 - 6.0 * rho2 + 1.0
    )
    base_field = simulator.base_amplitude * simulator.phase_factor

    def reconstruct(row: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
        phase_waves = (
            np.float32(row["z40_rms_waves"]) * z40_map
            + np.float32(row["z20_rms_waves"]) * z20_map
        )
        aberration = np.exp(
            1j * np.float32(2.0 * np.pi) * phase_waves
        ).astype(np.complex64)
        image = intensity(
            forward_fft((base_field * aberration).astype(np.complex64), np),
            np,
        ).astype(np.float32)
        normalized, metrics, _ = evaluate_image(image, simulator, experiment)
        return normalized, metrics

    baseline_image = simulator.reconstruct({})
    baseline, baseline_metrics, _ = evaluate_image(
        baseline_image, simulator, experiment
    )
    center_image, center_metrics = reconstruct(center_match)
    morphology_image, morphology_metrics = reconstruct(morphology_match)
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
            "Ideal V2 baseline",
            baseline,
            simulator.base.x_um,
            simulator.base.y_um,
            baseline_metrics,
            None,
        ),
        (
            "Rectangle-first center match",
            center_image,
            simulator.base.x_um,
            simulator.base.y_um,
            center_metrics,
            center_match,
        ),
        (
            "Closest internal morphology",
            morphology_image,
            simulator.base.x_um,
            simulator.base.y_um,
            morphology_metrics,
            morphology_match,
        ),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(20, 9), constrained_layout=True)
    colors = ["black", "#7570b3", "#1b9e77", "#d95f02"]
    for column, (title, image, x_um, y_um, metrics, row) in enumerate(panels):
        ax = axes[0, column]
        xmask = np.abs(x_um) <= 230.0
        ymask = np.abs(y_um) <= 105.0
        roi = image[np.ix_(ymask, xmask)]
        xa = x_um[xmask]
        ya = y_um[ymask]
        im = ax.imshow(
            roi,
            origin="upper",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            cmap="turbo",
            vmin=0.30,
            vmax=1.55,
            aspect="equal",
        )
        parameter = (
            ""
            if row is None
            else (
                f"\nZ40={row['z40_rms_waves']:+.4f}, "
                f"Z20={row['z20_rms_waves']:+.4f}, "
                f"IoU={row['footprint_iou']:.3f}"
            )
        )
        ax.set_title(
            f"{title}{parameter}\n"
            f"center={metrics['center_window_over_core_mean']:.3f}, "
            f"middle/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"RMS={100.0 * metrics['core_rms_fraction']:.1f}%"
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

        profile_ax = axes[1, column]
        px, profile_x = band_profile(
            image, x_um, y_um, "x", experiment["center_window_um"][1]
        )
        py, profile_y = band_profile(
            image, x_um, y_um, "y", experiment["center_window_um"][0]
        )
        profile_ax.plot(px, profile_x, color=colors[column], label="x profile")
        profile_ax.plot(
            py, profile_y, color=colors[column], ls="--", label="y profile"
        )
        profile_ax.set_xlim(-210.0, 210.0)
        profile_ax.set_ylim(0.0, 1.65)
        profile_ax.set_xlabel("position / um")
        profile_ax.set_ylabel("normalized intensity")
        profile_ax.grid(True, alpha=0.25)
        profile_ax.legend()
    fig.suptitle(
        "Zernike spherical + defocus reproduces an internal hollow "
        "without changing the Gaussian amplitude"
    )
    fig.savefig(RESULT_DIR / "KEY_ZERNIKE_COMPARISON.png", dpi=180)
    plt.close(fig)

    payload = {
        "experiment_metrics": {
            key: value
            for key, value in target.items()
            if not key.startswith("unit_profile")
        },
        "rectangle_first_center_match": center_match,
        "closest_internal_morphology_under_rectangle_constraint": morphology_match,
        "notes": [
            "Ideal 6.5 mm Gaussian amplitude and saved V2 phase are unchanged.",
            "Z40 and Z20 coefficients are RMS waves on the 15 mm clear pupil.",
            "The center-match case is selected only after choosing the best-envelope Z20 for each Z40.",
        ],
    }
    (RESULT_DIR / "KEY_ZERNIKE_CASES.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
