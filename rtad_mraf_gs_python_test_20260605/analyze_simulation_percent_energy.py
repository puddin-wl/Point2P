"""Analyze a refined simulation with the onsite 86.5% PercentEnergy metric.

Only ``reconstruction_refined.npy`` is analyzed. Installation compensation
such as X+5/Y+5 belongs to the real optical alignment and is deliberately not
forward-propagated with an ideal centered input beam.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PERCENT_ENERGY = 0.865
CAMERA_PIXEL_UM = 7.38
MANUAL_APERTURE_WIDTH_UM = 75.54526078852439 * CAMERA_PIXEL_UM
MANUAL_APERTURE_HEIGHT_UM = 28.18567502597709 * CAMERA_PIXEL_UM


def percent_energy_bbox(
    intensity: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
) -> dict[str, Any]:
    aperture = (
        (np.abs(x_um)[None, :] <= MANUAL_APERTURE_WIDTH_UM / 2.0)
        & (np.abs(y_um)[:, None] <= MANUAL_APERTURE_HEIGHT_UM / 2.0)
    )
    values = np.clip(intensity[aperture], 0.0, None)
    sorted_values = np.sort(values)[::-1]
    cumulative = np.cumsum(sorted_values, dtype=np.float64)
    threshold_index = int(np.searchsorted(cumulative, PERCENT_ENERGY * cumulative[-1]))
    threshold = float(sorted_values[min(threshold_index, sorted_values.size - 1)])
    selected = aperture & (intensity >= threshold)
    yy, xx = np.where(selected)
    xmin, xmax = int(xx.min()), int(xx.max())
    ymin, ymax = int(yy.min()), int(yy.max())
    dx_um = float(np.median(np.diff(x_um)))
    dy_um = float(np.median(np.diff(y_um)))
    return {
        "percent_energy": PERCENT_ENERGY,
        "aperture_center_um": [0.0, 0.0],
        "aperture_size_um": [MANUAL_APERTURE_WIDTH_UM, MANUAL_APERTURE_HEIGHT_UM],
        "sampling_um": [dx_um, dy_um],
        "threshold": threshold,
        "selected_pixel_count": int(selected.sum()),
        "selected_energy_fraction": float(
            intensity[selected].sum(dtype=np.float64)
            / intensity[aperture].sum(dtype=np.float64)
        ),
        "bbox_index_inclusive": [xmin, ymin, xmax, ymax],
        "bbox_count_px": [xmax - xmin + 1, ymax - ymin + 1],
        "bbox_center_span_px": [xmax - xmin, ymax - ymin],
        "bbox_center_span_um": [
            float(x_um[xmax] - x_um[xmin]),
            float(y_um[ymax] - y_um[ymin]),
        ],
        "bbox_pixel_coverage_um": [
            float((xmax - xmin + 1) * dx_um),
            float((ymax - ymin + 1) * dy_um),
        ],
    }


def analyze_run(run_dir: Path) -> dict[str, Any]:
    run_dir = run_dir.resolve()
    config = json.loads((run_dir / "config_used.json").read_text(encoding="utf-8"))
    target = np.load(run_dir / "target.npz")
    intensity_path = run_dir / "reconstruction_refined.npy"
    intensity = np.load(intensity_path).astype(np.float64)
    x_um = np.asarray(target["x_um"], dtype=np.float64)
    y_um = np.asarray(target["y_um"], dtype=np.float64)
    metric = percent_energy_bbox(intensity, x_um, y_um)

    outdir = run_dir / "percent_energy_analysis"
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "simulation_percent_energy_86p5_summary.json"
    csv_path = outdir / "simulation_percent_energy_86p5_summary.csv"
    figure_path = outdir / "simulation_percent_energy_86p5_diagnostics.png"

    result = {
        "run_dir": str(run_dir),
        "source": str(intensity_path),
        "run_parameters": {
            "target": config["target"],
            "physical": config["physical"],
        },
        "onsite_analysis_convention": {
            "percent_energy": PERCENT_ENERGY,
            "manual_aperture_size_um": [
                MANUAL_APERTURE_WIDTH_UM,
                MANUAL_APERTURE_HEIGHT_UM,
            ],
            "definition": "Sort 2D intensity pixels inside the equivalent onsite manual aperture, retain pixels until cumulative energy reaches 86.5%, and measure the retained-pixel bounding box.",
            "scope": "Analyze reconstruction_refined.npy only. X/Y installation compensation is for the real optical path and is not a second forward-propagation case.",
            "width_conventions": {
                "bbox_center_span_um": "Distance between the centers of the first and last retained samples; this is the convention used for the reported simulation size.",
                "bbox_pixel_coverage_um": "Full physical footprint of the retained discrete samples, retained only as a sampling diagnostic.",
            },
        },
        "simulation_86p5_percent_energy": metric,
        "outputs": {
            "summary_json": str(summary_path),
            "summary_csv": str(csv_path),
            "diagnostic_png": str(figure_path),
        },
    }
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8-sig") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "source",
                "span_x_um",
                "span_y_um",
                "coverage_x_um",
                "coverage_y_um",
                "sampling_x_um",
                "sampling_y_um",
            ]
        )
        writer.writerow(
            [
                "reconstruction_refined.npy",
                *metric["bbox_center_span_um"],
                *metric["bbox_pixel_coverage_um"],
                *metric["sampling_um"],
            ]
        )

    apw, aph = metric["aperture_size_um"]
    xmask = np.abs(x_um) <= apw / 2.0
    ymask = np.abs(y_um) <= aph / 2.0
    roi = intensity[np.ix_(ymask, xmask)]
    xa = x_um[xmask]
    ya = y_um[ymask]
    xmin, ymin, xmax, ymax = metric["bbox_index_inclusive"]
    span_x, span_y = metric["bbox_center_span_um"]
    cover_x, cover_y = metric["bbox_pixel_coverage_um"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6), constrained_layout=True)
    axes[0].imshow(
        roi,
        origin="upper",
        cmap="turbo",
        extent=[xa[0], xa[-1], ya[-1], ya[0]],
        vmin=0.0,
        vmax=float(np.percentile(roi, 99.5)),
        aspect="equal",
    )
    axes[0].add_patch(
        plt.Rectangle(
            (x_um[xmin], y_um[ymin]),
            x_um[xmax] - x_um[xmin],
            y_um[ymax] - y_um[ymin],
            fill=False,
            ec="magenta",
            lw=1.6,
        )
    )
    axes[0].set_title(f"86.5% PercentEnergy bbox\nX={span_x:.2f} um, Y={span_y:.2f} um")
    axes[0].set_xlabel("x / um")
    axes[0].set_ylabel("y / um")
    axes[0].text(
        0.02,
        0.04,
        f"sampling footprint={cover_x:.2f} x {cover_y:.2f} um",
        transform=axes[0].transAxes,
        color="white",
        fontsize=9,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
    )

    projection_x = roi.sum(axis=0)
    projection_y = roi.sum(axis=1)
    axes[1].plot(xa, projection_x / projection_x.max(), label="X projection")
    axes[1].plot(ya, projection_y / projection_y.max(), label="Y projection")
    axes[1].set_title("Intensity projections inside onsite aperture")
    axes[1].set_xlabel("position / um")
    axes[1].set_ylabel("normalized projected intensity")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend()
    fig.suptitle(f"Onsite-style simulation analysis\n{run_dir.name}")
    fig.savefig(figure_path, dpi=190)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    result = analyze_run(args.run_dir)
    metric = result["simulation_86p5_percent_energy"]
    print(
        "Simulation 86.5% size: "
        f"{metric['bbox_center_span_um'][0]:.2f} x "
        f"{metric['bbox_center_span_um'][1]:.2f} um"
    )
    print(result["outputs"]["summary_json"])


if __name__ == "__main__":
    main()
