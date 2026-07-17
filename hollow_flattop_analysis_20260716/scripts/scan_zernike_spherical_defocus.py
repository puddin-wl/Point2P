"""Scan Zernike primary spherical aberration together with extra defocus.

The ideal 6.5 mm Gaussian input and the fixed V2 phase are unchanged.

Noll-normalized forms on the physical 15 mm clear pupil are used:

    Z20 = sqrt(3) * (2 rho^2 - 1)
    Z40 = sqrt(5) * (6 rho^4 - 6 rho^2 + 1)

For every spherical coefficient, extra defocus is scanned.  Rectangle
preservation is evaluated first from the 50% footprint and outer-band
similarity to the ideal V2 result.  Center-depression metrics are evaluated
separately so a good rectangle is not automatically called a good hollow.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from dataclasses import replace
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
    evaluate_image,
    load_experimental_target,
)
from src.diagnostics import compute_diagnostics
from src.propagation import forward_fft, intensity


DEFAULT_OUTPUT = ANALYSIS_ROOT / "results" / "05_zernike_spherical_defocus"


def reconstruct_zernike(
    simulator: Simulator,
    z40_rms_waves: float,
    z20_rms_waves: float,
) -> np.ndarray:
    """Forward propagate V2 with Noll-normalized Z40 and Z20 wavefront terms."""
    pupil_radius_m = float(simulator.clear_aperture_m / 2.0)
    rho2 = simulator.R2 / np.float32(pupil_radius_m * pupil_radius_m)
    z20 = np.float32(math.sqrt(3.0)) * (2.0 * rho2 - 1.0)
    z40 = np.float32(math.sqrt(5.0)) * (
        6.0 * rho2 * rho2 - 6.0 * rho2 + 1.0
    )
    phase_waves = (
        np.float32(z40_rms_waves) * z40 + np.float32(z20_rms_waves) * z20
    )
    aberration = np.exp(1j * np.float32(2.0 * np.pi) * phase_waves).astype(
        np.complex64
    )
    field = simulator.base_amplitude * simulator.phase_factor * aberration
    return intensity(forward_fft(field.astype(np.complex64), np), np).astype(
        np.float32
    )


def envelope_metrics(
    normalized: np.ndarray,
    baseline_normalized: np.ndarray,
    simulator: Simulator,
    baseline_size50_x_um: float,
    baseline_size50_y_um: float,
) -> dict[str, float]:
    """Measure rectangular-envelope preservation independently of the center."""
    x_um = simulator.base.x_um
    y_um = simulator.base.y_um
    X, Y = np.meshgrid(x_um, y_um)
    analysis_window = (np.abs(X) <= 240.0) & (np.abs(Y) <= 110.0)
    baseline_binary = (baseline_normalized >= 0.5) & analysis_window
    candidate_binary = (normalized >= 0.5) & analysis_window
    intersection = int(np.count_nonzero(baseline_binary & candidate_binary))
    union = int(np.count_nonzero(baseline_binary | candidate_binary))
    iou = float(intersection / union) if union else 0.0

    rectangle_coordinate = np.maximum(
        np.abs(X) / (baseline_size50_x_um / 2.0),
        np.abs(Y) / (baseline_size50_y_um / 2.0),
    )
    outer_band = (
        (rectangle_coordinate >= 0.68)
        & (rectangle_coordinate <= 1.02)
        & analysis_window
    )
    outer_band_rmse = float(
        np.sqrt(
            np.mean(
                (normalized[outer_band] - baseline_normalized[outer_band]) ** 2
            )
        )
    )

    target_rectangle = (
        (np.abs(X) <= baseline_size50_x_um / 2.0)
        & (np.abs(Y) <= baseline_size50_y_um / 2.0)
        & analysis_window
    )
    inside_fill_fraction = float(
        np.count_nonzero(candidate_binary & target_rectangle)
        / max(1, np.count_nonzero(target_rectangle))
    )
    outside_fraction = float(
        np.count_nonzero(candidate_binary & analysis_window & ~target_rectangle)
        / max(1, np.count_nonzero(candidate_binary & analysis_window))
    )
    envelope_score = float(
        (1.0 - iou) / 0.04
        + outer_band_rmse / 0.12
        + (1.0 - inside_fill_fraction) / 0.08
        + outside_fraction / 0.08
    )
    return {
        "footprint_iou_vs_baseline": iou,
        "outer_band_rmse_vs_baseline": outer_band_rmse,
        "inside_fill_fraction": inside_fill_fraction,
        "outside_fraction": outside_fraction,
        "envelope_score": envelope_score,
    }


def build_row(
    simulator: Simulator,
    experiment: dict[str, Any],
    baseline_normalized: np.ndarray,
    baseline_size50_x_um: float,
    baseline_size50_y_um: float,
    z40: float,
    z20: float,
) -> tuple[dict[str, Any], np.ndarray]:
    image = reconstruct_zernike(simulator, z40, z20)
    normalized, hollow_metrics, experimental_score = evaluate_image(
        image, simulator, experiment
    )
    envelope = envelope_metrics(
        normalized,
        baseline_normalized,
        simulator,
        baseline_size50_x_um,
        baseline_size50_y_um,
    )
    data = replace(
        simulator.base,
        intensity_raw=image,
        intensity_source=f"V2 + Z40={z40:+.4f} + Z20={z20:+.4f} RMS waves",
    )
    diagnostics, _, warnings = compute_diagnostics(data)
    target = experiment["metrics"]
    hole_target_score = float(
        math.sqrt(
            (
                (
                    hollow_metrics["center_window_over_core_mean"]
                    - target["center_window_over_core_mean"]
                )
                / 0.04
            )
            ** 2
            + (
                (
                    hollow_metrics["middle_third_over_side_thirds"]
                    - target["middle_third_over_side_thirds"]
                )
                / 0.04
            )
            ** 2
            + (
                (
                    hollow_metrics["core_rms_fraction"]
                    - target["core_rms_fraction"]
                )
                / 0.05
            )
            ** 2
        )
    )
    row = {
        "family": "zernike_Z40_plus_Z20",
        "name": f"Z40{z40:+.3f}_Z20{z20:+.3f}_rmswaves",
        "params": {
            "z40_primary_spherical_rms_waves": z40,
            "z20_extra_defocus_rms_waves": z20,
            "pupil_radius_mm": simulator.clear_aperture_m * 0.5e3,
        },
        "envelope": envelope,
        "hollow_metrics": hollow_metrics,
        "hole_target_score": hole_target_score,
        "full_experimental_score": experimental_score["score"],
        "diagnostics": {
            "size50_x_um": diagnostics["size50_x_um"],
            "size50_y_um": diagnostics["size50_y_um"],
            "size90_x_um": diagnostics["size90_x_um"],
            "size90_y_um": diagnostics["size90_y_um"],
            "size13p5_x_um": diagnostics["size13p5_x_um"],
            "size13p5_y_um": diagnostics["size13p5_y_um"],
            "efficiency_e2_percent": diagnostics["efficiency_e2_percent"],
        },
        "warnings": warnings,
    }
    return row, image


def accepted_rectangle(row: dict[str, Any]) -> bool:
    envelope = row["envelope"]
    diagnostics = row["diagnostics"]
    return bool(
        envelope["footprint_iou_vs_baseline"] >= 0.82
        and envelope["inside_fill_fraction"] >= 0.82
        and envelope["outside_fraction"] <= 0.20
        and np.isfinite(diagnostics["size50_x_um"])
        and np.isfinite(diagnostics["size50_y_um"])
    )


def plot_selected_cases(
    output_path: Path,
    experiment: dict[str, Any],
    simulator: Simulator,
    baseline_normalized: np.ndarray,
    selected: list[dict[str, Any]],
) -> None:
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]] = [
        (
            "Experiment",
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            experiment["metrics"],
        ),
        (
            "Ideal V2 baseline",
            baseline_normalized,
            simulator.base.x_um,
            simulator.base.y_um,
            selected[0]["baseline_hollow_metrics"],
        ),
    ]
    for row in selected:
        image = reconstruct_zernike(
            simulator,
            row["params"]["z40_primary_spherical_rms_waves"],
            row["params"]["z20_extra_defocus_rms_waves"],
        )
        normalized, metrics, _ = evaluate_image(image, simulator, experiment)
        panels.append(
            (
                f"{row['name']}\nIoU={row['envelope']['footprint_iou_vs_baseline']:.3f}",
                normalized,
                simulator.base.x_um,
                simulator.base.y_um,
                metrics,
            )
        )
        if len(panels) >= 6:
            break

    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    for ax, (title, image, x_um, y_um, metrics) in zip(axes.ravel(), panels):
        xmask = np.abs(x_um) <= 230
        ymask = np.abs(y_um) <= 105
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
        ax.set_title(
            f"{title}\ncenter={metrics['center_window_over_core_mean']:.3f}, "
            f"middle/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"RMS={100.0 * metrics['core_rms_fraction']:.1f}%"
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(
        "Ideal Gaussian + fixed V2: Zernike primary spherical with tuned defocus"
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_spherical_trend(
    output_path: Path,
    best_envelope_by_spherical: list[dict[str, Any]],
    target: dict[str, Any],
) -> None:
    rows = sorted(
        best_envelope_by_spherical,
        key=lambda row: row["params"]["z40_primary_spherical_rms_waves"],
    )
    spherical = np.asarray(
        [row["params"]["z40_primary_spherical_rms_waves"] for row in rows]
    )
    defocus = np.asarray(
        [row["params"]["z20_extra_defocus_rms_waves"] for row in rows]
    )
    center = np.asarray(
        [row["hollow_metrics"]["center_window_over_core_mean"] for row in rows]
    )
    middle = np.asarray(
        [row["hollow_metrics"]["middle_third_over_side_thirds"] for row in rows]
    )
    iou = np.asarray(
        [row["envelope"]["footprint_iou_vs_baseline"] for row in rows]
    )
    rms = np.asarray(
        [100.0 * row["hollow_metrics"]["core_rms_fraction"] for row in rows]
    )

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    axes[0, 0].plot(spherical, defocus, "o-")
    axes[0, 0].set_ylabel("selected extra Z20 / RMS waves")
    axes[0, 0].set_title("Defocus chosen for best rectangular envelope")
    axes[0, 1].plot(spherical, iou, "o-")
    axes[0, 1].axhline(0.82, color="gray", ls=":", label="acceptance")
    axes[0, 1].set_ylabel("50% footprint IoU vs ideal V2")
    axes[0, 1].set_title("Rectangle preservation")
    axes[0, 1].legend()
    axes[1, 0].plot(spherical, center, "o-", label="simulation")
    axes[1, 0].axhline(
        target["center_window_over_core_mean"],
        color="black",
        ls=":",
        label="experiment",
    )
    axes[1, 0].set_ylabel("center window / core mean")
    axes[1, 0].set_title("Internal center depression")
    axes[1, 0].legend()
    axes[1, 1].plot(spherical, middle, "o-", label="middle / sides")
    axes[1, 1].axhline(
        target["middle_third_over_side_thirds"],
        color="black",
        ls=":",
        label="experiment",
    )
    axes[1, 1].plot(spherical, rms / 100.0, "s-", label="core RMS fraction")
    axes[1, 1].set_title("Internal morphology")
    axes[1, 1].legend()
    for ax in axes.ravel():
        ax.set_xlabel("Z40 primary spherical / RMS waves")
        ax.grid(True, alpha=0.25)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "rank_envelope",
        "name",
        "z40_rms_waves",
        "z20_rms_waves",
        "envelope_score",
        "footprint_iou",
        "inside_fill_fraction",
        "outside_fraction",
        "outer_band_rmse",
        "rectangle_accepted",
        "center_ratio",
        "middle_over_sides",
        "core_rms_fraction",
        "hole_target_score",
        "size50_x_um",
        "size50_y_um",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow(
                {
                    "rank_envelope": rank,
                    "name": row["name"],
                    "z40_rms_waves": row["params"][
                        "z40_primary_spherical_rms_waves"
                    ],
                    "z20_rms_waves": row["params"][
                        "z20_extra_defocus_rms_waves"
                    ],
                    "envelope_score": row["envelope"]["envelope_score"],
                    "footprint_iou": row["envelope"][
                        "footprint_iou_vs_baseline"
                    ],
                    "inside_fill_fraction": row["envelope"][
                        "inside_fill_fraction"
                    ],
                    "outside_fraction": row["envelope"]["outside_fraction"],
                    "outer_band_rmse": row["envelope"][
                        "outer_band_rmse_vs_baseline"
                    ],
                    "rectangle_accepted": accepted_rectangle(row),
                    "center_ratio": row["hollow_metrics"][
                        "center_window_over_core_mean"
                    ],
                    "middle_over_sides": row["hollow_metrics"][
                        "middle_third_over_side_thirds"
                    ],
                    "core_rms_fraction": row["hollow_metrics"][
                        "core_rms_fraction"
                    ],
                    "hole_target_score": row["hole_target_score"],
                    "size50_x_um": row["diagnostics"]["size50_x_um"],
                    "size50_y_um": row["diagnostics"]["size50_y_um"],
                }
            )


def write_summary(
    path: Path,
    rows: list[dict[str, Any]],
    best_by_spherical: list[dict[str, Any]],
    best_rectangular_hollow: dict[str, Any] | None,
    experiment: dict[str, Any],
    elapsed_seconds: float,
) -> None:
    target = experiment["metrics"]
    lines = [
        "# Zernike 球差 + 离焦扫描",
        "",
        "## 模型",
        "",
        "- 输入振幅保持理想 6.5 mm Gaussian，不加入 Gaussian 缺陷。",
        "- V2 `phase_refined.npy` 保持不变。",
        "- 物理 pupil 直径使用配置中的 15 mm clear aperture。",
        "- 使用 Noll 归一化 `Z40` 主球差和独立 `Z20` 离焦。",
        "- 对每个 `Z40` 系数，先按矩形外轮廓选择最佳 `Z20`，再检查内部。",
        "",
        f"- 扫描案例：{len(rows)}",
        f"- 耗时：{elapsed_seconds:.1f} s",
        "",
        "## 实验内部指标",
        "",
        f"- 中心比：`{target['center_window_over_core_mean']:.4f}`",
        f"- 中部/两侧：`{target['middle_third_over_side_thirds']:.4f}`",
        f"- 核心 RMS：`{100.0 * target['core_rms_fraction']:.2f}%`",
        "",
        "## 每个球差量对应的最佳矩形离焦",
        "",
        "| Z40 RMS waves | 选择的 Z20 | footprint IoU | 中心比 | 中部/两侧 | RMS |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(
        best_by_spherical,
        key=lambda item: item["params"]["z40_primary_spherical_rms_waves"],
    ):
        lines.append(
            f"| {row['params']['z40_primary_spherical_rms_waves']:+.3f} | "
            f"{row['params']['z20_extra_defocus_rms_waves']:+.3f} | "
            f"{row['envelope']['footprint_iou_vs_baseline']:.3f} | "
            f"{row['hollow_metrics']['center_window_over_core_mean']:.3f} | "
            f"{row['hollow_metrics']['middle_third_over_side_thirds']:.3f} | "
            f"{100.0 * row['hollow_metrics']['core_rms_fraction']:.1f}% |"
        )
    lines.extend(["", "## 矩形保持条件下最接近实验空洞的案例", ""])
    if best_rectangular_hollow is None:
        lines.append("没有案例通过当前矩形外轮廓筛选。")
    else:
        row = best_rectangular_hollow
        lines.extend(
            [
                f"- 名称：`{row['name']}`",
                f"- Z40：`{row['params']['z40_primary_spherical_rms_waves']:+.4f}` RMS waves",
                f"- 额外 Z20：`{row['params']['z20_extra_defocus_rms_waves']:+.4f}` RMS waves",
                f"- footprint IoU：`{row['envelope']['footprint_iou_vs_baseline']:.4f}`",
                f"- 中心比：`{row['hollow_metrics']['center_window_over_core_mean']:.4f}`",
                f"- 中部/两侧：`{row['hollow_metrics']['middle_third_over_side_thirds']:.4f}`",
                f"- 核心 RMS：`{100.0 * row['hollow_metrics']['core_rms_fraction']:.2f}%`",
            ]
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        args.bmdata, args.percent_summary, args.hollow_summary
    )
    simulator = Simulator(args.case_dir)
    baseline_image = simulator.reconstruct({})
    baseline_normalized, baseline_hollow_metrics, _ = evaluate_image(
        baseline_image, simulator, experiment
    )
    baseline_data = replace(
        simulator.base,
        intensity_raw=baseline_image,
        intensity_source="recomputed ideal V2 baseline",
    )
    baseline_diagnostics, _, _ = compute_diagnostics(baseline_data)
    baseline_size50_x_um = float(baseline_diagnostics["size50_x_um"])
    baseline_size50_y_um = float(baseline_diagnostics["size50_y_um"])

    spherical_values = np.linspace(
        args.spherical_min, args.spherical_max, args.spherical_count
    )
    defocus_values = np.linspace(
        args.defocus_min, args.defocus_max, args.defocus_count
    )
    total = len(spherical_values) * len(defocus_values)
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    index = 0
    for spherical in spherical_values:
        for defocus in defocus_values:
            index += 1
            row, _ = build_row(
                simulator,
                experiment,
                baseline_normalized,
                baseline_size50_x_um,
                baseline_size50_y_um,
                float(spherical),
                float(defocus),
            )
            row["baseline_hollow_metrics"] = baseline_hollow_metrics
            rows.append(row)
            if index == 1 or index % 40 == 0 or index == total:
                print(
                    f"[{index:3d}/{total}] {row['name']} "
                    f"envelope={row['envelope']['envelope_score']:.3f}, "
                    f"center={row['hollow_metrics']['center_window_over_core_mean']:.3f}",
                    flush=True,
                )
    elapsed = time.perf_counter() - started
    rows_by_envelope = sorted(
        rows, key=lambda row: row["envelope"]["envelope_score"]
    )
    best_by_spherical: list[dict[str, Any]] = []
    for spherical in spherical_values:
        group = [
            row
            for row in rows
            if math.isclose(
                row["params"]["z40_primary_spherical_rms_waves"],
                float(spherical),
                abs_tol=1e-12,
            )
        ]
        best_by_spherical.append(
            min(group, key=lambda row: row["envelope"]["envelope_score"])
        )

    accepted = [row for row in rows if accepted_rectangle(row)]
    best_rectangular_hollow = (
        min(accepted, key=lambda row: row["hole_target_score"])
        if accepted
        else None
    )
    selected_for_plot = sorted(
        accepted, key=lambda row: row["hole_target_score"]
    )[:4]
    if not selected_for_plot:
        selected_for_plot = rows_by_envelope[:4]

    serializable = {
        "phase_source": str(args.case_dir / "phase_refined.npy"),
        "input_model": "ideal 6.5 mm Gaussian from V2 config",
        "zernike_definition": {
            "pupil_radius_mm": simulator.clear_aperture_m * 0.5e3,
            "Z20": "sqrt(3) * (2*rho^2 - 1)",
            "Z40": "sqrt(5) * (6*rho^4 - 6*rho^2 + 1)",
            "coefficient_unit": "RMS waves on the 15 mm clear pupil",
        },
        "baseline": {
            "size50_x_um": baseline_size50_x_um,
            "size50_y_um": baseline_size50_y_um,
            "hollow_metrics": baseline_hollow_metrics,
        },
        "experiment_metrics": experiment["metrics"],
        "candidate_count": len(rows),
        "elapsed_seconds": elapsed,
        "rectangle_acceptance": {
            "footprint_iou_min": 0.82,
            "inside_fill_min": 0.82,
            "outside_fraction_max": 0.20,
            "accepted_count": len(accepted),
        },
        "best_envelope_by_spherical": best_by_spherical,
        "best_rectangular_hollow": best_rectangular_hollow,
        "candidates": rows_by_envelope,
    }
    (args.outdir / "zernike_spherical_defocus_results.json").write_text(
        json.dumps(serializable, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_csv(args.outdir / "zernike_spherical_defocus_results.csv", rows_by_envelope)
    plot_selected_cases(
        args.outdir / "zernike_spherical_defocus_overview.png",
        experiment,
        simulator,
        baseline_normalized,
        selected_for_plot,
    )
    plot_spherical_trend(
        args.outdir / "zernike_spherical_defocus_trend.png",
        best_by_spherical,
        experiment["metrics"],
    )
    write_summary(
        args.outdir / "SUMMARY_ZERNIKE_SPHERICAL_DEFOCUS.md",
        rows,
        best_by_spherical,
        best_rectangular_hollow,
        experiment,
        elapsed,
    )
    print(
        json.dumps(
            {
                "candidate_count": len(rows),
                "accepted_rectangles": len(accepted),
                "elapsed_seconds": elapsed,
                "best_rectangular_hollow": (
                    {
                        "name": best_rectangular_hollow["name"],
                        "params": best_rectangular_hollow["params"],
                        "envelope": best_rectangular_hollow["envelope"],
                        "hollow_metrics": {
                            key: value
                            for key, value in best_rectangular_hollow[
                                "hollow_metrics"
                            ].items()
                            if not key.startswith("unit_profile")
                        },
                        "hole_target_score": best_rectangular_hollow[
                            "hole_target_score"
                        ],
                    }
                    if best_rectangular_hollow
                    else None
                ),
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
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--spherical-min", type=float, default=-12.0)
    parser.add_argument("--spherical-max", type=float, default=12.0)
    parser.add_argument("--spherical-count", type=int, default=9)
    parser.add_argument("--defocus-min", type=float, default=-6.0)
    parser.add_argument("--defocus-max", type=float, default=6.0)
    parser.add_argument("--defocus-count", type=int, default=33)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
