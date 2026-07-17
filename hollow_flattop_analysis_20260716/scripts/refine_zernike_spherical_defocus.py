"""Fine scan of Noll-normalized Z40 spherical aberration plus Z20 defocus.

The scan follows the requested two-stage logic:

1. For every fixed Z40 coefficient, select the Z20 coefficient that best
   preserves the ideal V2 rectangular 50% footprint and outer band.
2. Only after that selection, inspect the internal hollow metrics.

The input stays the ideal 6.5 mm Gaussian and the saved V2 phase is unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
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
from src.propagation import forward_fft, intensity


DEFAULT_OUTPUT = (
    ANALYSIS_ROOT / "results" / "06_zernike_spherical_defocus_fine"
)


def footprint_size(
    binary: np.ndarray, x_um: np.ndarray, y_um: np.ndarray
) -> tuple[float, float]:
    yy, xx = np.nonzero(binary)
    if xx.size == 0:
        return float("nan"), float("nan")
    return (
        float(x_um[xx.max()] - x_um[xx.min()]),
        float(y_um[yy.max()] - y_um[yy.min()]),
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "z40_rms_waves",
        "z20_rms_waves",
        "envelope_score",
        "footprint_iou",
        "outer_band_rmse",
        "inside_fill_fraction",
        "outside_fraction",
        "size50_x_um",
        "size50_y_um",
        "rectangle_accepted",
        "center_ratio",
        "middle_over_sides",
        "core_rms_fraction",
        "hole_target_score",
    ]
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fields})


def plot_trend(
    path: Path, selected: list[dict[str, Any]], target: dict[str, float]
) -> None:
    z40 = np.asarray([row["z40_rms_waves"] for row in selected])
    z20 = np.asarray([row["z20_rms_waves"] for row in selected])
    iou = np.asarray([row["footprint_iou"] for row in selected])
    center = np.asarray([row["center_ratio"] for row in selected])
    middle = np.asarray([row["middle_over_sides"] for row in selected])
    rms = np.asarray([row["core_rms_fraction"] for row in selected])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    axes[0, 0].plot(z40, z20, "o-")
    axes[0, 0].set_title("Defocus selected to preserve rectangle")
    axes[0, 0].set_ylabel("selected Z20 / RMS waves")
    axes[0, 1].plot(z40, iou, "o-")
    axes[0, 1].axhline(0.82, color="gray", ls=":", label="acceptance")
    axes[0, 1].set_title("Rectangular 50% footprint")
    axes[0, 1].set_ylabel("IoU vs ideal V2")
    axes[0, 1].legend()
    axes[1, 0].plot(z40, center, "o-", label="simulation")
    axes[1, 0].axhline(
        target["center_window_over_core_mean"],
        color="black",
        ls=":",
        label="experiment",
    )
    axes[1, 0].set_title("Center depression after envelope tuning")
    axes[1, 0].set_ylabel("center window / core mean")
    axes[1, 0].legend()
    axes[1, 1].plot(z40, middle, "o-", label="middle / sides")
    axes[1, 1].plot(z40, rms, "s-", label="core RMS fraction")
    axes[1, 1].axhline(
        target["middle_third_over_side_thirds"],
        color="black",
        ls=":",
        label="experiment middle/sides",
    )
    axes[1, 1].set_title("Internal morphology")
    axes[1, 1].legend()
    for ax in axes.ravel():
        ax.set_xlabel("Z40 primary spherical / RMS waves")
        ax.grid(True, alpha=0.25)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_montage(
    path: Path,
    experiment: dict[str, Any],
    simulator: Simulator,
    baseline: np.ndarray,
    selected: list[dict[str, Any]],
    reconstruct: Any,
) -> None:
    # Five representative tuned cases spanning the Z40 range.
    indices = np.unique(
        np.rint(np.linspace(0, len(selected) - 1, 5)).astype(int)
    )
    panels: list[tuple[str, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]] = [
        (
            "Experiment",
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            experiment["metrics"],
        ),
        (
            "Ideal V2 baseline",
            baseline,
            simulator.base.x_um,
            simulator.base.y_um,
            evaluate_image(
                simulator.reconstruct({}), simulator, experiment
            )[1],
        ),
    ]
    for index in indices:
        row = selected[int(index)]
        image = reconstruct(row["z40_rms_waves"], row["z20_rms_waves"])
        normalized, metrics, _ = evaluate_image(image, simulator, experiment)
        panels.append(
            (
                f"Z40={row['z40_rms_waves']:+.3f}, "
                f"Z20={row['z20_rms_waves']:+.3f}\n"
                f"IoU={row['footprint_iou']:.3f}",
                normalized,
                simulator.base.x_um,
                simulator.base.y_um,
                metrics,
            )
        )

    fig, axes = plt.subplots(2, 4, figsize=(19, 9), constrained_layout=True)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (title, image, x_um, y_um, metrics) in zip(axes.ravel(), panels):
        ax.axis("on")
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
        ax.set_title(
            f"{title}\ncenter={metrics['center_window_over_core_mean']:.3f}, "
            f"middle/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"RMS={100.0 * metrics['core_rms_fraction']:.1f}%"
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(
        "Ideal Gaussian + fixed V2: Z40 spherical with envelope-tuned Z20 defocus"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(
    path: Path,
    selected: list[dict[str, Any]],
    accepted: list[dict[str, Any]],
    target: dict[str, float],
    count: int,
    elapsed: float,
) -> None:
    nonzero = [
        row
        for row in selected
        if abs(row["z40_rms_waves"]) > 1e-12 and row["rectangle_accepted"]
    ]
    deepest = min(nonzero, key=lambda row: row["center_ratio"]) if nonzero else None
    lines = [
        "# 球差 + 离焦细扫结果",
        "",
        "## 模型与判定顺序",
        "",
        "- 入射振幅保持理想 6.5 mm Gaussian，V2 相位保持不变。",
        "- 15 mm clear pupil 上采用 Noll 归一化 `Z40` 球差与 `Z20` 离焦，系数单位为 RMS waves。",
        "- 对每个固定 `Z40`，先选使矩形外轮廓评分最好的 `Z20`；再读取内部空洞指标。",
        "- 矩形接受条件：50% footprint IoU ≥ 0.82、内部填充率 ≥ 0.82、外溢比例 ≤ 0.20。",
        "",
        f"- 总案例数：{count}",
        f"- 计算耗时：{elapsed:.1f} s",
        f"- 调焦后仍通过矩形条件的球差点：{len(accepted)}/{len(selected)}",
        "",
        "## 实验空洞指标",
        "",
        f"- 中心/核心均值：`{target['center_window_over_core_mean']:.4f}`",
        f"- 中部/两侧：`{target['middle_third_over_side_thirds']:.4f}`",
        f"- 核心 RMS：`{100.0 * target['core_rms_fraction']:.2f}%`",
        "",
        "## 每个球差量经离焦调矩形后的结果",
        "",
        "| Z40 | 选中的 Z20 | IoU | 50%宽×高 (µm) | 中心比 | 中部/两侧 | RMS | 矩形通过 |",
        "|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in selected:
        lines.append(
            f"| {row['z40_rms_waves']:+.3f} | "
            f"{row['z20_rms_waves']:+.4f} | "
            f"{row['footprint_iou']:.3f} | "
            f"{row['size50_x_um']:.1f}×{row['size50_y_um']:.1f} | "
            f"{row['center_ratio']:.3f} | "
            f"{row['middle_over_sides']:.3f} | "
            f"{100.0 * row['core_rms_fraction']:.1f}% | "
            f"{'是' if row['rectangle_accepted'] else '否'} |"
        )
    lines.extend(["", "## 当前结论", ""])
    if deepest is None:
        lines.append(
            "- 在本次细扫范围内，非零球差即使配合离焦，也没有案例通过矩形外轮廓条件。"
        )
    else:
        lines.extend(
            [
                "- 在保持矩形的非零球差点中，中心最深的案例为：",
                f"  `Z40={deepest['z40_rms_waves']:+.3f}`、"
                f"`Z20={deepest['z20_rms_waves']:+.4f}` RMS waves；"
                f"IoU=`{deepest['footprint_iou']:.3f}`，"
                f"中心比=`{deepest['center_ratio']:.3f}`。",
                f"- 与实验中心比 `{target['center_window_over_core_mean']:.3f}` 比较，"
                + (
                    "该组合已达到或超过实验空洞深度。"
                    if deepest["center_ratio"]
                    <= target["center_window_over_core_mean"]
                    else "该组合尚未达到实验空洞深度。"
                ),
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
    baseline, _, _ = evaluate_image(baseline_image, simulator, experiment)

    x_um = simulator.base.x_um
    y_um = simulator.base.y_um
    X, Y = np.meshgrid(x_um, y_um)
    window = (np.abs(X) <= 240.0) & (np.abs(Y) <= 110.0)
    baseline_binary = (baseline >= 0.5) & window
    baseline_size_x, baseline_size_y = footprint_size(
        baseline_binary, x_um, y_um
    )
    rect_coordinate = np.maximum(
        np.abs(X) / (baseline_size_x / 2.0),
        np.abs(Y) / (baseline_size_y / 2.0),
    )
    outer_band = (rect_coordinate >= 0.68) & (rect_coordinate <= 1.02) & window
    target_rectangle = (
        (np.abs(X) <= baseline_size_x / 2.0)
        & (np.abs(Y) <= baseline_size_y / 2.0)
        & window
    )
    target_rectangle_count = max(1, int(np.count_nonzero(target_rectangle)))

    pupil_radius_m = float(simulator.clear_aperture_m / 2.0)
    rho2 = simulator.R2 / np.float32(pupil_radius_m * pupil_radius_m)
    z20_map = np.float32(math.sqrt(3.0)) * (2.0 * rho2 - 1.0)
    z40_map = np.float32(math.sqrt(5.0)) * (
        6.0 * rho2 * rho2 - 6.0 * rho2 + 1.0
    )
    base_field = simulator.base_amplitude * simulator.phase_factor

    def reconstruct(z40: float, z20: float) -> np.ndarray:
        phase_waves = np.float32(z40) * z40_map + np.float32(z20) * z20_map
        aberration = np.exp(
            1j * np.float32(2.0 * np.pi) * phase_waves
        ).astype(np.complex64)
        field = base_field * aberration
        return intensity(forward_fft(field.astype(np.complex64), np), np).astype(
            np.float32
        )

    spherical_values = np.linspace(
        args.spherical_min, args.spherical_max, args.spherical_count
    )
    defocus_values = np.linspace(
        args.defocus_min, args.defocus_max, args.defocus_count
    )
    total = int(spherical_values.size * defocus_values.size)
    target = experiment["metrics"]
    rows: list[dict[str, Any]] = []
    started = time.perf_counter()
    index = 0
    for z40 in spherical_values:
        for z20 in defocus_values:
            index += 1
            image = reconstruct(float(z40), float(z20))
            normalized, hollow, _ = evaluate_image(
                image, simulator, experiment
            )
            candidate_binary = (normalized >= 0.5) & window
            intersection = int(
                np.count_nonzero(baseline_binary & candidate_binary)
            )
            union = int(np.count_nonzero(baseline_binary | candidate_binary))
            iou = float(intersection / union) if union else 0.0
            outer_rmse = float(
                np.sqrt(
                    np.mean(
                        (normalized[outer_band] - baseline[outer_band]) ** 2
                    )
                )
            )
            fill = float(
                np.count_nonzero(candidate_binary & target_rectangle)
                / target_rectangle_count
            )
            candidate_count = max(1, int(np.count_nonzero(candidate_binary)))
            outside = float(
                np.count_nonzero(candidate_binary & ~target_rectangle)
                / candidate_count
            )
            size_x, size_y = footprint_size(candidate_binary, x_um, y_um)
            envelope_score = float(
                (1.0 - iou) / 0.04
                + outer_rmse / 0.12
                + (1.0 - fill) / 0.08
                + outside / 0.08
            )
            hole_score = float(
                math.sqrt(
                    (
                        (
                            hollow["center_window_over_core_mean"]
                            - target["center_window_over_core_mean"]
                        )
                        / 0.04
                    )
                    ** 2
                    + (
                        (
                            hollow["middle_third_over_side_thirds"]
                            - target["middle_third_over_side_thirds"]
                        )
                        / 0.04
                    )
                    ** 2
                    + (
                        (
                            hollow["core_rms_fraction"]
                            - target["core_rms_fraction"]
                        )
                        / 0.05
                    )
                    ** 2
                )
            )
            accepted = bool(iou >= 0.82 and fill >= 0.82 and outside <= 0.20)
            rows.append(
                {
                    "z40_rms_waves": float(z40),
                    "z20_rms_waves": float(z20),
                    "envelope_score": envelope_score,
                    "footprint_iou": iou,
                    "outer_band_rmse": outer_rmse,
                    "inside_fill_fraction": fill,
                    "outside_fraction": outside,
                    "size50_x_um": size_x,
                    "size50_y_um": size_y,
                    "rectangle_accepted": accepted,
                    "center_ratio": float(
                        hollow["center_window_over_core_mean"]
                    ),
                    "middle_over_sides": float(
                        hollow["middle_third_over_side_thirds"]
                    ),
                    "core_rms_fraction": float(hollow["core_rms_fraction"]),
                    "hole_target_score": hole_score,
                }
            )
            if index == 1 or index % 80 == 0 or index == total:
                print(
                    f"[{index}/{total}] Z40={z40:+.3f}, Z20={z20:+.4f}, "
                    f"IoU={iou:.3f}, center="
                    f"{hollow['center_window_over_core_mean']:.3f}",
                    flush=True,
                )

    elapsed = time.perf_counter() - started
    selected: list[dict[str, Any]] = []
    for z40 in spherical_values:
        group = [
            row
            for row in rows
            if math.isclose(row["z40_rms_waves"], float(z40), abs_tol=1e-12)
        ]
        selected.append(min(group, key=lambda row: row["envelope_score"]))
    accepted_selected = [
        row for row in selected if row["rectangle_accepted"]
    ]

    write_csv(args.outdir / "all_cases.csv", rows)
    write_csv(args.outdir / "best_defocus_for_each_spherical.csv", selected)
    payload = {
        "model": {
            "input": "ideal 6.5 mm Gaussian",
            "phase": str(args.case_dir / "phase_refined.npy"),
            "pupil_radius_mm": pupil_radius_m * 1e3,
            "Z20": "sqrt(3)*(2*rho^2-1)",
            "Z40": "sqrt(5)*(6*rho^4-6*rho^2+1)",
            "coefficient_unit": "RMS waves on 15 mm clear pupil",
        },
        "scan": {
            "z40": spherical_values.tolist(),
            "z20": defocus_values.tolist(),
            "case_count": total,
            "elapsed_seconds": elapsed,
        },
        "baseline_size50_um": [baseline_size_x, baseline_size_y],
        "experiment_metrics": target,
        "best_defocus_for_each_spherical": selected,
    }
    (args.outdir / "fine_scan_results.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_trend(
        args.outdir / "fine_scan_trend.png",
        selected,
        target,
    )
    plot_montage(
        args.outdir / "fine_scan_montage.png",
        experiment,
        simulator,
        baseline,
        selected,
        reconstruct,
    )
    write_summary(
        args.outdir / "SUMMARY_FINE_ZERNIKE.md",
        selected,
        accepted_selected,
        target,
        total,
        elapsed,
    )
    print(
        json.dumps(
            {
                "case_count": total,
                "elapsed_seconds": elapsed,
                "accepted_selected": len(accepted_selected),
                "selected": selected,
            },
            ensure_ascii=False,
            indent=2,
        ),
        flush=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--bmdata", type=Path, default=DEFAULT_BMDATA)
    parser.add_argument(
        "--percent-summary", type=Path, default=DEFAULT_PERCENT_SUMMARY
    )
    parser.add_argument(
        "--hollow-summary", type=Path, default=DEFAULT_HOLLOW_SUMMARY
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--spherical-min", type=float, default=-1.0)
    parser.add_argument("--spherical-max", type=float, default=1.0)
    parser.add_argument("--spherical-count", type=int, default=17)
    parser.add_argument("--defocus-min", type=float, default=-0.5)
    parser.add_argument("--defocus-max", type=float, default=0.5)
    parser.add_argument("--defocus-count", type=int, default=33)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
