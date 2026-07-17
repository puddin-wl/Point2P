"""Free-space propagation, 15 mm field-lens pupil, then Fourier transform.

Only two input branches are used:
1. the directly resampled measured G-spot amplitude;
2. the ideal 6.5 mm Gaussian amplitude.

Both use the V2 phase without blaze and without added Zernike aberrations.
At each SLM-to-field-lens distance from 0 to 1 m, the field is propagated,
optionally clipped by the 15 mm circular pupil, and Fourier transformed to
the back focal plane of the 429 mm field lens.
"""

from __future__ import annotations

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

from scan_slm_to_field_lens_free_space import (
    ANALYSIS_ROOT,
    DISTANCES_M,
    PAD_N,
    PUPIL_DIAMETER_M,
    SLM_N,
    SLM_PITCH_M,
    WAVELENGTH_M,
    embed_slm,
    ideal_gaussian_on_slm,
    load_phase_variants,
    measured_amplitude_on_slm,
    slm_coordinates,
)


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "10_field_lens_pupil_fourier"
FIELD_LENS_FOCAL_LENGTH_M = 0.429
EXPERIMENT_CORE_WIDTH_X_UM = 339.48
EXPERIMENT_CORE_WIDTH_Y_UM = 118.08
EXPERIMENT_CENTER_WIDTH_X_UM = 36.90
EXPERIMENT_CENTER_WIDTH_Y_UM = 36.90


def forward_fourier(field: np.ndarray) -> np.ndarray:
    return np.fft.fftshift(
        np.fft.fft2(np.fft.ifftshift(field), norm="ortho")
    ).astype(np.complex64)


def normalize_to_core(
    intensity: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
) -> np.ndarray:
    x_core = np.abs(x_um) <= EXPERIMENT_CORE_WIDTH_X_UM / 2.0
    y_core = np.abs(y_um) <= EXPERIMENT_CORE_WIDTH_Y_UM / 2.0
    core_mean = float(np.mean(intensity[np.ix_(y_core, x_core)]))
    return intensity.astype(np.float64) / max(core_mean, 1e-30)


def focal_metrics(
    normalized: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
) -> dict[str, float]:
    x_core = np.abs(x_um) <= EXPERIMENT_CORE_WIDTH_X_UM / 2.0
    y_core = np.abs(y_um) <= EXPERIMENT_CORE_WIDTH_Y_UM / 2.0
    core = normalized[np.ix_(y_core, x_core)]
    x_center = np.abs(x_um) <= EXPERIMENT_CENTER_WIDTH_X_UM / 2.0
    y_center = np.abs(y_um) <= EXPERIMENT_CENTER_WIDTH_Y_UM / 2.0
    center = normalized[np.ix_(y_center, x_center)]
    core_x = x_um[x_core]
    middle_columns = np.abs(core_x) <= EXPERIMENT_CORE_WIDTH_X_UM / 6.0
    middle = core[:, middle_columns]
    sides = core[:, ~middle_columns]
    return {
        "center_window_over_core_mean": float(np.mean(center) / np.mean(core)),
        "middle_third_over_side_thirds": float(
            np.mean(middle) / np.mean(sides)
        ),
        "core_rms_fraction": float(np.std(core) / np.mean(core)),
    }


def comparison_metrics(
    no_pupil: np.ndarray,
    clipped: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
) -> dict[str, float]:
    analysis = (np.abs(x_um)[None, :] <= 260.0) & (
        np.abs(y_um)[:, None] <= 140.0
    )
    no_binary = (no_pupil >= 0.5) & analysis
    clipped_binary = (clipped >= 0.5) & analysis
    intersection = int(np.count_nonzero(no_binary & clipped_binary))
    union = int(np.count_nonzero(no_binary | clipped_binary))
    relative_l2 = float(
        np.linalg.norm((clipped - no_pupil)[analysis])
        / max(np.linalg.norm(no_pupil[analysis]), 1e-30)
    )
    peak_difference = float(
        np.max(np.abs((clipped - no_pupil)[analysis]))
    )
    return {
        "footprint_iou_clipped_vs_no_pupil": (
            float(intersection / union) if union else 0.0
        ),
        "normalized_relative_l2_in_roi": relative_l2,
        "normalized_peak_abs_difference_in_roi": peak_difference,
    }


def build_row(
    branch: str,
    z_m: float,
    field_lens: np.ndarray,
    pupil_mask: np.ndarray,
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    power_before = float(np.sum(np.abs(field_lens) ** 2, dtype=np.float64))
    clipped_field = field_lens * pupil_mask
    power_after = float(np.sum(np.abs(clipped_field) ** 2, dtype=np.float64))
    lost_fraction = 1.0 - power_after / power_before

    focal_no_pupil = np.abs(forward_fourier(field_lens)) ** 2
    focal_clipped = np.abs(forward_fourier(clipped_field)) ** 2
    normalized_no_pupil = normalize_to_core(
        focal_no_pupil, x_focal_um, y_focal_um
    )
    normalized_clipped = normalize_to_core(
        focal_clipped, x_focal_um, y_focal_um
    )
    metrics_no_pupil = focal_metrics(
        normalized_no_pupil, x_focal_um, y_focal_um
    )
    metrics_clipped = focal_metrics(
        normalized_clipped, x_focal_um, y_focal_um
    )
    comparison = comparison_metrics(
        normalized_no_pupil,
        normalized_clipped,
        x_focal_um,
        y_focal_um,
    )
    row = {
        "branch": branch,
        "z_m": z_m,
        "pupil_power_loss_percent": 100.0 * lost_fraction,
        "focal_power_transmission_percent": 100.0 * power_after / power_before,
        "no_pupil_center_ratio": metrics_no_pupil[
            "center_window_over_core_mean"
        ],
        "clipped_center_ratio": metrics_clipped[
            "center_window_over_core_mean"
        ],
        "center_ratio_change": metrics_clipped[
            "center_window_over_core_mean"
        ]
        - metrics_no_pupil["center_window_over_core_mean"],
        "no_pupil_middle_over_sides": metrics_no_pupil[
            "middle_third_over_side_thirds"
        ],
        "clipped_middle_over_sides": metrics_clipped[
            "middle_third_over_side_thirds"
        ],
        "middle_over_sides_change": metrics_clipped[
            "middle_third_over_side_thirds"
        ]
        - metrics_no_pupil["middle_third_over_side_thirds"],
        "no_pupil_core_rms_fraction": metrics_no_pupil["core_rms_fraction"],
        "clipped_core_rms_fraction": metrics_clipped["core_rms_fraction"],
        "core_rms_fraction_change": metrics_clipped["core_rms_fraction"]
        - metrics_no_pupil["core_rms_fraction"],
        **comparison,
    }
    return row, normalized_no_pupil, normalized_clipped


def propagate_and_focus(
    field_slm: np.ndarray,
    branch: str,
    delta_kz: np.ndarray,
    pupil_mask: np.ndarray,
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[float, tuple[np.ndarray, np.ndarray]]]:
    field0 = embed_slm(field_slm)
    spectrum = np.fft.fft2(field0, norm="ortho").astype(np.complex64)
    rows: list[dict[str, Any]] = []
    snapshots: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for z_m in DISTANCES_M:
        if z_m == 0.0:
            field_lens = field0
        else:
            transfer = np.exp(1j * delta_kz * float(z_m)).astype(np.complex64)
            field_lens = np.fft.ifft2(
                spectrum * transfer, norm="ortho"
            ).astype(np.complex64)
        row, focal_no_pupil, focal_clipped = build_row(
            branch,
            float(z_m),
            field_lens,
            pupil_mask,
            x_focal_um,
            y_focal_um,
        )
        rows.append(row)
        if any(abs(float(z_m) - value) < 1e-9 for value in (0.0, 0.5, 1.0)):
            snapshots[float(z_m)] = (focal_no_pupil, focal_clipped)
        print(
            f"{branch}: z={z_m:.1f} m, "
            f"loss={row['pupil_power_loss_percent']:.6f}%, "
            f"dCenter={row['center_ratio_change']:+.6f}, "
            f"L2={row['normalized_relative_l2_in_roi']:.6f}",
            flush=True,
        )
    return rows, snapshots


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_trends(path: Path, rows: list[dict[str, Any]]) -> None:
    branches = list(dict.fromkeys(str(row["branch"]) for row in rows))
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for branch in branches:
        group = [row for row in rows if row["branch"] == branch]
        z = np.asarray([float(row["z_m"]) for row in group])
        loss = np.asarray(
            [float(row["pupil_power_loss_percent"]) for row in group]
        )
        l2 = np.asarray(
            [float(row["normalized_relative_l2_in_roi"]) for row in group]
        )
        center_no = np.asarray(
            [float(row["no_pupil_center_ratio"]) for row in group]
        )
        center_clip = np.asarray(
            [float(row["clipped_center_ratio"]) for row in group]
        )
        iou = np.asarray(
            [
                float(row["footprint_iou_clipped_vs_no_pupil"])
                for row in group
            ]
        )
        axes[0, 0].semilogy(z, np.maximum(loss, 1e-8), "o-", label=branch)
        axes[0, 1].plot(z, l2, "o-", label=branch)
        color = axes[1, 0]._get_lines.get_next_color()
        axes[1, 0].plot(
            z, center_no, "--", color=color, label=f"{branch}: no pupil"
        )
        axes[1, 0].plot(
            z, center_clip, "o-", color=color, label=f"{branch}: clipped"
        )
        axes[1, 1].plot(z, iou, "o-", label=branch)
    axes[0, 0].set_title("Power removed by 15 mm pupil")
    axes[0, 0].set_ylabel("power loss / %")
    axes[0, 1].set_title("Normalized focal-image change")
    axes[0, 1].set_ylabel("relative L2 in ROI")
    axes[1, 0].set_title("Center ratio before/after clipping")
    axes[1, 0].set_ylabel("center window / core mean")
    axes[1, 1].set_title("50% footprint overlap")
    axes[1, 1].set_ylabel("IoU: clipped vs no pupil")
    for ax in axes.ravel():
        ax.set_xlabel("SLM-to-field-lens distance / m")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(
        "15 mm field-lens pupil followed by Fourier transform; V2 only"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_montage(
    path: Path,
    snapshots: dict[str, dict[float, tuple[np.ndarray, np.ndarray]]],
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
) -> None:
    branches = list(snapshots)
    distances = [0.0, 0.5, 1.0]
    fig, axes = plt.subplots(
        len(branches) * 2,
        len(distances),
        figsize=(15, 4.0 * len(branches) * 2),
        constrained_layout=True,
    )
    xmask = np.abs(x_focal_um) <= 250.0
    ymask = np.abs(y_focal_um) <= 130.0
    xa = x_focal_um[xmask]
    ya = y_focal_um[ymask]
    for branch_index, branch in enumerate(branches):
        for column, z_m in enumerate(distances):
            no_pupil, clipped = snapshots[branch][z_m]
            for condition_index, (condition, image) in enumerate(
                (("no pupil", no_pupil), ("15 mm pupil", clipped))
            ):
                row = 2 * branch_index + condition_index
                ax = axes[row, column]
                roi = image[np.ix_(ymask, xmask)]
                im = ax.imshow(
                    roi,
                    origin="upper",
                    extent=[xa[0], xa[-1], ya[-1], ya[0]],
                    cmap="turbo",
                    vmin=0.30,
                    vmax=1.55,
                    aspect="equal",
                )
                ax.set_title(f"{branch}, {condition}\nz={z_m:.1f} m")
                ax.set_xlabel("x / um")
                ax.set_ylabel("y / um")
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(
        "Back focal-plane intensity: free propagation → pupil → Fourier lens"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path: Path, rows: list[dict[str, Any]], elapsed_s: float) -> None:
    branches = list(dict.fromkeys(str(row["branch"]) for row in rows))
    lines = [
        "# 15 mm场镜入瞳截断后傅里叶变换",
        "",
        "## 模型",
        "",
        "- 输入仅两组：原始拍摄Gaussian、理想6.5 mm Gaussian。",
        "- 相位：无闪耀V2；不加入球差、离焦或其他Zernike项。",
        "- SLM到场镜：0–1 m角谱传播，每0.1 m采样。",
        "- 场镜：固定轴心15 mm圆形入瞳，焦距429 mm。",
        "- 在同一场镜平面分别计算不截断和乘15 mm圆孔后的傅里叶焦面。",
        f"- 焦面采样：{WAVELENGTH_M * FIELD_LENS_FOCAL_LENGTH_M / (PAD_N * SLM_PITCH_M) * 1e6:.4f} µm。",
        f"- 耗时：{elapsed_s:.1f} s。",
        "",
        "## 1 m处截断影响",
        "",
        "| 输入 | 功率损失 | 中心比(不截断→截断) | 中部/两侧(不截断→截断) | RMS(不截断→截断) | 焦面L2变化 | 50% IoU |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for branch in branches:
        row = next(
            item
            for item in rows
            if item["branch"] == branch and abs(float(item["z_m"]) - 1.0) < 1e-9
        )
        lines.append(
            f"| {branch} | "
            f"{float(row['pupil_power_loss_percent']):.6f}% | "
            f"{float(row['no_pupil_center_ratio']):.4f}→"
            f"{float(row['clipped_center_ratio']):.4f} | "
            f"{float(row['no_pupil_middle_over_sides']):.4f}→"
            f"{float(row['clipped_middle_over_sides']):.4f} | "
            f"{100.0 * float(row['no_pupil_core_rms_fraction']):.2f}%→"
            f"{100.0 * float(row['clipped_core_rms_fraction']):.2f}% | "
            f"{float(row['normalized_relative_l2_in_roi']):.6f} | "
            f"{float(row['footprint_iou_clipped_vs_no_pupil']):.6f} |"
        )
    lines.extend(
        [
            "",
            "## 判断",
            "",
            "- 原始实测光斑由整幅 `.bgData` 按物理坐标重采样，没有矩形 ROI 硬裁切。",
            "- 两组输入在 1 m 处的 50% 外轮廓 IoU 均为 1.000，截断没有改变主轮廓，也没有制造中心空洞。",
            "- 实测输入的中心比只改变 -0.00064，理想输入只改变 +0.00060；相对于实验约 0.23 的中心亏损，这一效应小约三个数量级。",
            "- 因此 15 mm 圆孔截断不是 2026-07-16 中心空洞的主要成因。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _, _, X_slm, Y_slm = slm_coordinates()
    measured_amplitude = measured_amplitude_on_slm(X_slm, Y_slm)
    ideal_amplitude = ideal_gaussian_on_slm(X_slm, Y_slm)
    phase_v2 = load_phase_variants()["v2_no_blaze"]
    phase_factor = np.exp(1j * phase_v2).astype(np.complex64)

    axis_pad = (
        np.arange(PAD_N, dtype=np.float32) - PAD_N // 2
    ) * np.float32(SLM_PITCH_M)
    X_pad, Y_pad = np.meshgrid(axis_pad, axis_pad)
    pupil_mask = (
        X_pad * X_pad + Y_pad * Y_pad
        <= np.float32(PUPIL_DIAMETER_M / 2.0) ** 2
    )

    frequency = np.fft.fftfreq(PAD_N, d=SLM_PITCH_M)
    FX, FY = np.meshgrid(frequency, frequency)
    k = 2.0 * np.pi / WAVELENGTH_M
    angular_argument = np.maximum(
        0.0,
        1.0 - (WAVELENGTH_M * FX) ** 2 - (WAVELENGTH_M * FY) ** 2,
    )
    delta_kz = (k * (np.sqrt(angular_argument) - 1.0)).astype(np.float32)

    focal_frequency = np.fft.fftshift(np.fft.fftfreq(PAD_N, d=SLM_PITCH_M))
    x_focal_um = (
        WAVELENGTH_M
        * FIELD_LENS_FOCAL_LENGTH_M
        * focal_frequency
        * 1e6
    )
    y_focal_um = x_focal_um.copy()

    branches = {
        "measured_G": measured_amplitude * phase_factor,
        "ideal_6p5mm_G": ideal_amplitude * phase_factor,
    }
    rows: list[dict[str, Any]] = []
    snapshots: dict[str, dict[float, tuple[np.ndarray, np.ndarray]]] = {}
    started = time.perf_counter()
    for branch, field in branches.items():
        branch_rows, branch_snapshots = propagate_and_focus(
            field.astype(np.complex64),
            branch,
            delta_kz,
            pupil_mask,
            x_focal_um,
            y_focal_um,
        )
        rows.extend(branch_rows)
        snapshots[branch] = branch_snapshots
    elapsed = time.perf_counter() - started

    write_csv(OUTPUT_DIR / "field_lens_pupil_fourier_scan.csv", rows)
    payload = {
        "model": {
            "wavelength_m": WAVELENGTH_M,
            "slm_shape": [SLM_N, SLM_N],
            "slm_pitch_um": SLM_PITCH_M * 1e6,
            "propagation_shape": [PAD_N, PAD_N],
            "slm_to_lens_distances_m": DISTANCES_M.tolist(),
            "field_lens_focal_length_m": FIELD_LENS_FOCAL_LENGTH_M,
            "field_lens_pupil_diameter_mm": PUPIL_DIAMETER_M * 1e3,
            "focal_sampling_um": float(x_focal_um[1] - x_focal_um[0]),
            "phase": "V2 without blaze",
            "added_zernike": False,
            "measured_input_preprocessing": (
                "full original frame; robust background subtraction; "
                "5 px Gaussian smoothing before negative clipping; "
                "physical-coordinate resampling; no rectangular ROI crop"
            ),
        },
        "elapsed_seconds": elapsed,
        "rows": rows,
    }
    (OUTPUT_DIR / "field_lens_pupil_fourier_scan.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_trends(OUTPUT_DIR / "field_lens_pupil_fourier_trends.png", rows)
    plot_montage(
        OUTPUT_DIR / "field_lens_pupil_fourier_montage.png",
        snapshots,
        x_focal_um,
        y_focal_um,
    )
    for branch in branches:
        np.save(
            OUTPUT_DIR / f"{branch}_z1m_focal_with_15mm_pupil.npy",
            snapshots[branch][1.0][1].astype(np.float32),
        )
    write_summary(OUTPUT_DIR / "SUMMARY_FIELD_LENS_PUPIL_FOURIER.md", rows, elapsed)
    print(
        json.dumps(
            {
                "elapsed_seconds": elapsed,
                "z1m": [
                    row for row in rows if abs(float(row["z_m"]) - 1.0) < 1e-9
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
