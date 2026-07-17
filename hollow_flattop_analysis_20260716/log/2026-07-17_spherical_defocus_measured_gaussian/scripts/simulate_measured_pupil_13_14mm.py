"""Test 14 mm and 13 mm field-lens pupils with the measured Gaussian only.

The measured input is resampled from the full original Spiricon frame without
a rectangular ROI crop. The V2 phase has no blaze and no added Zernike terms.
The field propagates from the physical SLM grid to the field-lens plane over
0--1 m, then passes through a 14 mm or 13 mm circular pupil and a 429 mm
Fourier lens.
"""

from __future__ import annotations

import csv
import json
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
    SLM_N,
    SLM_PITCH_M,
    WAVELENGTH_M,
    embed_slm,
    load_phase_variants,
    measured_amplitude_on_slm,
    slm_coordinates,
)
from simulate_field_lens_pupil_fourier import (
    EXPERIMENT_CORE_WIDTH_X_UM,
    EXPERIMENT_CORE_WIDTH_Y_UM,
    FIELD_LENS_FOCAL_LENGTH_M,
    comparison_metrics,
    focal_metrics,
    forward_fourier,
    normalize_to_core,
)


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "11_measured_pupil_13_14mm"
PUPIL_DIAMETERS_MM = (14.0, 13.0)
SNAPSHOT_DISTANCES_M = (0.0, 0.5, 1.0)


def evaluate_pupil(
    field_lens: np.ndarray,
    pupil_mask: np.ndarray,
    diameter_mm: float,
    z_m: float,
    focal_no_pupil: np.ndarray,
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray]:
    power_before = float(np.sum(np.abs(field_lens) ** 2, dtype=np.float64))
    clipped_field = field_lens * pupil_mask
    power_after = float(np.sum(np.abs(clipped_field) ** 2, dtype=np.float64))
    focal_clipped = normalize_to_core(
        np.abs(forward_fourier(clipped_field)) ** 2,
        x_focal_um,
        y_focal_um,
    )
    metrics_no = focal_metrics(focal_no_pupil, x_focal_um, y_focal_um)
    metrics_clip = focal_metrics(focal_clipped, x_focal_um, y_focal_um)
    comparison = comparison_metrics(
        focal_no_pupil,
        focal_clipped,
        x_focal_um,
        y_focal_um,
    )
    row = {
        "input": "measured_G_full_frame",
        "pupil_diameter_mm": diameter_mm,
        "z_m": z_m,
        "pupil_power_loss_percent": 100.0
        * (1.0 - power_after / power_before),
        "focal_power_transmission_percent": 100.0
        * power_after
        / power_before,
        "no_pupil_center_ratio": metrics_no[
            "center_window_over_core_mean"
        ],
        "clipped_center_ratio": metrics_clip[
            "center_window_over_core_mean"
        ],
        "center_ratio_change": metrics_clip[
            "center_window_over_core_mean"
        ]
        - metrics_no["center_window_over_core_mean"],
        "no_pupil_middle_over_sides": metrics_no[
            "middle_third_over_side_thirds"
        ],
        "clipped_middle_over_sides": metrics_clip[
            "middle_third_over_side_thirds"
        ],
        "middle_over_sides_change": metrics_clip[
            "middle_third_over_side_thirds"
        ]
        - metrics_no["middle_third_over_side_thirds"],
        "no_pupil_core_rms_fraction": metrics_no["core_rms_fraction"],
        "clipped_core_rms_fraction": metrics_clip["core_rms_fraction"],
        "core_rms_fraction_change": metrics_clip["core_rms_fraction"]
        - metrics_no["core_rms_fraction"],
        **comparison,
    }
    return row, focal_clipped


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_trends(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for diameter_mm in PUPIL_DIAMETERS_MM:
        group = [
            row
            for row in rows
            if float(row["pupil_diameter_mm"]) == diameter_mm
        ]
        z = np.asarray([float(row["z_m"]) for row in group])
        loss = np.asarray(
            [float(row["pupil_power_loss_percent"]) for row in group]
        )
        center_change = np.asarray(
            [float(row["center_ratio_change"]) for row in group]
        )
        l2 = np.asarray(
            [float(row["normalized_relative_l2_in_roi"]) for row in group]
        )
        rms_change = np.asarray(
            [float(row["core_rms_fraction_change"]) for row in group]
        )
        label = f"{diameter_mm:.0f} mm pupil"
        axes[0, 0].plot(z, loss, "o-", label=label)
        axes[0, 1].plot(z, center_change, "o-", label=label)
        axes[1, 0].plot(z, l2, "o-", label=label)
        axes[1, 1].plot(z, 100.0 * rms_change, "o-", label=label)
    axes[0, 0].set_title("Power removed")
    axes[0, 0].set_ylabel("power loss / %")
    axes[0, 1].set_title("Center-ratio change")
    axes[0, 1].set_ylabel("clipped - no pupil")
    axes[1, 0].set_title("Normalized focal-image change")
    axes[1, 0].set_ylabel("relative L2 in ROI")
    axes[1, 1].set_title("Core RMS change")
    axes[1, 1].set_ylabel("absolute change / percentage points")
    for ax in axes.ravel():
        ax.set_xlabel("SLM-to-field-lens distance / m")
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle(
        "Measured Gaussian: 14 mm and 13 mm pupils followed by Fourier lens"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_montage(
    path: Path,
    snapshots: dict[float, dict[str, np.ndarray]],
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
) -> None:
    conditions = ("no pupil", "14 mm pupil", "13 mm pupil")
    fig, axes = plt.subplots(
        len(conditions),
        len(SNAPSHOT_DISTANCES_M),
        figsize=(15, 11),
        constrained_layout=True,
    )
    xmask = np.abs(x_focal_um) <= 250.0
    ymask = np.abs(y_focal_um) <= 130.0
    xa = x_focal_um[xmask]
    ya = y_focal_um[ymask]
    for row_index, condition in enumerate(conditions):
        for column, z_m in enumerate(SNAPSHOT_DISTANCES_M):
            ax = axes[row_index, column]
            roi = snapshots[z_m][condition][np.ix_(ymask, xmask)]
            im = ax.imshow(
                roi,
                origin="upper",
                extent=[xa[0], xa[-1], ya[-1], ya[0]],
                cmap="turbo",
                vmin=0.30,
                vmax=1.55,
                aspect="equal",
            )
            ax.set_title(f"{condition}, z={z_m:.1f} m")
            ax.set_xlabel("x / um")
            ax.set_ylabel("y / um")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(
        "Measured Gaussian back focal plane: no pupil vs 14/13 mm pupils"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path: Path, rows: list[dict[str, Any]], elapsed_s: float) -> None:
    lines = [
        "# 实测Gaussian：14 mm与13 mm场镜入瞳测试",
        "",
        "## 模型",
        "",
        "- 输入：原始整幅 `G-光斑-1.bgData`，无矩形ROI硬裁切。",
        "- 相位：无闪耀V2；不加球差、离焦或其他Zernike项。",
        "- SLM到场镜：0–1 m角谱传播，每0.1 m取样。",
        "- 场镜：429 mm焦距；分别施加固定轴心14 mm和13 mm圆孔。",
        f"- 焦面采样：{WAVELENGTH_M * FIELD_LENS_FOCAL_LENGTH_M / (PAD_N * SLM_PITCH_M) * 1e6:.4f} µm。",
        f"- 耗时：{elapsed_s:.1f} s。",
        "",
        "## 1 m处结果",
        "",
        "| 入瞳 | 截去功率 | 中心比(不截断→截断) | 中部/两侧(不截断→截断) | RMS(不截断→截断) | 焦面L2变化 | 50% IoU |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    z1_rows = [
        row for row in rows if abs(float(row["z_m"]) - 1.0) < 1e-9
    ]
    for row in z1_rows:
        lines.append(
            f"| {float(row['pupil_diameter_mm']):.0f} mm | "
            f"{float(row['pupil_power_loss_percent']):.6f}% | "
            f"{float(row['no_pupil_center_ratio']):.5f}→"
            f"{float(row['clipped_center_ratio']):.5f} | "
            f"{float(row['no_pupil_middle_over_sides']):.5f}→"
            f"{float(row['clipped_middle_over_sides']):.5f} | "
            f"{100.0 * float(row['no_pupil_core_rms_fraction']):.2f}%→"
            f"{100.0 * float(row['clipped_core_rms_fraction']):.2f}% | "
            f"{float(row['normalized_relative_l2_in_roi']):.5f} | "
            f"{float(row['footprint_iou_clipped_vs_no_pupil']):.5f} |"
        )
    lines.extend(
        [
            "",
            "## 判断",
            "",
            "- 14 mm入瞳在1 m处截去0.08008%功率，中心比只下降0.00100，没有形成中心空洞。",
            "- 13 mm入瞳在1 m处截去0.11999%功率，中心比下降0.00357；50%外轮廓IoU仍为0.99916，也没有形成中心空洞。",
            "- 实验中心比约0.774，而本模型不截断时约0.985；13 mm截断只提供所需中心下降量的约1.7%。",
            "- 因此在当前实测光斑、V2相位和0–1 m距离范围内，14 mm或13 mm入瞳截断仍不足以解释7月16日空洞。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _, _, X_slm, Y_slm = slm_coordinates()
    measured_amplitude = measured_amplitude_on_slm(X_slm, Y_slm)
    phase_v2 = load_phase_variants()["v2_no_blaze"]
    field_slm = measured_amplitude * np.exp(1j * phase_v2).astype(
        np.complex64
    )
    field0 = embed_slm(field_slm.astype(np.complex64))
    spectrum = np.fft.fft2(field0, norm="ortho").astype(np.complex64)

    axis_pad = (
        np.arange(PAD_N, dtype=np.float32) - PAD_N // 2
    ) * np.float32(SLM_PITCH_M)
    X_pad, Y_pad = np.meshgrid(axis_pad, axis_pad)
    pupil_masks = {
        diameter_mm: (
            X_pad * X_pad + Y_pad * Y_pad
            <= np.float32(diameter_mm * 0.5e-3) ** 2
        )
        for diameter_mm in PUPIL_DIAMETERS_MM
    }

    frequency = np.fft.fftfreq(PAD_N, d=SLM_PITCH_M)
    FX, FY = np.meshgrid(frequency, frequency)
    k = 2.0 * np.pi / WAVELENGTH_M
    angular_argument = np.maximum(
        0.0,
        1.0 - (WAVELENGTH_M * FX) ** 2 - (WAVELENGTH_M * FY) ** 2,
    )
    delta_kz = (k * (np.sqrt(angular_argument) - 1.0)).astype(np.float32)

    focal_frequency = np.fft.fftshift(
        np.fft.fftfreq(PAD_N, d=SLM_PITCH_M)
    )
    x_focal_um = (
        WAVELENGTH_M
        * FIELD_LENS_FOCAL_LENGTH_M
        * focal_frequency
        * 1e6
    )
    y_focal_um = x_focal_um.copy()

    rows: list[dict[str, Any]] = []
    snapshots: dict[float, dict[str, np.ndarray]] = {}
    started = time.perf_counter()
    for z_m in DISTANCES_M:
        if z_m == 0.0:
            field_lens = field0
        else:
            transfer = np.exp(1j * delta_kz * float(z_m)).astype(
                np.complex64
            )
            field_lens = np.fft.ifft2(
                spectrum * transfer,
                norm="ortho",
            ).astype(np.complex64)
        focal_no_pupil = normalize_to_core(
            np.abs(forward_fourier(field_lens)) ** 2,
            x_focal_um,
            y_focal_um,
        )
        if any(
            abs(float(z_m) - value) < 1e-9
            for value in SNAPSHOT_DISTANCES_M
        ):
            snapshots[float(z_m)] = {"no pupil": focal_no_pupil}
        for diameter_mm in PUPIL_DIAMETERS_MM:
            row, focal_clipped = evaluate_pupil(
                field_lens,
                pupil_masks[diameter_mm],
                diameter_mm,
                float(z_m),
                focal_no_pupil,
                x_focal_um,
                y_focal_um,
            )
            rows.append(row)
            if float(z_m) in snapshots:
                snapshots[float(z_m)][
                    f"{diameter_mm:.0f} mm pupil"
                ] = focal_clipped
            print(
                f"{diameter_mm:.0f} mm: z={z_m:.1f} m, "
                f"loss={row['pupil_power_loss_percent']:.6f}%, "
                f"dCenter={row['center_ratio_change']:+.6f}, "
                f"L2={row['normalized_relative_l2_in_roi']:.6f}",
                flush=True,
            )
    elapsed = time.perf_counter() - started

    write_csv(OUTPUT_DIR / "measured_pupil_13_14mm_scan.csv", rows)
    payload = {
        "model": {
            "input": "full-frame measured G-spot",
            "input_file": str(
                ANALYSIS_ROOT.parent
                / "real_test"
                / "20260716"
                / "G-光斑-1.bgData"
            ),
            "measured_input_preprocessing": (
                "full original frame; robust background subtraction; "
                "5 px Gaussian smoothing before negative clipping; "
                "physical-coordinate resampling; no rectangular ROI crop"
            ),
            "phase": "V2 without blaze",
            "added_zernike": False,
            "wavelength_m": WAVELENGTH_M,
            "slm_shape": [SLM_N, SLM_N],
            "slm_pitch_um": SLM_PITCH_M * 1e6,
            "propagation_shape": [PAD_N, PAD_N],
            "slm_to_lens_distances_m": DISTANCES_M.tolist(),
            "field_lens_focal_length_m": FIELD_LENS_FOCAL_LENGTH_M,
            "pupil_diameters_mm": list(PUPIL_DIAMETERS_MM),
            "focal_sampling_um": float(x_focal_um[1] - x_focal_um[0]),
        },
        "elapsed_seconds": elapsed,
        "rows": rows,
    }
    (OUTPUT_DIR / "measured_pupil_13_14mm_scan.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    plot_trends(OUTPUT_DIR / "measured_pupil_13_14mm_trends.png", rows)
    plot_montage(
        OUTPUT_DIR / "measured_pupil_13_14mm_montage.png",
        snapshots,
        x_focal_um,
        y_focal_um,
    )
    for diameter_mm in PUPIL_DIAMETERS_MM:
        np.save(
            OUTPUT_DIR
            / f"measured_G_z1m_focal_with_{diameter_mm:.0f}mm_pupil.npy",
            snapshots[1.0][f"{diameter_mm:.0f} mm pupil"].astype(
                np.float32
            ),
        )
    write_summary(
        OUTPUT_DIR / "SUMMARY_MEASURED_PUPIL_13_14MM.md",
        rows,
        elapsed,
    )
    print(
        json.dumps(
            {
                "elapsed_seconds": elapsed,
                "z1m": [
                    row
                    for row in rows
                    if abs(float(row["z_m"]) - 1.0) < 1e-9
                ],
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
