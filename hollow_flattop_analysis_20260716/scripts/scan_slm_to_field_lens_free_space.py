"""Propagate the V2 SLM field through free space and test a 15 mm lens pupil.

The simulation uses the physical SLM grid (1024 x 1024 at 17 um), embeds it
in a zero-padded 2048 x 2048 grid, and propagates with the angular-spectrum
method from z=0 to 1 m in 0.1 m increments.

No Zernike spherical aberration, defocus, or other wavefront error is added.
The 15 mm circular pupil is *measured* at every plane; it is not applied
during propagation.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
ANALYSIS_ROOT = SCRIPT_DIR.parent
POINT2P_ROOT = ANALYSIS_ROOT.parent
PROJECT_ROOT = POINT2P_ROOT / "rtad_mraf_gs_python_test_20260605"
REAL_TEST_ROOT = POINT2P_ROOT / "real_test"
CASE_DIR = (
    PROJECT_ROOT
    / "artifacts"
    / "run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm"
)
INPUT_BGDATA = REAL_TEST_ROOT / "20260716" / "G-光斑-1.bgData"
INPUT_SUMMARY = (
    REAL_TEST_ROOT
    / "20260716"
    / "analysis_G_spot_1"
    / "G-光斑-1_beam_size_summary.json"
)
OUTPUT_DIR = ANALYSIS_ROOT / "results" / "09_slm_to_field_lens_free_space"

for module_root in (PROJECT_ROOT, REAL_TEST_ROOT):
    value = str(module_root)
    if value not in sys.path:
        sys.path.insert(0, value)

from analyze_rect_flattop_size import load_spiricon_frame, robust_background
from convert_to_slm import convert_to_slm


WAVELENGTH_M = 532e-9
SLM_N = 1024
SLM_PITCH_M = 17e-6
PAD_N = 2048
PUPIL_DIAMETER_M = 15e-3
IDEAL_GAUSSIAN_1E2_DIAMETER_M = 6.5e-3
DISTANCES_M = np.linspace(0.0, 1.0, 11)


def normalize_power(amplitude: np.ndarray) -> np.ndarray:
    power = float(np.sum(np.abs(amplitude) ** 2, dtype=np.float64))
    if power <= 0.0 or not np.isfinite(power):
        raise ValueError(f"Invalid power: {power}")
    return (amplitude / math.sqrt(power)).astype(amplitude.dtype, copy=False)


def slm_coordinates() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    axis = (
        np.arange(SLM_N, dtype=np.float32) - SLM_N // 2
    ) * np.float32(SLM_PITCH_M)
    X, Y = np.meshgrid(axis, axis)
    return axis, axis.copy(), X, Y


def measured_amplitude_on_slm(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Measured amplitude from the full original frame without a box crop.

    The signed background-subtracted image is lightly smoothed before
    clipping negative values. This suppresses camera noise without imposing
    the artificial rectangular boundary used by the legacy component crop.
    """
    summary = json.loads(INPUT_SUMMARY.read_text(encoding="utf-8"))
    image, _ = load_spiricon_frame(INPUT_BGDATA)
    background, _ = robust_background(image, corner_px=200)
    signed_signal = image.astype(np.float32) - np.float32(background)
    smoothed_signal = gaussian_filter(
        signed_signal,
        sigma=5.0,
        mode="nearest",
    ).astype(np.float32)
    full_signal = np.clip(smoothed_signal, 0.0, None)
    center_x_px, center_y_px = (
        float(value) for value in summary["moments"]["center_px_global"]
    )
    sx_um = float(summary["pixel_scale_x_um"])
    sy_um = float(summary["pixel_scale_y_um"])
    source_x_px = center_x_px + X * np.float32(1e6 / sx_um)
    source_y_px = center_y_px + Y * np.float32(1e6 / sy_um)
    sampled_intensity = map_coordinates(
        full_signal,
        [source_y_px, source_x_px],
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32)
    return normalize_power(np.sqrt(sampled_intensity).astype(np.float32))


def ideal_gaussian_on_slm(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    w = np.float32(IDEAL_GAUSSIAN_1E2_DIAMETER_M / 2.0)
    amplitude = np.exp(-(X * X + Y * Y) / (w * w)).astype(np.float32)
    return normalize_power(amplitude)


def load_phase_variants() -> dict[str, np.ndarray]:
    phase_refined = np.load(CASE_DIR / "phase_refined.npy").astype(np.float32)
    config = json.loads((CASE_DIR / "config_used.json").read_text(encoding="utf-8"))
    dx_doe_um = float(config["grid"]["dx_doe_m"]) * 1e6
    no_blaze = convert_to_slm(
        phase_refined,
        dx_doe_um=dx_doe_um,
        N=phase_refined.shape[0],
        lambda_m=WAVELENGTH_M,
        f_m=float(config["physical"]["focal_length_m"]),
        focal_dx_um=float(config["grid"]["focal_dx_um"]),
        slm_res=SLM_N,
        slm_pitch_um=SLM_PITCH_M * 1e6,
    ).astype(np.float32)
    return {"v2_no_blaze": no_blaze}


def embed_slm(field_slm: np.ndarray) -> np.ndarray:
    padded = np.zeros((PAD_N, PAD_N), dtype=np.complex64)
    start = (PAD_N - SLM_N) // 2
    padded[start : start + SLM_N, start : start + SLM_N] = field_slm
    return padded


def radial_energy_diameter(
    intensity: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    center_x_m: float,
    center_y_m: float,
    fraction: float,
) -> float:
    radius = np.sqrt((X - center_x_m) ** 2 + (Y - center_y_m) ** 2)
    radial_index = np.floor(radius / SLM_PITCH_M).astype(np.int32)
    energy_by_radius = np.bincount(
        radial_index.ravel(),
        weights=intensity.ravel(),
    )
    cumulative = np.cumsum(energy_by_radius, dtype=np.float64)
    target = fraction * cumulative[-1]
    index = int(np.searchsorted(cumulative, target, side="left"))
    return float(2.0 * (index + 0.5) * SLM_PITCH_M)


def measure_plane(
    intensity: np.ndarray,
    X: np.ndarray,
    Y: np.ndarray,
    z_m: float,
    branch: str,
) -> dict[str, float | str]:
    total = float(np.sum(intensity, dtype=np.float64))
    centroid_x = float(np.sum(intensity * X, dtype=np.float64) / total)
    centroid_y = float(np.sum(intensity * Y, dtype=np.float64) / total)
    pupil_radius = PUPIL_DIAMETER_M / 2.0
    fixed_mask = X * X + Y * Y <= pupil_radius * pupil_radius
    recentered_mask = (
        (X - centroid_x) ** 2 + (Y - centroid_y) ** 2
        <= pupil_radius * pupil_radius
    )
    inside_fixed = float(np.sum(intensity[fixed_mask], dtype=np.float64) / total)
    inside_recentered = float(
        np.sum(intensity[recentered_mask], dtype=np.float64) / total
    )
    return {
        "branch": branch,
        "z_m": z_m,
        "centroid_x_mm": centroid_x * 1e3,
        "centroid_y_mm": centroid_y * 1e3,
        "centroid_radius_mm": math.hypot(centroid_x, centroid_y) * 1e3,
        "power_inside_fixed_15mm": inside_fixed,
        "power_outside_fixed_15mm_percent": 100.0 * (1.0 - inside_fixed),
        "power_inside_recentered_15mm": inside_recentered,
        "power_outside_recentered_15mm_percent": 100.0
        * (1.0 - inside_recentered),
        "diameter_99_percent_fixed_mm": radial_energy_diameter(
            intensity, X, Y, 0.0, 0.0, 0.99
        )
        * 1e3,
        "diameter_99p9_percent_fixed_mm": radial_energy_diameter(
            intensity, X, Y, 0.0, 0.0, 0.999
        )
        * 1e3,
        "diameter_99_percent_recentered_mm": radial_energy_diameter(
            intensity, X, Y, centroid_x, centroid_y, 0.99
        )
        * 1e3,
        "diameter_99p9_percent_recentered_mm": radial_energy_diameter(
            intensity, X, Y, centroid_x, centroid_y, 0.999
        )
        * 1e3,
    }


def propagate_branch(
    field_slm: np.ndarray,
    branch: str,
    X_pad: np.ndarray,
    Y_pad: np.ndarray,
    delta_kz: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[float, np.ndarray]]:
    field0 = embed_slm(field_slm)
    spectrum = np.fft.fft2(field0, norm="ortho").astype(np.complex64)
    rows: list[dict[str, Any]] = []
    snapshots: dict[float, np.ndarray] = {}
    for z_m in DISTANCES_M:
        if z_m == 0.0:
            field = field0
        else:
            transfer = np.exp(1j * delta_kz * float(z_m)).astype(np.complex64)
            field = np.fft.ifft2(spectrum * transfer, norm="ortho").astype(
                np.complex64
            )
        intensity = (np.abs(field) ** 2).astype(np.float32)
        rows.append(measure_plane(intensity, X_pad, Y_pad, float(z_m), branch))
        if any(abs(float(z_m) - value) < 1e-9 for value in (0.0, 0.5, 1.0)):
            snapshots[float(z_m)] = intensity
        print(
            f"{branch}: z={z_m:.1f} m, "
            f"outside15={rows[-1]['power_outside_fixed_15mm_percent']:.5f}%, "
            f"D99.9={rows[-1]['diameter_99p9_percent_fixed_mm']:.3f} mm",
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
        outside = np.asarray(
            [float(row["power_outside_fixed_15mm_percent"]) for row in group]
        )
        d999 = np.asarray(
            [float(row["diameter_99p9_percent_fixed_mm"]) for row in group]
        )
        d99 = np.asarray(
            [float(row["diameter_99_percent_fixed_mm"]) for row in group]
        )
        centroid = np.asarray(
            [float(row["centroid_radius_mm"]) for row in group]
        )
        axes[0, 0].semilogy(z, np.maximum(outside, 1e-7), "o-", label=branch)
        axes[0, 1].plot(z, d999, "o-", label=branch)
        axes[1, 0].plot(z, d99, "o-", label=branch)
        axes[1, 1].plot(z, centroid, "o-", label=branch)
    axes[0, 0].set_title("Power outside 15 mm circular pupil")
    axes[0, 0].set_ylabel("outside power / %")
    axes[0, 1].set_title("Fixed-axis 99.9% energy diameter")
    axes[0, 1].set_ylabel("diameter / mm")
    axes[0, 1].axhline(15.0, color="black", ls=":", label="15 mm pupil")
    axes[1, 0].set_title("Fixed-axis 99% energy diameter")
    axes[1, 0].set_ylabel("diameter / mm")
    axes[1, 0].axhline(15.0, color="black", ls=":", label="15 mm pupil")
    axes[1, 1].set_title("Beam energy centroid displacement")
    axes[1, 1].set_ylabel("centroid radius / mm")
    for ax in axes.ravel():
        ax.set_xlabel("SLM-to-field-lens distance / m")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle(
        "Free-space propagation from physical SLM grid; no Zernike aberration"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_snapshots(
    path: Path,
    snapshots_by_branch: dict[str, dict[float, np.ndarray]],
    axis_pad_mm: np.ndarray,
) -> None:
    branches = list(snapshots_by_branch)
    distances = [0.0, 0.5, 1.0]
    fig, axes = plt.subplots(
        len(branches),
        len(distances),
        figsize=(14, 4.3 * len(branches)),
        constrained_layout=True,
    )
    if len(branches) == 1:
        axes = np.asarray([axes])
    plot_mask = np.abs(axis_pad_mm) <= 10.0
    xa = axis_pad_mm[plot_mask]
    for row_index, branch in enumerate(branches):
        for column, z_m in enumerate(distances):
            ax = axes[row_index, column]
            intensity = snapshots_by_branch[branch][z_m]
            roi = intensity[np.ix_(plot_mask, plot_mask)].astype(np.float64)
            db = 10.0 * np.log10(roi / max(float(np.max(roi)), 1e-30) + 1e-12)
            im = ax.imshow(
                db,
                origin="upper",
                extent=[xa[0], xa[-1], xa[-1], xa[0]],
                cmap="turbo",
                vmin=-40.0,
                vmax=0.0,
                aspect="equal",
            )
            ax.add_patch(
                plt.Circle((0.0, 0.0), 7.5, fill=False, color="white", ls="--")
            )
            ax.set_title(f"{branch}\nz={z_m:.1f} m")
            ax.set_xlabel("x / mm")
            ax.set_ylabel("y / mm")
            fig.colorbar(im, ax=ax, label="relative intensity / dB")
    fig.suptitle("SLM-field free-space intensity; dashed circle = 15 mm pupil")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path: Path, rows: list[dict[str, Any]], elapsed_s: float) -> None:
    branches = list(dict.fromkeys(str(row["branch"]) for row in rows))
    lines = [
        "# SLM 到场镜自由空间传播：15 mm 入瞳检查",
        "",
        "## 模型",
        "",
        "- SLM：1024×1024，17 µm 像素，物理尺寸 17.408 mm。",
        "- 波长：532 nm。",
        "- 输入仅两组：原始拍摄 `G-光斑-1.bgData` 与理想 6.5 mm Gaussian。",
        "- 实测光斑从原始整幅相机数据按物理坐标重采样；没有使用矩形 ROI 硬裁切。",
        "- 原始有符号背景扣除结果先作 5 px Gaussian 平滑，再将负值裁为零，以抑制相机噪声而不制造方形边界。",
        "- 相位：V2，无闪耀；不加入球差、离焦或其他 Zernike 项。",
        "- 传播：2048×2048 零填充角谱法，0–1.0 m，每 0.1 m 取样。",
        "- 15 mm 圆形入瞳仅用于每个传播面的能量统计，不参与前面的传播。",
        f"- 耗时：{elapsed_s:.1f} s。",
        "",
        "## 1 m 处结果",
        "",
        "| 分支 | 固定15mm外功率 | D99 | D99.9 | 质心偏移 |",
        "|---|---:|---:|---:|---:|",
    ]
    for branch in branches:
        row = next(
            item
            for item in rows
            if item["branch"] == branch and abs(float(item["z_m"]) - 1.0) < 1e-9
        )
        lines.append(
            f"| {branch} | "
            f"{float(row['power_outside_fixed_15mm_percent']):.6f}% | "
            f"{float(row['diameter_99_percent_fixed_mm']):.3f} mm | "
            f"{float(row['diameter_99p9_percent_fixed_mm']):.3f} mm | "
            f"{float(row['centroid_radius_mm']):.3f} mm |"
        )
    lines.extend(
        [
            "",
            "## 说明",
            "",
            "- 旧版本图中的方形硬边来自后处理：只保留检测到的矩形包围框，并把框外强制置零；它不是原始 `.bgData` 的边界。",
            "- 本结果已取消该矩形裁切，因此用于判断 15 mm 入瞳的实测分支不再包含这项人为衍射源。",
            "- 主要物理判据是固定场镜轴心的 15 mm 圆孔外功率；D99/D99.9 作为能量包络的辅助指标。",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _, _, X_slm, Y_slm = slm_coordinates()
    measured_amplitude = measured_amplitude_on_slm(X_slm, Y_slm)
    ideal_amplitude = ideal_gaussian_on_slm(X_slm, Y_slm)
    phases = load_phase_variants()

    axis_pad = (
        np.arange(PAD_N, dtype=np.float32) - PAD_N // 2
    ) * np.float32(SLM_PITCH_M)
    X_pad, Y_pad = np.meshgrid(axis_pad, axis_pad)
    frequency = np.fft.fftfreq(PAD_N, d=SLM_PITCH_M)
    FX, FY = np.meshgrid(frequency, frequency)
    k = 2.0 * np.pi / WAVELENGTH_M
    argument = np.maximum(
        0.0,
        1.0 - (WAVELENGTH_M * FX) ** 2 - (WAVELENGTH_M * FY) ** 2,
    )
    # Remove the irrelevant global exp(i*k*z) phase to retain numeric precision.
    delta_kz = (k * (np.sqrt(argument) - 1.0)).astype(np.float32)

    branches = {
        "measured_G_v2_no_blaze": measured_amplitude
        * np.exp(1j * phases["v2_no_blaze"]).astype(np.complex64),
        "ideal_G_v2_no_blaze_reference": ideal_amplitude
        * np.exp(1j * phases["v2_no_blaze"]).astype(np.complex64),
    }
    started = time.perf_counter()
    rows: list[dict[str, Any]] = []
    snapshots_by_branch: dict[str, dict[float, np.ndarray]] = {}
    for branch, field in branches.items():
        branch_rows, snapshots = propagate_branch(
            field.astype(np.complex64),
            branch,
            X_pad,
            Y_pad,
            delta_kz,
        )
        rows.extend(branch_rows)
        snapshots_by_branch[branch] = snapshots
    elapsed = time.perf_counter() - started

    write_csv(OUTPUT_DIR / "slm_to_field_lens_scan.csv", rows)
    payload = {
        "model": {
            "wavelength_m": WAVELENGTH_M,
            "slm_shape": [SLM_N, SLM_N],
            "slm_pitch_um": SLM_PITCH_M * 1e6,
            "slm_active_size_mm": SLM_N * SLM_PITCH_M * 1e3,
            "propagation": "zero-padded angular spectrum",
            "padded_shape": [PAD_N, PAD_N],
            "padded_width_mm": PAD_N * SLM_PITCH_M * 1e3,
            "field_lens_pupil_diameter_mm": PUPIL_DIAMETER_M * 1e3,
            "distance_m": DISTANCES_M.tolist(),
            "zernike_aberration_added": False,
        },
        "sources": {
            "v2_phase": str(CASE_DIR / "phase_refined.npy"),
            "measured_gaussian": str(INPUT_BGDATA),
        },
        "elapsed_seconds": elapsed,
        "rows": rows,
    }
    (OUTPUT_DIR / "slm_to_field_lens_scan.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_trends(OUTPUT_DIR / "slm_to_field_lens_trends.png", rows)
    plot_snapshots(
        OUTPUT_DIR / "slm_to_field_lens_snapshots.png",
        snapshots_by_branch,
        axis_pad.astype(np.float64) * 1e3,
    )
    write_summary(OUTPUT_DIR / "SUMMARY_SLM_TO_FIELD_LENS.md", rows, elapsed)
    print(
        json.dumps(
            {
                "elapsed_seconds": elapsed,
                "output_dir": str(OUTPUT_DIR),
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
