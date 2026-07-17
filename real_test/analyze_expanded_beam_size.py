"""Analyze the expanded input beam saved as Spiricon HDF5 data.

Default target:
  real_test/20260605-1/*.bgData

The script reads the raw /BG_DATA/1/DATA int32 array, uses the WIDTH/HEIGHT
and pixel-scale metadata from RAWFRAME, and reports beam size metrics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, uniform_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LEVEL_E2 = float(np.exp(-2.0))


def _scalar(group: h5py.Group, name: str) -> Any:
    value = group[name][()]
    if getattr(value, "shape", None) == (1,):
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value.item() if hasattr(value, "item") else value


def load_spiricon_frame(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    with h5py.File(path, "r") as h5:
        frame = h5["BG_DATA"]["1"]
        raw = frame["RAWFRAME"]
        width = int(_scalar(raw, "WIDTH"))
        height = int(_scalar(raw, "HEIGHT"))
        data_1d = frame["DATA"][()]
        if data_1d.size != width * height:
            raise ValueError(f"{path.name}: DATA length {data_1d.size} != WIDTH*HEIGHT {width * height}")
        image = data_1d.reshape((height, width)).astype(np.float64)
        meta = {
            "source_file": str(path),
            "width_px": width,
            "height_px": height,
            "raw_pixel_scale_x_um": float(_scalar(raw, "PIXELSCALEXUM")),
            "raw_pixel_scale_y_um": float(_scalar(raw, "PIXELSCALEYUM")),
            "binning_x": int(_scalar(raw, "BINNINGX")),
            "binning_y": int(_scalar(raw, "BINNINGY")),
            "exposure_stamp": float(_scalar(raw, "EXPOSURESTAMP")),
            "gain_stamp": float(_scalar(raw, "GAINSTAMP")),
            "timestamp": _scalar(raw, "TIMESTAMP"),
            "data_dtype": str(data_1d.dtype),
        }
        meta["pixel_scale_x_um"] = meta["raw_pixel_scale_x_um"] * meta["binning_x"]
        meta["pixel_scale_y_um"] = meta["raw_pixel_scale_y_um"] * meta["binning_y"]
    return image, meta


def robust_background(image: np.ndarray, corner_px: int) -> tuple[float, float]:
    ny, nx = image.shape
    c = min(corner_px, ny // 4, nx // 4)
    corners = np.concatenate(
        [
            image[:c, :c].ravel(),
            image[:c, -c:].ravel(),
            image[-c:, :c].ravel(),
            image[-c:, -c:].ravel(),
        ]
    )
    bg = float(np.median(corners))
    mad = float(np.median(np.abs(corners - bg)))
    sigma = 1.4826 * mad
    return bg, sigma


def largest_component(mask: np.ndarray) -> np.ndarray:
    labels, count = ndi.label(mask)
    if count == 0:
        raise RuntimeError("No connected beam region found.")
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    return labels == int(np.argmax(sizes))


def expand_bbox(mask: np.ndarray, margin_px: int, shape: tuple[int, int]) -> tuple[int, int, int, int]:
    ys, xs = np.where(mask)
    if len(xs) == 0:
        raise RuntimeError("Cannot compute bbox for an empty mask.")
    ny, nx = shape
    xmin = max(0, int(xs.min()) - margin_px)
    xmax = min(nx - 1, int(xs.max()) + margin_px)
    ymin = max(0, int(ys.min()) - margin_px)
    ymax = min(ny - 1, int(ys.max()) + margin_px)
    return xmin, ymin, xmax, ymax


def weighted_moments(signal: np.ndarray, sx_um: float, sy_um: float) -> dict[str, Any]:
    yy, xx = np.indices(signal.shape)
    weights = np.clip(signal, 0.0, None)
    total = float(weights.sum())
    if total <= 0:
        raise RuntimeError("No positive signal for moment calculation.")
    cx = float((xx * weights).sum() / total)
    cy = float((yy * weights).sum() / total)
    var_x = float(((xx - cx) ** 2 * weights).sum() / total)
    var_y = float(((yy - cy) ** 2 * weights).sum() / total)
    cov_xy = float(((xx - cx) * (yy - cy) * weights).sum() / total)
    cov = np.array([[var_x, cov_xy], [cov_xy, var_y]], dtype=np.float64)
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    angle_deg = float(np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0])))
    return {
        "center_px": [cx, cy],
        "d4sigma_x_um": float(4.0 * np.sqrt(max(var_x, 0.0)) * sx_um),
        "d4sigma_y_um": float(4.0 * np.sqrt(max(var_y, 0.0)) * sy_um),
        "major_d4sigma_um": float(4.0 * np.sqrt(max(float(eigvals[0]), 0.0)) * sx_um),
        "minor_d4sigma_um": float(4.0 * np.sqrt(max(float(eigvals[1]), 0.0)) * sy_um),
        "major_axis_angle_deg": angle_deg,
    }


def crossing_width(axis_um: np.ndarray, profile: np.ndarray, center_idx: int, level: float) -> dict[str, float]:
    peak = float(profile[center_idx])
    local_peak = float(np.max(profile[max(0, center_idx - 200) : min(len(profile), center_idx + 201)]))
    peak = max(peak, local_peak)
    threshold = level * peak

    left = np.nan
    for i in range(center_idx, 0, -1):
        if profile[i] >= threshold and profile[i - 1] < threshold:
            denom = profile[i] - profile[i - 1]
            left = i - 1 + (threshold - profile[i - 1]) / denom if abs(denom) > 1e-12 else float(i)
            break

    right = np.nan
    for i in range(center_idx, len(profile) - 1):
        if profile[i] >= threshold and profile[i + 1] < threshold:
            denom = profile[i + 1] - profile[i]
            right = i + (threshold - profile[i]) / denom if abs(denom) > 1e-12 else float(i)
            break

    width_px = float(right - left) if np.isfinite(left) and np.isfinite(right) else float("nan")
    width_um = float(width_px * np.median(np.diff(axis_um))) if np.isfinite(width_px) else float("nan")
    return {"left_px": float(left), "right_px": float(right), "width_px": width_px, "width_um": width_um}


def analyze(path: Path, outdir: Path, threshold_frac: float, sigma_px: float, corner_px: int) -> dict[str, Any]:
    image, meta = load_spiricon_frame(path)
    sx_um = float(meta["pixel_scale_x_um"])
    sy_um = float(meta["pixel_scale_y_um"])
    bg, bg_sigma = robust_background(image, corner_px=corner_px)
    signal = np.clip(image - bg, 0.0, None)

    clipped = np.clip(signal, 0.0, np.percentile(signal, 99.9))
    smoothed = gaussian_filter(clipped, sigma=sigma_px)
    component = largest_component(smoothed > threshold_frac * float(smoothed.max()))
    bbox_margin = max(50, int(0.08 * max(image.shape)))
    xmin, ymin, xmax, ymax = expand_bbox(component, bbox_margin, image.shape)
    roi_signal = signal[ymin : ymax + 1, xmin : xmax + 1]
    roi_smoothed = smoothed[ymin : ymax + 1, xmin : xmax + 1]

    moments = weighted_moments(roi_signal, sx_um=sx_um, sy_um=sy_um)
    center_x = moments["center_px"][0] + xmin
    center_y = moments["center_px"][1] + ymin
    moments["center_px_global"] = [center_x, center_y]
    moments["center_um_from_image_origin"] = [center_x * sx_um, center_y * sy_um]

    cxi = int(round(center_x))
    cyi = int(round(center_y))
    band = 15
    prof_x = signal[max(0, cyi - band) : min(signal.shape[0], cyi + band + 1), :].mean(axis=0)
    prof_y = signal[:, max(0, cxi - band) : min(signal.shape[1], cxi + band + 1)].mean(axis=1)
    prof_x = uniform_filter1d(prof_x, 21)
    prof_y = uniform_filter1d(prof_y, 21)
    x_um = np.arange(signal.shape[1], dtype=np.float64) * sx_um
    y_um = np.arange(signal.shape[0], dtype=np.float64) * sy_um

    widths = {
        "x_fwhm": crossing_width(x_um, prof_x, cxi, 0.5),
        "y_fwhm": crossing_width(y_um, prof_y, cyi, 0.5),
        "x_1e2": crossing_width(x_um, prof_x, cxi, LEVEL_E2),
        "y_1e2": crossing_width(y_um, prof_y, cyi, LEVEL_E2),
    }

    result = {
        **meta,
        "background": {"median": bg, "robust_sigma": bg_sigma},
        "preprocessing": {
            "component_threshold_fraction_of_smoothed_peak": threshold_frac,
            "gaussian_sigma_px": sigma_px,
            "corner_px": corner_px,
        },
        "beam_region_bbox_px": [xmin, ymin, xmax - xmin + 1, ymax - ymin + 1],
        "beam_region_bbox_um": [
            xmin * sx_um,
            ymin * sy_um,
            (xmax - xmin + 1) * sx_um,
            (ymax - ymin + 1) * sy_um,
        ],
        "moments": moments,
        "profile_widths": widths,
        "data_percentiles": {
            f"p{p:g}": float(v) for p, v in zip([0, 1, 5, 50, 95, 99, 99.5, 99.9, 100], np.percentile(image, [0, 1, 5, 50, 95, 99, 99.5, 99.9, 100]))
        },
    }

    outdir.mkdir(parents=True, exist_ok=True)
    stem = path.stem
    json_path = outdir / f"{stem}_beam_size_summary.json"
    fig_path = outdir / f"{stem}_beam_size_diagnostics.png"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    preview = np.clip(signal, 0.0, np.percentile(signal, 99.5))
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    axes[0, 0].imshow(preview, cmap="gray", origin="upper")
    axes[0, 0].plot(center_x, center_y, "r+", ms=12, mew=2)
    axes[0, 0].add_patch(plt.Rectangle((xmin, ymin), xmax - xmin + 1, ymax - ymin + 1, fill=False, ec="lime", lw=1.5))
    axes[0, 0].set_title("Background-subtracted preview")
    axes[0, 1].imshow(roi_smoothed, cmap="magma", origin="upper")
    axes[0, 1].set_title("Smoothed beam ROI")

    axes[1, 0].plot(x_um / 1000.0, prof_x, "k-", lw=1.0)
    axes[1, 0].axvline(center_x * sx_um / 1000.0, color="r", ls="--", lw=0.9)
    axes[1, 0].set_title("X profile")
    axes[1, 0].set_xlabel("x / mm")
    axes[1, 0].set_ylabel("background-subtracted intensity")

    axes[1, 1].plot(y_um / 1000.0, prof_y, "k-", lw=1.0)
    axes[1, 1].axvline(center_y * sy_um / 1000.0, color="r", ls="--", lw=0.9)
    axes[1, 1].set_title("Y profile")
    axes[1, 1].set_xlabel("y / mm")
    axes[1, 1].set_ylabel("background-subtracted intensity")
    fig.suptitle(f"Expanded beam size: {path.name}")
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    result["outputs"] = {"summary_json": str(json_path), "diagnostic_png": str(fig_path)}
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def default_input() -> Path:
    data_dir = Path(__file__).resolve().parent / "20260605-1"
    matches = sorted(data_dir.glob("*.bgData"))
    if not matches:
        raise FileNotFoundError(f"No .bgData file found in {data_dir}")
    return matches[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze expanded beam size from .bgData/.bmData HDF5 data.")
    parser.add_argument("--input", type=Path, default=None, help="Input .bgData/.bmData file. Defaults to real_test/20260605-1/*.bgData.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory.")
    parser.add_argument("--threshold-frac", type=float, default=0.05, help="Fraction of smoothed peak used to find the main beam component.")
    parser.add_argument("--sigma-px", type=float, default=5.0, help="Gaussian sigma for beam-region localization.")
    parser.add_argument("--corner-px", type=int, default=200, help="Corner size for robust background estimation.")
    args = parser.parse_args()

    input_path = args.input if args.input is not None else default_input()
    outdir = args.outdir if args.outdir is not None else input_path.parent / "beam_size_analysis"
    result = analyze(input_path, outdir, args.threshold_frac, args.sigma_px, args.corner_px)

    print(f"Input: {input_path}")
    print(f"Output: {outdir}")
    print(f"Center: ({result['moments']['center_px_global'][0]:.1f}, {result['moments']['center_px_global'][1]:.1f}) px")
    print(
        "D4sigma: "
        f"x={result['moments']['d4sigma_x_um']/1000:.3f} mm, "
        f"y={result['moments']['d4sigma_y_um']/1000:.3f} mm"
    )
    print(
        "Profile 1/e^2: "
        f"x={result['profile_widths']['x_1e2']['width_um']/1000:.3f} mm, "
        f"y={result['profile_widths']['y_1e2']['width_um']/1000:.3f} mm"
    )
    print(
        "Profile FWHM: "
        f"x={result['profile_widths']['x_fwhm']['width_um']/1000:.3f} mm, "
        f"y={result['profile_widths']['y_fwhm']['width_um']/1000:.3f} mm"
    )


if __name__ == "__main__":
    main()
