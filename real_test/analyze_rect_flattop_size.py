"""Analyze the rectangular flat-top spot saved as Spiricon HDF5 data.

Default target:
  real_test/20260605-1/20260605-5.bmData

This script is separate from expanded-beam analysis because flat-top size is
measured from gradient-bracketed edges and relative flat levels, not Gaussian
beam moments.
"""

from __future__ import annotations

import argparse
import json
from io import BytesIO
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, uniform_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


LEVEL_E2 = float(np.exp(-2.0))


def _scalar(group: h5py.Group, name: str) -> Any:
    value = group[name][()]
    if getattr(value, "shape", None) == (1,):
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value.item() if hasattr(value, "item") else value


def load_spiricon_frame(path: Path, frame_key: str = "1") -> tuple[np.ndarray, dict[str, Any]]:
    with h5py.File(path, "r") as h5:
        frame = h5["BG_DATA"][frame_key]
        raw = frame["RAWFRAME"]
        width = int(_scalar(raw, "WIDTH"))
        height = int(_scalar(raw, "HEIGHT"))
        data_1d = frame["DATA"][()]
        data_type = _scalar(frame, "DATA_TYPE") if "DATA_TYPE" in frame else "Raw"
        if str(data_type).lower() == "tiff":
            with Image.open(BytesIO(data_1d.tobytes())) as tiff:
                decoded = np.asarray(tiff)
            # BeamGage's compressed float TIFF omits the SampleFormat tag.
            # Pillow consequently exposes the float32 bit pattern as int32.
            if decoded.dtype == np.int32 and np.max(np.abs(decoded.astype(np.int64))) > 10_000_000:
                decoded = decoded.view(np.float32)
            image = np.asarray(decoded, dtype=np.float64)
            if image.shape != (height, width):
                raise ValueError(f"{path.name} frame {frame_key}: decoded TIFF shape {image.shape} != {(height, width)}")
        else:
            if data_1d.size != width * height:
                raise ValueError(f"{path.name}: DATA length {data_1d.size} != WIDTH*HEIGHT {width * height}")
            image = data_1d.reshape((height, width)).astype(np.float64)
        meta = {
            "source_file": str(path),
            "frame_key": str(frame_key),
            "stored_data_type": str(data_type),
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
        raise RuntimeError("No bright flat-top region found.")
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    return labels == int(np.argmax(sizes))


def component_center(mask: np.ndarray, signal: np.ndarray) -> tuple[float, float, tuple[int, int, int, int]]:
    ys, xs = np.where(mask)
    vals = np.clip(signal[mask], 0.0, None)
    total = float(vals.sum())
    if total > 0:
        cx = float((xs * vals).sum() / total)
        cy = float((ys * vals).sum() / total)
    else:
        cx = float(xs.mean())
        cy = float(ys.mean())
    bbox = (int(xs.min()), int(ys.min()), int(xs.max() - xs.min() + 1), int(ys.max() - ys.min() + 1))
    return cx, cy, bbox


def find_edge_pair(gradient: np.ndarray, center: int, search_radius: int) -> tuple[int, int, list[tuple[int, float]], list[tuple[int, float]]]:
    start = max(2, center - search_radius)
    stop = min(len(gradient) - 2, center + search_radius)
    threshold = 0.15 * float(np.max(gradient[start:stop])) if stop > start else 0.15 * float(np.max(gradient))
    peaks: list[tuple[int, float]] = []
    for idx in range(start, stop):
        if (
            gradient[idx] > threshold
            and gradient[idx] >= gradient[idx - 1]
            and gradient[idx] >= gradient[idx - 2]
            and gradient[idx] > gradient[idx + 1]
            and gradient[idx] > gradient[idx + 2]
        ):
            peaks.append((idx, float(gradient[idx])))
    left_peaks = sorted([p for p in peaks if p[0] < center], key=lambda p: -p[1])
    right_peaks = sorted([p for p in peaks if p[0] > center], key=lambda p: -p[1])
    left_edge = left_peaks[0][0] if left_peaks else max(0, center - 25)
    right_edge = right_peaks[0][0] if right_peaks else min(len(gradient) - 1, center + 25)
    return left_edge, right_edge, left_peaks[:5], right_peaks[:5]


def crossing_width(profile: np.ndarray, center: int, bg: float, flat_level: float, frac: float) -> dict[str, float]:
    threshold = bg + frac * (flat_level - bg)
    left = np.nan
    for i in range(center, 0, -1):
        if profile[i] >= threshold and profile[i - 1] < threshold:
            denom = profile[i] - profile[i - 1]
            left = i - 1 + (threshold - profile[i - 1]) / denom if abs(denom) > 1e-12 else float(i)
            break
    right = np.nan
    for i in range(center, len(profile) - 1):
        if profile[i] >= threshold and profile[i + 1] < threshold:
            denom = profile[i + 1] - profile[i]
            right = i + (threshold - profile[i]) / denom if abs(denom) > 1e-12 else float(i)
            break
    width = float(right - left) if np.isfinite(left) and np.isfinite(right) else float("nan")
    return {"left_px": float(left), "right_px": float(right), "width_px": width, "threshold": float(threshold)}


def analyze(path: Path, outdir: Path, threshold_frac: float, sigma_px: float, corner_px: int) -> dict[str, Any]:
    image, meta = load_spiricon_frame(path)
    sx_um = float(meta["pixel_scale_x_um"])
    sy_um = float(meta["pixel_scale_y_um"])
    bg, bg_sigma = robust_background(image, corner_px=corner_px)
    signal = np.clip(image - bg, 0.0, None)
    locate = np.clip(signal, 0.0, np.percentile(signal, 99.9))
    smoothed = gaussian_filter(locate, sigma=sigma_px)
    reference = float(np.percentile(smoothed, 99.9))
    component = largest_component(smoothed > threshold_frac * reference)
    cx, cy, bbox = component_center(component, smoothed)

    cxi = int(round(cx))
    cyi = int(round(cy))
    band = 2
    prof_x = smoothed[max(0, cyi - band) : min(smoothed.shape[0], cyi + band + 1), :].mean(axis=0)
    prof_y = smoothed[:, max(0, cxi - band) : min(smoothed.shape[1], cxi + band + 1)].mean(axis=1)
    prof_x = uniform_filter1d(prof_x, 7)
    prof_y = uniform_filter1d(prof_y, 7)
    gx = np.abs(np.gradient(prof_x))
    gy = np.abs(np.gradient(prof_y))
    search_radius = 180
    lx, rx, left_x_peaks, right_x_peaks = find_edge_pair(gx, cxi, search_radius)
    ly, ry, left_y_peaks, right_y_peaks = find_edge_pair(gy, cyi, search_radius)

    margin = 4
    flat_x = prof_x[max(0, lx + margin) : min(len(prof_x), rx - margin)]
    flat_y = prof_y[max(0, ly + margin) : min(len(prof_y), ry - margin)]
    flat_level_x = float(np.median(flat_x)) if len(flat_x) else float(np.median(prof_x[component[cyi, :]]))
    flat_level_y = float(np.median(flat_y)) if len(flat_y) else float(np.median(prof_y[component[:, cxi]]))
    bg_x = float(np.percentile(prof_x, 5))
    bg_y = float(np.percentile(prof_y, 5))

    x90 = crossing_width(prof_x, cxi, bg_x, flat_level_x, 0.9)
    x50 = crossing_width(prof_x, cxi, bg_x, flat_level_x, 0.5)
    x13 = crossing_width(prof_x, cxi, bg_x, flat_level_x, LEVEL_E2)
    y90 = crossing_width(prof_y, cyi, bg_y, flat_level_y, 0.9)
    y50 = crossing_width(prof_y, cyi, bg_y, flat_level_y, 0.5)
    y13 = crossing_width(prof_y, cyi, bg_y, flat_level_y, LEVEL_E2)

    widths = {
        "size90": {"x_px": x90["width_px"], "y_px": y90["width_px"], "x_um": x90["width_px"] * sx_um, "y_um": y90["width_px"] * sy_um},
        "size50": {"x_px": x50["width_px"], "y_px": y50["width_px"], "x_um": x50["width_px"] * sx_um, "y_um": y50["width_px"] * sy_um},
        "size13p5": {"x_px": x13["width_px"], "y_px": y13["width_px"], "x_um": x13["width_px"] * sx_um, "y_um": y13["width_px"] * sy_um},
    }
    transition_x_um = 0.5 * (widths["size13p5"]["x_um"] - widths["size90"]["x_um"])
    transition_y_um = 0.5 * (widths["size13p5"]["y_um"] - widths["size90"]["y_um"])

    flat_mask = np.zeros_like(component, dtype=bool)
    fx0 = max(0, int(np.floor(x90["left_px"])))
    fx1 = min(image.shape[1], int(np.ceil(x90["right_px"])))
    fy0 = max(0, int(np.floor(y90["left_px"])))
    fy1 = min(image.shape[0], int(np.ceil(y90["right_px"])))
    if fx1 > fx0 and fy1 > fy0:
        flat_mask[fy0:fy1, fx0:fx1] = True
    flat_vals = signal[flat_mask]
    flat_mean = float(np.mean(flat_vals)) if flat_vals.size else float("nan")
    flat_std = float(np.std(flat_vals)) if flat_vals.size else float("nan")

    result = {
        **meta,
        "background": {"median": bg, "robust_sigma": bg_sigma},
        "preprocessing": {
            "component_threshold_fraction_of_p99_9_smoothed": threshold_frac,
            "gaussian_sigma_px": sigma_px,
            "corner_px": corner_px,
        },
        "flat_center_px": [cx, cy],
        "bright_component_bbox_px": list(bbox),
        "gradient_edges_x_px": [lx, rx],
        "gradient_edges_y_px": [ly, ry],
        "gradient_edges_x_um": [lx * sx_um, rx * sx_um],
        "gradient_edges_y_um": [ly * sy_um, ry * sy_um],
        "flat_level_x": flat_level_x,
        "flat_level_y": flat_level_y,
        "profile_background_x": bg_x,
        "profile_background_y": bg_y,
        "edge_peak_candidates": {
            "x_left": left_x_peaks,
            "x_right": right_x_peaks,
            "y_left": left_y_peaks,
            "y_right": right_y_peaks,
        },
        "crossings": {"x90": x90, "x50": x50, "x13p5": x13, "y90": y90, "y50": y50, "y13p5": y13},
        "widths": widths,
        "transition_width_um": {"x": float(transition_x_um), "y": float(transition_y_um)},
        "aspect_ratio_size50": float(widths["size50"]["x_um"] / widths["size50"]["y_um"]) if widths["size50"]["y_um"] > 0 else float("nan"),
        "flat_region_stats_from_size90_box": {
            "pixel_count": int(flat_vals.size),
            "mean": flat_mean,
            "std": flat_std,
            "rms_percent": float(flat_std / flat_mean * 100.0) if flat_mean > 0 else float("nan"),
        },
        "data_percentiles": {
            f"p{p:g}": float(v) for p, v in zip([0, 1, 5, 50, 95, 99, 99.5, 99.9, 100], np.percentile(image, [0, 1, 5, 50, 95, 99, 99.5, 99.9, 100]))
        },
    }

    outdir.mkdir(parents=True, exist_ok=True)
    stem = path.stem
    json_path = outdir / f"{stem}_flattop_size_summary.json"
    fig_path = outdir / f"{stem}_flattop_size_diagnostics.png"

    preview = np.clip(signal, 0.0, np.percentile(signal, 99.5))
    x_um = np.arange(image.shape[1]) * sx_um
    y_um = np.arange(image.shape[0]) * sy_um
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    axes[0, 0].imshow(preview, cmap="gray", origin="upper")
    axes[0, 0].plot(cx, cy, "r+", ms=12, mew=2)
    axes[0, 0].add_patch(plt.Rectangle((bbox[0], bbox[1]), bbox[2], bbox[3], fill=False, ec="lime", lw=1.5))
    axes[0, 0].set_title("Background-subtracted preview")
    zoom = 120
    axes[0, 1].imshow(preview, cmap="gray", origin="upper")
    axes[0, 1].plot(cx, cy, "r+", ms=12, mew=2)
    axes[0, 1].axvline(lx, color="orange", ls=":", lw=1)
    axes[0, 1].axvline(rx, color="orange", ls=":", lw=1)
    axes[0, 1].axhline(ly, color="orange", ls=":", lw=1)
    axes[0, 1].axhline(ry, color="orange", ls=":", lw=1)
    axes[0, 1].set_xlim(max(0, cx - zoom), min(image.shape[1], cx + zoom))
    axes[0, 1].set_ylim(min(image.shape[0], cy + zoom), max(0, cy - zoom))
    axes[0, 1].set_title("Flat-top zoom and gradient edges")

    axes[1, 0].plot(x_um, prof_x, "k-", lw=1.0)
    for item, color in [(x90, "tab:green"), (x50, "tab:red"), (x13, "tab:blue")]:
        axes[1, 0].axvline(item["left_px"] * sx_um, color=color, ls="--", lw=0.9)
        axes[1, 0].axvline(item["right_px"] * sx_um, color=color, ls="--", lw=0.9)
    axes[1, 0].axvline(lx * sx_um, color="orange", ls=":", lw=1.2)
    axes[1, 0].axvline(rx * sx_um, color="orange", ls=":", lw=1.2)
    axes[1, 0].set_xlim((cx - zoom) * sx_um, (cx + zoom) * sx_um)
    axes[1, 0].set_title("X profile")
    axes[1, 0].set_xlabel("x / um")

    axes[1, 1].plot(y_um, prof_y, "k-", lw=1.0)
    for item, color in [(y90, "tab:green"), (y50, "tab:red"), (y13, "tab:blue")]:
        axes[1, 1].axvline(item["left_px"] * sy_um, color=color, ls="--", lw=0.9)
        axes[1, 1].axvline(item["right_px"] * sy_um, color=color, ls="--", lw=0.9)
    axes[1, 1].axvline(ly * sy_um, color="orange", ls=":", lw=1.2)
    axes[1, 1].axvline(ry * sy_um, color="orange", ls=":", lw=1.2)
    axes[1, 1].set_xlim((cy - zoom) * sy_um, (cy + zoom) * sy_um)
    axes[1, 1].set_title("Y profile")
    axes[1, 1].set_xlabel("y / um")
    fig.suptitle(f"Rectangular flat-top size: {path.name}")
    fig.savefig(fig_path, dpi=180)
    plt.close(fig)

    result["outputs"] = {"summary_json": str(json_path), "diagnostic_png": str(fig_path)}
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def default_input() -> Path:
    path = Path(__file__).resolve().parent / "20260605-1" / "20260605-5.bmData"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze rectangular flat-top size from .bmData HDF5 data.")
    parser.add_argument("--input", type=Path, default=None, help="Input .bmData file. Defaults to real_test/20260605-1/20260605-5.bmData.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory.")
    parser.add_argument("--threshold-frac", type=float, default=0.10, help="Fraction of p99.9 smoothed signal used to find the main flat-top component.")
    parser.add_argument("--sigma-px", type=float, default=1.5, help="Gaussian sigma for bright-region localization.")
    parser.add_argument("--corner-px", type=int, default=30, help="Corner size for robust background estimation.")
    args = parser.parse_args()

    input_path = args.input if args.input is not None else default_input()
    outdir = args.outdir if args.outdir is not None else input_path.parent / "flattop_size_analysis"
    result = analyze(input_path, outdir, args.threshold_frac, args.sigma_px, args.corner_px)

    print(f"Input: {input_path}")
    print(f"Output: {outdir}")
    print(f"Center: ({result['flat_center_px'][0]:.1f}, {result['flat_center_px'][1]:.1f}) px")
    print(
        "size50: "
        f"x={result['widths']['size50']['x_um']:.1f} um, "
        f"y={result['widths']['size50']['y_um']:.1f} um, "
        f"aspect={result['aspect_ratio_size50']:.3f}"
    )
    print(
        "size90: "
        f"x={result['widths']['size90']['x_um']:.1f} um, "
        f"y={result['widths']['size90']['y_um']:.1f} um"
    )
    print(
        "size13.5: "
        f"x={result['widths']['size13p5']['x_um']:.1f} um, "
        f"y={result['widths']['size13p5']['y_um']:.1f} um"
    )


if __name__ == "__main__":
    main()
