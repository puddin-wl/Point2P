"""Make readable visualization figures for the real-test beam data.

This script creates colorbar-based figures that are easier to inspect than the
plain extracted PNG previews:
  - rectangular flat-top crop normalized by its flat-region median
  - expanded beam pseudocolor and log-intensity views
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "20260605-1"
OUT_DIR = DATA_DIR / "visualizations"


def _scalar(group: h5py.Group, name: str) -> Any:
    value = group[name][()]
    if getattr(value, "shape", None) == (1,):
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value.item() if hasattr(value, "item") else value


def load_frame(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    with h5py.File(path, "r") as h5:
        frame = h5["BG_DATA"]["1"]
        raw = frame["RAWFRAME"]
        width = int(_scalar(raw, "WIDTH"))
        height = int(_scalar(raw, "HEIGHT"))
        data = frame["DATA"][()]
        if data.size != width * height:
            raise ValueError(f"{path.name}: DATA length {data.size} != WIDTH*HEIGHT {width * height}")
        bin_x = int(_scalar(raw, "BINNINGX"))
        bin_y = int(_scalar(raw, "BINNINGY"))
        raw_sx = float(_scalar(raw, "PIXELSCALEXUM"))
        raw_sy = float(_scalar(raw, "PIXELSCALEYUM"))
        meta = {
            "width_px": width,
            "height_px": height,
            "raw_pixel_scale_x_um": raw_sx,
            "raw_pixel_scale_y_um": raw_sy,
            "binning_x": bin_x,
            "binning_y": bin_y,
            "pixel_scale_x_um": raw_sx * bin_x,
            "pixel_scale_y_um": raw_sy * bin_y,
            "exposure_stamp": float(_scalar(raw, "EXPOSURESTAMP")),
            "gain_stamp": float(_scalar(raw, "GAINSTAMP")),
        }
    return data.reshape((height, width)).astype(np.float64), meta


def robust_background(image: np.ndarray, corner_px: int) -> float:
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
    return float(np.median(corners))


def crop_bounds(left: float, right: float, top: float, bottom: float, shape: tuple[int, int], margin_px: int) -> tuple[int, int, int, int]:
    ny, nx = shape
    x0 = max(0, int(np.floor(left)) - margin_px)
    x1 = min(nx, int(np.ceil(right)) + margin_px)
    y0 = max(0, int(np.floor(top)) - margin_px)
    y1 = min(ny, int(np.ceil(bottom)) + margin_px)
    return x0, x1, y0, y1


def add_colorbar(fig: plt.Figure, ax: plt.Axes, image: Any, label: str) -> None:
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)


def make_flattop_visualization() -> Path:
    src = DATA_DIR / "20260605-5.bmData"
    summary_path = DATA_DIR / "flattop_size_analysis" / "20260605-5_flattop_size_summary.json"
    image, meta = load_frame(src)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    bg = robust_background(image, corner_px=30)
    signal = np.clip(image - bg, 0.0, None)

    x13 = summary["crossings"]["x13p5"]
    y13 = summary["crossings"]["y13p5"]
    x90 = summary["crossings"]["x90"]
    y90 = summary["crossings"]["y90"]
    x0, x1, y0, y1 = crop_bounds(x13["left_px"], x13["right_px"], y13["left_px"], y13["right_px"], signal.shape, margin_px=18)

    crop = signal[y0:y1, x0:x1]
    flat_x0 = max(0, int(np.floor(x90["left_px"])) - x0)
    flat_x1 = min(crop.shape[1], int(np.ceil(x90["right_px"])) - x0)
    flat_y0 = max(0, int(np.floor(y90["left_px"])) - y0)
    flat_y1 = min(crop.shape[0], int(np.ceil(y90["right_px"])) - y0)
    flat_vals = crop[flat_y0:flat_y1, flat_x0:flat_x1]
    flat_level = float(np.median(flat_vals[flat_vals > 0])) if np.any(flat_vals > 0) else float(np.median(crop))
    norm = crop / flat_level if flat_level > 0 else crop
    norm_smooth = gaussian_filter(norm, sigma=0.8)

    sx = float(meta["pixel_scale_x_um"])
    sy = float(meta["pixel_scale_y_um"])
    extent_um = [x0 * sx, x1 * sx, y1 * sy, y0 * sy]
    center = summary["flat_center_px"]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "20260605-5_flattop_crop_uniformity.png"
    profile_out = OUT_DIR / "20260605-5_flattop_raw_profiles.png"

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), constrained_layout=True)

    im0 = axes[0].imshow(crop, cmap="turbo", origin="upper", extent=extent_um, vmin=np.percentile(crop, 2), vmax=np.percentile(crop, 99.5))
    axes[0].set_title("Flat-top ROI, background-subtracted")
    axes[0].set_xlabel("x / um")
    axes[0].set_ylabel("y / um")
    add_colorbar(fig, axes[0], im0, "intensity")

    im1 = axes[1].imshow(norm, cmap="turbo", origin="upper", extent=extent_um, vmin=0.55, vmax=1.35)
    axes[1].set_title("Relative intensity: I / median(flat)")
    axes[1].set_xlabel("x / um")
    axes[1].plot(center[0] * sx, center[1] * sy, "k+", ms=10, mew=1.6)
    add_colorbar(fig, axes[1], im1, "relative intensity")

    im2 = axes[2].imshow(norm_smooth, cmap="RdBu_r", origin="upper", extent=extent_um, vmin=0.75, vmax=1.25)
    axes[2].set_title("Uniformity deviation view")
    axes[2].set_xlabel("x / um")
    add_colorbar(fig, axes[2], im2, "I / median(flat)")

    for ax in axes:
        ax.axvline(x90["left_px"] * sx, color="white", ls="--", lw=1)
        ax.axvline(x90["right_px"] * sx, color="white", ls="--", lw=1)
        ax.axhline(y90["left_px"] * sy, color="white", ls="--", lw=1)
        ax.axhline(y90["right_px"] * sy, color="white", ls="--", lw=1)
        ax.set_aspect("equal")

    fig.suptitle(
        "20260605-5 rectangular flat-top crop | "
        f"effective pixel = {sx:.2f} x {sy:.2f} um | "
        f"size50 = {summary['widths']['size50']['x_um']:.1f} x {summary['widths']['size50']['y_um']:.1f} um"
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)

    make_flattop_raw_profiles(
        signal=signal,
        summary=summary,
        sx=sx,
        sy=sy,
        flat_level=flat_level,
        out=profile_out,
    )
    return out


def make_flattop_raw_profiles(signal: np.ndarray, summary: dict[str, Any], sx: float, sy: float, flat_level: float, out: Path) -> None:
    center_x, center_y = summary["flat_center_px"]
    cx = int(round(center_x))
    cy = int(round(center_y))
    band_half_width = 0

    if band_half_width == 0:
        prof_x = signal[cy, :].astype(np.float64)
        prof_y = signal[:, cx].astype(np.float64)
        profile_label = "single-pixel raw profile"
    else:
        prof_x = signal[max(0, cy - band_half_width) : min(signal.shape[0], cy + band_half_width + 1), :].mean(axis=0)
        prof_y = signal[:, max(0, cx - band_half_width) : min(signal.shape[1], cx + band_half_width + 1)].mean(axis=1)
        profile_label = f"{2 * band_half_width + 1}-pixel band mean"

    prof_x_norm = prof_x / flat_level
    prof_y_norm = prof_y / flat_level
    prof_x_light = uniform_filter1d(prof_x_norm, 3)
    prof_y_light = uniform_filter1d(prof_y_norm, 3)

    x13 = summary["crossings"]["x13p5"]
    y13 = summary["crossings"]["y13p5"]
    x0, x1, _, _ = crop_bounds(x13["left_px"], x13["right_px"], y13["left_px"], y13["right_px"], signal.shape, margin_px=28)
    _, _, y0, y1 = crop_bounds(x13["left_px"], x13["right_px"], y13["left_px"], y13["right_px"], signal.shape, margin_px=28)

    levels = [
        ("90", 0.900, "tab:green", "--"),
        ("86.5", 0.865, "tab:purple", "-."),
        ("50", 0.500, "tab:red", "-"),
        ("13.5", float(np.exp(-2.0)), "tab:blue", "--"),
    ]

    def profile_outer_span(profile_norm: np.ndarray, level: float, search_start: int, search_stop: int) -> dict[str, float]:
        search_start = max(1, search_start)
        search_stop = min(len(profile_norm) - 1, search_stop)
        above_local = np.where(profile_norm[search_start:search_stop] >= level)[0]
        if above_local.size == 0:
            return {"left_px": float("nan"), "right_px": float("nan"), "width_px": float("nan")}

        first = int(search_start + above_local[0])
        last = int(search_start + above_local[-1])

        if first > 0:
            denom = profile_norm[first] - profile_norm[first - 1]
            left = first - 1 + (level - profile_norm[first - 1]) / denom if abs(denom) > 1e-12 else float(first)
        else:
            left = float(first)

        if last < len(profile_norm) - 1:
            denom = profile_norm[last + 1] - profile_norm[last]
            right = last + (level - profile_norm[last]) / denom if abs(denom) > 1e-12 else float(last)
        else:
            right = float(last)

        width_px = float(right - left) if np.isfinite(left) and np.isfinite(right) else float("nan")
        return {"left_px": left, "right_px": right, "width_px": width_px}

    def target_level_for_width(profile_norm: np.ndarray, target_width_px: float, search_start: int, search_stop: int) -> dict[str, float]:
        levels_grid = np.linspace(0.02, 1.45, 2500)
        candidates: list[dict[str, float]] = []
        for level in levels_grid:
            crossing = profile_outer_span(profile_norm, float(level), search_start, search_stop)
            if np.isfinite(crossing["width_px"]):
                candidates.append({"level": float(level), **crossing})
        if not candidates:
            return {"level": float("nan"), "left_px": float("nan"), "right_px": float("nan"), "width_px": float("nan")}
        best = min(candidates, key=lambda item: abs(item["width_px"] - target_width_px))
        return best

    raw_profile_sizes: dict[str, Any] = {
        "note": "Widths are computed from the plotted single-pixel raw profiles after background subtraction and normalization by median(flat). For each threshold, width is the outermost above-threshold span inside the plotted ROI, so internal dips do not split the measured width.",
        "x_profile_y_px": cy,
        "y_profile_x_px": cx,
        "flat_level_intensity": flat_level,
        "levels": {},
        "target_size_equivalent_levels": {},
        "flat_core_gradient_edges_px": {
            "x": summary["gradient_edges_x_px"],
            "y": summary["gradient_edges_y_px"],
        },
        "flat_core_gradient_edges_um": {
            "x": [summary["gradient_edges_x_px"][0] * sx, summary["gradient_edges_x_px"][1] * sx],
            "y": [summary["gradient_edges_y_px"][0] * sy, summary["gradient_edges_y_px"][1] * sy],
        },
    }
    for label, level, _, _ in levels:
        x_cross = profile_outer_span(prof_x_norm, level, x0, x1)
        y_cross = profile_outer_span(prof_y_norm, level, y0, y1)
        raw_profile_sizes["levels"][label] = {
            "fraction_of_flat_level": level,
            "x": {
                **x_cross,
                "left_um": x_cross["left_px"] * sx if np.isfinite(x_cross["left_px"]) else float("nan"),
                "right_um": x_cross["right_px"] * sx if np.isfinite(x_cross["right_px"]) else float("nan"),
                "width_um": x_cross["width_px"] * sx if np.isfinite(x_cross["width_px"]) else float("nan"),
            },
            "y": {
                **y_cross,
                "left_um": y_cross["left_px"] * sy if np.isfinite(y_cross["left_px"]) else float("nan"),
                "right_um": y_cross["right_px"] * sy if np.isfinite(y_cross["right_px"]) else float("nan"),
                "width_um": y_cross["width_px"] * sy if np.isfinite(y_cross["width_px"]) else float("nan"),
            },
        }

    target_specs = {
        "x_target_330um": {
            "axis": "x",
            "target_width_um": 330.0,
            "target_width_px": 330.0 / sx,
            "result": target_level_for_width(prof_x_norm, 330.0 / sx, x0, x1),
        },
        "y_target_120um": {
            "axis": "y",
            "target_width_um": 120.0,
            "target_width_px": 120.0 / sy,
            "result": target_level_for_width(prof_y_norm, 120.0 / sy, y0, y1),
        },
    }
    for key, spec in target_specs.items():
        scale = sx if spec["axis"] == "x" else sy
        result = spec["result"]
        raw_profile_sizes["target_size_equivalent_levels"][key] = {
            "axis": spec["axis"],
            "target_width_um": spec["target_width_um"],
            "target_width_px": spec["target_width_px"],
            "fraction_of_flat_level": result["level"],
            "left_px": result["left_px"],
            "right_px": result["right_px"],
            "width_px": result["width_px"],
            "left_um": result["left_px"] * scale if np.isfinite(result["left_px"]) else float("nan"),
            "right_um": result["right_px"] * scale if np.isfinite(result["right_px"]) else float("nan"),
            "width_um": result["width_px"] * scale if np.isfinite(result["width_px"]) else float("nan"),
        }

    x_um = np.arange(signal.shape[1]) * sx
    y_um = np.arange(signal.shape[0]) * sy

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.4), constrained_layout=True)

    axes[0].plot(x_um, prof_x_norm, color="0.2", lw=0.8, label=profile_label)
    axes[0].plot(x_um, prof_x_light, color="tab:red", lw=0.9, alpha=0.75, label="3-pixel moving average")
    axes[0].set_xlim(x0 * sx, x1 * sx)
    axes[0].set_title(f"X raw profile at y={cy}px")
    axes[0].set_xlabel("x / um")
    axes[0].set_ylabel("I / median(flat)")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(y_um, prof_y_norm, color="0.2", lw=0.8, label=profile_label)
    axes[1].plot(y_um, prof_y_light, color="tab:red", lw=0.9, alpha=0.75, label="3-pixel moving average")
    axes[1].set_xlim(y0 * sy, y1 * sy)
    axes[1].set_title(f"Y raw profile at x={cx}px")
    axes[1].set_xlabel("y / um")
    axes[1].set_ylabel("I / median(flat)")
    axes[1].grid(True, alpha=0.25)

    for ax, axis_name in [(axes[0], "x"), (axes[1], "y")]:
        edge_left, edge_right = raw_profile_sizes["flat_core_gradient_edges_px"][axis_name]
        scale = sx if axis_name == "x" else sy
        ax.axvspan(edge_left * scale, edge_right * scale, color="gold", alpha=0.16, label="gradient-edge flat core")
        for label, level, color, linestyle in levels:
            item = raw_profile_sizes["levels"][label][axis_name]
            if np.isfinite(item["left_px"]) and np.isfinite(item["right_px"]):
                ax.axvline(item["left_px"] * scale, color=color, ls=linestyle, lw=1.1)
                ax.axvline(item["right_px"] * scale, color=color, ls=linestyle, lw=1.1)
            ax.axhline(level, color=color, ls=":", lw=0.8)
        target_key = "x_target_330um" if axis_name == "x" else "y_target_120um"
        target = raw_profile_sizes["target_size_equivalent_levels"][target_key]
        if np.isfinite(target["fraction_of_flat_level"]):
            target_label = f"target {target['target_width_um']:.0f} um @ {target['fraction_of_flat_level']*100:.1f}%"
            ax.axhline(target["fraction_of_flat_level"], color="darkorange", ls="-", lw=1.4, label=target_label)
            ax.axvline(target["left_um"], color="darkorange", ls="-", lw=1.4)
            ax.axvline(target["right_um"], color="darkorange", ls="-", lw=1.4)
            ax.text(
                0.02,
                0.06,
                target_label,
                transform=ax.transAxes,
                color="darkorange",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 2},
            )
        ax.axhline(1.0, color="0.1", ls=":", lw=0.9)
        ax.set_ylim(0.0, max(1.55, float(np.nanpercentile(prof_x_norm if axis_name == "x" else prof_y_norm, 99)) * 1.05))
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "20260605-5 raw center profiles | "
        "threshold widths plus target-size equivalent levels"
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)

    json_out = out.with_name("20260605-5_flattop_raw_profile_sizes.json")
    json_out.write_text(json.dumps(raw_profile_sizes, ensure_ascii=False, indent=2), encoding="utf-8")


def make_expanded_beam_visualization() -> Path:
    matches = sorted(DATA_DIR.glob("*.bgData"))
    if not matches:
        raise FileNotFoundError(f"No .bgData file found in {DATA_DIR}")
    src = matches[0]
    summary_path = next((DATA_DIR / "beam_size_analysis").glob("*_beam_size_summary.json"))
    image, meta = load_frame(src)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    bg = robust_background(image, corner_px=200)
    signal = np.clip(image - bg, 0.0, None)
    sx = float(meta["pixel_scale_x_um"])
    sy = float(meta["pixel_scale_y_um"])
    bbox = summary["beam_region_bbox_px"]
    x0 = max(0, int(bbox[0]) - 120)
    y0 = max(0, int(bbox[1]) - 120)
    x1 = min(signal.shape[1], int(bbox[0] + bbox[2]) + 120)
    y1 = min(signal.shape[0], int(bbox[1] + bbox[3]) + 120)
    crop = signal[y0:y1, x0:x1]
    crop_smooth = gaussian_filter(crop, sigma=2.0)
    rel = crop_smooth / float(np.max(crop_smooth))
    log_rel = np.log10(np.clip(rel, 1e-4, None))
    extent_mm = [x0 * sx / 1000.0, x1 * sx / 1000.0, y1 * sy / 1000.0, y0 * sy / 1000.0]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "expanded_beam_colormap_log.png"

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.8), constrained_layout=True)

    im0 = axes[0].imshow(crop, cmap="turbo", origin="upper", extent=extent_mm, vmin=np.percentile(crop, 1), vmax=np.percentile(crop, 99.6))
    axes[0].set_title("Expanded beam crop")
    axes[0].set_xlabel("x / mm")
    axes[0].set_ylabel("y / mm")
    add_colorbar(fig, axes[0], im0, "intensity")

    im1 = axes[1].imshow(rel, cmap="viridis", origin="upper", extent=extent_mm, vmin=0.0, vmax=1.0)
    axes[1].set_title("Relative linear intensity")
    axes[1].set_xlabel("x / mm")
    add_colorbar(fig, axes[1], im1, "I / peak")

    im2 = axes[2].imshow(log_rel, cmap="magma", origin="upper", extent=extent_mm, vmin=-3.0, vmax=0.0)
    axes[2].set_title("Log intensity")
    axes[2].set_xlabel("x / mm")
    add_colorbar(fig, axes[2], im2, "log10(I / peak)")

    cx, cy = summary["moments"]["center_px_global"]
    for ax in axes:
        ax.plot(cx * sx / 1000.0, cy * sy / 1000.0, "w+", ms=11, mew=1.6)
        ax.set_aspect("equal")

    fig.suptitle(
        "Expanded beam visualization | "
        f"1/e^2 profile = {summary['profile_widths']['x_1e2']['width_um']/1000:.2f} x "
        f"{summary['profile_widths']['y_1e2']['width_um']/1000:.2f} mm"
    )
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def main() -> None:
    flat = make_flattop_visualization()
    beam = make_expanded_beam_visualization()
    print(f"Flat-top visualization: {flat}")
    print(f"Flat-top raw profiles: {OUT_DIR / '20260605-5_flattop_raw_profiles.png'}")
    print(f"Flat-top raw profile sizes: {OUT_DIR / '20260605-5_flattop_raw_profile_sizes.json'}")
    print(f"Expanded-beam visualization: {beam}")


if __name__ == "__main__":
    main()
