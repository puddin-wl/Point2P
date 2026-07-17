"""Analyze rectangular flat-top Spiricon bmData using energy and profile metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
from scipy.ndimage import uniform_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_rect_flattop_size import analyze, load_spiricon_frame, robust_background


def load_beamgage_percent_energy_setup(path: Path, shape: tuple[int, int]) -> dict[str, Any]:
    ny, nx = shape
    with h5py.File(path, "r") as h5:
        root = h5["BG_SETUP/RESULTS_MANAGER/RESULT_ENGINE"]
        aperture = root["MANUAL_APERTURE"]
        settings = root["PROGRAMMABLE_SETTINGS_MANAGER"]

        def scalar(group: h5py.Group, name: str) -> Any:
            value = group[name][()]
            if getattr(value, "shape", None) == (1,):
                value = value[0]
            if isinstance(value, bytes):
                return value.decode("utf-8", errors="replace")
            return value.item() if hasattr(value, "item") else value

        cx = float(scalar(aperture, "CenterXPercent")) * nx
        cy = float(scalar(aperture, "CenterYPercent")) * ny
        width = float(scalar(aperture, "WidthPercent")) * nx
        height = float(scalar(aperture, "HeightPercent")) * ny
        x0 = max(0, int(np.floor(cx - width / 2.0)))
        x1 = min(nx, int(np.ceil(cx + width / 2.0)))
        y0 = max(0, int(np.floor(cy - height / 2.0)))
        y1 = min(ny, int(np.ceil(cy + height / 2.0)))
        return {
            "beam_width_type": scalar(settings["BEAM_WIDTH_BASIS"], "BEAM_WIDTH_TYPE"),
            "percent_total": float(scalar(settings["PERCENT_TOTAL_CLIP"], "PERCENT_TOTAL")),
            "manual_aperture_enabled": bool(scalar(aperture, "IsEnabled")),
            "manual_aperture_shape": scalar(aperture, "Shape"),
            "manual_aperture_center_px": [cx, cy],
            "manual_aperture_size_px": [width, height],
            "manual_aperture_bounds_px": [x0, y0, x1 - x0, y1 - y0],
        }


def beamgage_percent_energy_bbox(signal: np.ndarray, setup: dict[str, Any], sx: float, sy: float) -> dict[str, Any]:
    x0, y0, width, height = setup["manual_aperture_bounds_px"]
    x1, y1 = x0 + width, y0 + height
    roi = np.clip(signal[y0:y1, x0:x1], 0.0, None)
    fraction = float(setup["percent_total"])
    sorted_values = np.sort(roi.ravel())[::-1]
    cumulative = np.cumsum(sorted_values)
    threshold_index = int(np.searchsorted(cumulative, fraction * cumulative[-1]))
    threshold = float(sorted_values[min(threshold_index, sorted_values.size - 1)])
    selected = roi >= threshold
    ys, xs = np.where(selected)
    gx0, gx1 = int(xs.min() + x0), int(xs.max() + x0)
    gy0, gy1 = int(ys.min() + y0), int(ys.max() + y0)
    return {
        "definition": "Sort all background-subtracted pixels inside the saved manual aperture by intensity, retain pixels until their sum reaches PercentEnergy, then report the coordinate span of the retained-pixel bounding box.",
        "percent_energy": fraction,
        "intensity_threshold": threshold,
        "selected_pixel_count": int(selected.sum()),
        "selected_energy_fraction": float(roi[selected].sum() / roi.sum()),
        "selected_bbox_global_px_inclusive": [gx0, gy0, gx1, gy1],
        "selected_bbox_count_px": [gx1 - gx0 + 1, gy1 - gy0 + 1],
        "selected_bbox_span_px": [gx1 - gx0, gy1 - gy0],
        "pixel_coverage_width_x_um": float((gx1 - gx0 + 1) * sx),
        "pixel_coverage_width_y_um": float((gy1 - gy0 + 1) * sy),
        "width_x_um": float((gx1 - gx0) * sx),
        "width_y_um": float((gy1 - gy0) * sy),
    }


def outer_span(profile: np.ndarray, level: float, start: int, stop: int) -> dict[str, float]:
    start = max(1, int(start))
    stop = min(len(profile) - 1, int(stop))
    indices = np.where(profile[start:stop] >= level)[0]
    if indices.size == 0:
        return {"left_px": float("nan"), "right_px": float("nan"), "width_px": float("nan")}
    first = start + int(indices[0])
    last = start + int(indices[-1])
    left = float(first)
    right = float(last)
    if first > 0:
        dy = profile[first] - profile[first - 1]
        if abs(dy) > 1e-12:
            left = first - 1 + (level - profile[first - 1]) / dy
    if last < len(profile) - 1:
        dy = profile[last + 1] - profile[last]
        if abs(dy) > 1e-12:
            right = last + (level - profile[last]) / dy
    return {"left_px": left, "right_px": right, "width_px": float(right - left)}


def minimum_energy_interval(profile: np.ndarray, fraction: float) -> dict[str, float]:
    values = np.clip(np.asarray(profile, dtype=np.float64), 0.0, None)
    total = float(values.sum())
    if total <= 0:
        return {"start_px": float("nan"), "stop_px": float("nan"), "count_px": float("nan"), "fraction": fraction}
    target = fraction * total
    cumulative = np.concatenate(([0.0], np.cumsum(values)))
    best: tuple[int, int, float] | None = None
    stop = 0
    for start in range(values.size):
        stop = max(stop, start)
        while stop < values.size and cumulative[stop + 1] - cumulative[start] < target:
            stop += 1
        if stop >= values.size:
            break
        energy = cumulative[stop + 1] - cumulative[start]
        width = stop - start + 1
        if best is None or width < best[1] - best[0] + 1 or (width == best[1] - best[0] + 1 and energy < best[2]):
            best = (start, stop, energy)
    if best is None:
        return {"start_px": 0.0, "stop_px": float(values.size - 1), "count_px": float(values.size), "fraction": 1.0}
    start, stop, energy = best
    return {
        "start_px": float(start),
        "stop_px": float(stop),
        "count_px": float(stop - start + 1),
        "fraction": float(energy / total),
    }


def half_stats(values: np.ndarray, center: int, axis: int) -> dict[str, float]:
    if axis == 0:
        first = values[:center, :]
        second = values[center:, :]
        first_name, second_name = "top", "bottom"
    else:
        first = values[:, :center]
        second = values[:, center:]
        first_name, second_name = "left", "right"
    mean_first = float(np.mean(first))
    mean_second = float(np.mean(second))
    denom = 0.5 * (mean_first + mean_second)
    return {
        f"{first_name}_mean": mean_first,
        f"{second_name}_mean": mean_second,
        f"{first_name}_over_{second_name}": float(mean_first / mean_second) if mean_second > 0 else float("nan"),
        "difference_percent_of_pair_mean": float((mean_first - mean_second) / denom * 100.0) if denom > 0 else float("nan"),
    }


def run(
    path: Path,
    outdir: Path,
    phase_label: str | None = None,
    reported_width_x_um: float | None = None,
    reported_width_y_um: float | None = None,
    main_spot_only: bool = False,
) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    baseline = analyze(path, outdir, threshold_frac=0.10, sigma_px=1.5, corner_px=30)
    image, meta = load_spiricon_frame(path)
    bg, bg_sigma = robust_background(image, corner_px=30)
    signal = np.clip(image - bg, 0.0, None)
    sx = float(meta["pixel_scale_x_um"])
    sy = float(meta["pixel_scale_y_um"])
    saved_beamgage_setup = load_beamgage_percent_energy_setup(path, signal.shape)
    cx, cy = baseline["flat_center_px"]
    cxi, cyi = int(round(cx)), int(round(cy))

    x13 = baseline["crossings"]["x13p5"]
    y13 = baseline["crossings"]["y13p5"]
    margin = 5 if main_spot_only else 28
    x0 = max(0, int(np.floor(x13["left_px"])) - margin)
    x1 = min(signal.shape[1], int(np.ceil(x13["right_px"])) + margin + 1)
    y0 = max(0, int(np.floor(y13["left_px"])) - margin)
    y1 = min(signal.shape[0], int(np.ceil(y13["right_px"])) + margin + 1)
    roi = signal[y0:y1, x0:x1]
    beamgage_setup = dict(saved_beamgage_setup)
    if main_spot_only:
        beamgage_setup["manual_aperture_enabled"] = True
        beamgage_setup["manual_aperture_shape"] = "Main rectangular spot only"
        beamgage_setup["manual_aperture_center_px"] = [0.5 * (x0 + x1), 0.5 * (y0 + y1)]
        beamgage_setup["manual_aperture_size_px"] = [x1 - x0, y1 - y0]
        beamgage_setup["manual_aperture_bounds_px"] = [x0, y0, x1 - x0, y1 - y0]
    beamgage_result = beamgage_percent_energy_bbox(signal, beamgage_setup, sx, sy)

    projection_x = roi.sum(axis=0)
    projection_y = roi.sum(axis=1)
    energy: dict[str, Any] = {}
    for label, fraction in (("50", 0.50), ("86.5", 0.865), ("90", 0.90)):
        ex = minimum_energy_interval(projection_x, fraction)
        ey = minimum_energy_interval(projection_y, fraction)
        ex["width_um"] = ex["count_px"] * sx
        ey["width_um"] = ey["count_px"] * sy
        ex["global_start_px"] = ex["start_px"] + x0
        ex["global_stop_px"] = ex["stop_px"] + x0
        ey["global_start_px"] = ey["start_px"] + y0
        ey["global_stop_px"] = ey["stop_px"] + y0
        energy[label] = {"x": ex, "y": ey}

    band = 2
    profile_x = signal[max(0, cyi - band) : min(signal.shape[0], cyi + band + 1), :].mean(axis=0)
    profile_y = signal[:, max(0, cxi - band) : min(signal.shape[1], cxi + band + 1)].mean(axis=1)
    profile_x = uniform_filter1d(profile_x, 3)
    profile_y = uniform_filter1d(profile_y, 3)

    lx, rx = baseline["gradient_edges_x_px"]
    ty, by = baseline["gradient_edges_y_px"]
    x90 = baseline["crossings"]["x90"]
    y90 = baseline["crossings"]["y90"]
    core_x0 = max(0, int(np.ceil(x90["left_px"])))
    core_x1 = min(signal.shape[1], int(np.floor(x90["right_px"])) + 1)
    core_y0 = max(0, int(np.ceil(y90["left_px"])))
    core_y1 = min(signal.shape[0], int(np.floor(y90["right_px"])) + 1)
    core = signal[core_y0:core_y1, core_x0:core_x1]
    flat_level = float(np.median(core))
    peak = float(np.max(core))
    profile_metrics: dict[str, Any] = {}
    for label, fraction in (("90", 0.90), ("86.5", 0.865), ("50", 0.50), ("13.5", float(np.exp(-2.0)))):
        xf = outer_span(profile_x / flat_level, fraction, x0, x1)
        yf = outer_span(profile_y / flat_level, fraction, y0, y1)
        xp = outer_span(profile_x / peak, fraction, x0, x1)
        yp = outer_span(profile_y / peak, fraction, y0, y1)
        for item, scale in ((xf, sx), (yf, sy), (xp, sx), (yp, sy)):
            item["width_um"] = item["width_px"] * scale
        profile_metrics[label] = {
            "fraction": fraction,
            "relative_to_flat": {"x": xf, "y": yf},
            "relative_to_peak": {"x": xp, "y": yp},
        }

    core_cx = int(np.clip(round(cx) - core_x0, 1, core.shape[1] - 1))
    core_cy = int(np.clip(round(cy) - core_y0, 1, core.shape[0] - 1))
    asymmetry = {
        "top_bottom": half_stats(core, core_cy, axis=0),
        "left_right": half_stats(core, core_cx, axis=1),
        "size90_box_px": [core_x0, core_y0, core_x1 - core_x0, core_y1 - core_y0],
        "core_mean": float(np.mean(core)),
        "core_std": float(np.std(core)),
        "core_rms_percent": float(np.std(core) / np.mean(core) * 100.0),
    }

    result = {
        "source_file": str(path.resolve()),
        "phase_label": phase_label,
        "onsite_reported_size_um": {
            "x": reported_width_x_um,
            "y": reported_width_y_um,
            "note": "User-confirmed onsite BeamGage readout; kept separate from raw selected-pixel bbox conventions.",
        },
        "metadata": meta,
        "background": {"median": bg, "robust_sigma": bg_sigma},
        "roi_px": [x0, y0, x1 - x0, y1 - y0],
        "center_px": [cx, cy],
        "baseline_profile_widths": baseline["widths"],
        "beamgage_saved_setup": saved_beamgage_setup,
        "percent_energy_analysis_aperture": beamgage_setup,
        "main_spot_only": bool(main_spot_only),
        "beamgage_percent_energy_reproduction": beamgage_result,
        "minimum_interval_percent_energy": energy,
        "profile_widths": profile_metrics,
        "flat_core_asymmetry": asymmetry,
    }

    x_axis = np.arange(signal.shape[1]) * sx
    y_axis = np.arange(signal.shape[0]) * sy
    extent = [x0 * sx, x1 * sx, y1 * sy, y0 * sy]
    norm_roi = roi / flat_level
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    im = axes[0, 0].imshow(norm_roi, cmap="turbo", origin="upper", extent=extent, vmin=0, vmax=1.35)
    axes[0, 0].set_title("Measured morphology (background-subtracted)")
    axes[0, 0].set_xlabel("x / um")
    axes[0, 0].set_ylabel("y / um")
    fig.colorbar(im, ax=axes[0, 0], label="relative intensity")

    axes[0, 1].imshow(norm_roi, cmap="turbo", origin="upper", extent=extent, vmin=0, vmax=1.35)
    axes[0, 1].set_title("BeamGage PercentEnergy size (86.5%)")
    axes[0, 1].set_xlabel("x / um")
    axes[0, 1].set_ylabel("y / um")
    ax0, ay0, aw, ah = beamgage_setup["manual_aperture_bounds_px"]
    axes[0, 1].add_patch(
        plt.Rectangle(
            (ax0 * sx, ay0 * sy),
            aw * sx,
            ah * sy,
            fill=False,
            ec="white",
            lw=1.3,
            ls="--",
            label="main-spot analysis aperture" if main_spot_only else "saved manual aperture",
        )
    )
    bx0, by0, bx1, by1 = beamgage_result["selected_bbox_global_px_inclusive"]
    box_x0, box_x1 = bx0 * sx, bx1 * sx
    box_y0, box_y1 = by0 * sy, by1 * sy
    display_width_x = beamgage_result["width_x_um"] if reported_width_x_um is None else reported_width_x_um
    display_width_y = beamgage_result["width_y_um"] if reported_width_y_um is None else reported_width_y_um
    box_cx = 0.5 * (box_x0 + box_x1)
    box_cy = 0.5 * (box_y0 + box_y1)
    display_x0 = box_cx - 0.5 * display_width_x
    display_x1 = box_cx + 0.5 * display_width_x
    display_y0 = box_cy - 0.5 * display_width_y
    display_y1 = box_cy + 0.5 * display_width_y
    axes[0, 1].add_patch(
        plt.Rectangle(
            (display_x0, display_y0),
            display_width_x,
            display_width_y,
            fill=False,
            ec="magenta",
            lw=1.5,
            label="onsite BeamGage 86.5% size",
        )
    )
    arrow_y = display_y0 - 14.0
    arrow_x = display_x1 + 14.0
    axes[0, 1].annotate("", xy=(display_x1, arrow_y), xytext=(display_x0, arrow_y), arrowprops={"arrowstyle": "<->", "color": "white", "lw": 1.6})
    axes[0, 1].text(
        0.5 * (display_x0 + display_x1),
        arrow_y - 5.0,
        f"X = {display_width_x:.2f} um",
        ha="center",
        va="bottom",
        color="white",
        fontsize=10,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
    )
    axes[0, 1].annotate("", xy=(arrow_x, display_y1), xytext=(arrow_x, display_y0), arrowprops={"arrowstyle": "<->", "color": "white", "lw": 1.6})
    axes[0, 1].text(
        arrow_x + 5.0,
        0.5 * (display_y0 + display_y1),
        f"Y = {display_width_y:.2f} um",
        ha="left",
        va="center",
        rotation=90,
        color="white",
        fontsize=10,
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 2},
    )
    axes[0, 1].legend(loc="lower right", fontsize=8)

    axes[1, 0].plot(x_axis, profile_x / flat_level, color="black", lw=1)
    axes[1, 0].axhline(1.0, color="gray", ls=":")
    axes[1, 0].set_xlim(x0 * sx, x1 * sx)
    axes[1, 0].set_title("X center profile (5-pixel band mean)")
    axes[1, 0].set_xlabel("x / um")
    axes[1, 0].set_ylabel("I / median(flat core)")

    axes[1, 1].plot(y_axis, profile_y / flat_level, color="black", lw=1)
    axes[1, 1].axhline(1.0, color="gray", ls=":")
    axes[1, 1].set_xlim(y0 * sy, y1 * sy)
    axes[1, 1].set_title("Y center profile (5-pixel band mean)")
    axes[1, 1].set_xlabel("y / um")
    axes[1, 1].set_ylabel("I / median(flat core)")

    axes[0, 0].set_aspect("equal")
    axes[0, 1].set_aspect("equal")

    title = path.name if not phase_label else f"{path.name}\nPhase: {phase_label}"
    fig.suptitle(title)
    figure_path = outdir / f"{path.stem}_percent_energy_diagnostics.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)
    json_path = outdir / f"{path.stem}_percent_energy_summary.json"
    result["outputs"] = {"summary_json": str(json_path), "diagnostic_png": str(figure_path)}
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--phase-label", default=None, help="Optional phase/run label shown in the diagnostic title and JSON output.")
    parser.add_argument("--reported-width-x", type=float, default=None, help="Optional onsite-reported BeamGage X size in um.")
    parser.add_argument("--reported-width-y", type=float, default=None, help="Optional onsite-reported BeamGage Y size in um.")
    parser.add_argument("--main-spot-only", action="store_true", help="Exclude zero-order and other disconnected spots from the analysis aperture.")
    args = parser.parse_args()
    result = run(
        args.input,
        args.outdir,
        phase_label=args.phase_label,
        reported_width_x_um=args.reported_width_x,
        reported_width_y_um=args.reported_width_y,
        main_spot_only=args.main_spot_only,
    )
    vendor = result["beamgage_percent_energy_reproduction"]
    e = result["minimum_interval_percent_energy"]["86.5"]
    a = result["flat_core_asymmetry"]["top_bottom"]
    print(f"Raw selected-mask center span: X={vendor['width_x_um']:.1f} um, Y={vendor['width_y_um']:.1f} um")
    if args.reported_width_x is not None and args.reported_width_y is not None:
        print(f"Onsite reported: X={args.reported_width_x:.1f} um, Y={args.reported_width_y:.1f} um")
    print(f"86.5% energy: X={e['x']['width_um']:.1f} um, Y={e['y']['width_um']:.1f} um")
    print(f"Top/bottom mean ratio: {a['top_over_bottom']:.4f}")
    print(f"Top-bottom difference: {a['difference_percent_of_pair_mean']:+.2f}%")


if __name__ == "__main__":
    main()
