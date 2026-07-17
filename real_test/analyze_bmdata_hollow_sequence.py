"""Quantify center depression across every frame of a BeamGage bmData sequence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analyze_rect_flattop_size import load_spiricon_frame, robust_background


def stats(values: list[float]) -> dict[str, float]:
    data = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(data)),
        "std": float(np.std(data)),
        "min": float(np.min(data)),
        "max": float(np.max(data)),
    }


def run(path: Path, percent_summary: Path, outdir: Path) -> dict[str, Any]:
    setup = json.loads(percent_summary.read_text(encoding="utf-8"))
    cx, cy = (int(round(v)) for v in setup["center_px"])
    x50 = setup["profile_widths"]["50"]["relative_to_flat"]["x"]
    y50 = setup["profile_widths"]["50"]["relative_to_flat"]["y"]
    x0 = int(np.ceil(x50["left_px"]))
    x1 = int(np.floor(x50["right_px"])) + 1
    y0 = int(np.ceil(y50["left_px"]))
    y1 = int(np.floor(y50["right_px"])) + 1

    with h5py.File(path, "r") as h5:
        frame_keys = sorted(h5["BG_DATA"].keys(), key=int)

    records: list[dict[str, Any]] = []
    signals: list[np.ndarray] = []
    for key in frame_keys:
        image, meta = load_spiricon_frame(path, frame_key=key)
        background, background_sigma = robust_background(image, corner_px=30)
        signal = np.clip(image - background, 0.0, None)
        signals.append(signal)
        core = signal[y0:y1, x0:x1]
        core_mean = float(np.mean(core))
        width = core.shape[1]
        third = max(1, width // 3)
        middle = core[:, third : width - third]
        sides = np.concatenate((core[:, :third].ravel(), core[:, width - third :].ravel()))
        records.append(
            {
                "frame_key": key,
                "timestamp": meta["timestamp"],
                "background": background,
                "background_robust_sigma": background_sigma,
                "center_5x5_over_size50_box_mean": float(np.mean(signal[cy - 2 : cy + 3, cx - 2 : cx + 3]) / core_mean),
                "center_15x15_over_size50_box_mean": float(np.mean(signal[cy - 7 : cy + 8, cx - 7 : cx + 8]) / core_mean),
                "middle_third_over_side_thirds": float(np.mean(middle) / np.mean(sides)),
                "size50_box_rms_percent": float(np.std(core) / core_mean * 100.0),
            }
        )

    average = np.mean(np.stack(signals, axis=0), axis=0)
    sx = float(setup["metadata"]["pixel_scale_x_um"])
    sy = float(setup["metadata"]["pixel_scale_y_um"])
    center5_values = [r["center_5x5_over_size50_box_mean"] for r in records]
    center15_values = [r["center_15x15_over_size50_box_mean"] for r in records]
    middle_values = [r["middle_third_over_side_thirds"] for r in records]
    rms_values = [r["size50_box_rms_percent"] for r in records]
    result = {
        "source_file": str(path.resolve()),
        "frame_count": len(records),
        "analysis_center_px": [cx, cy],
        "size50_analysis_box_px": [x0, y0, x1 - x0, y1 - y0],
        "definitions": {
            "center_5x5": "Mean of the 5x5 pixels centered at the flat-top center divided by the mean inside the frame-1 size50 box.",
            "middle_third": "Mean of the middle x-third of the size50 box divided by the combined left/right thirds.",
        },
        "aggregate": {
            "center_5x5_over_size50_box_mean": stats(center5_values),
            "center_15x15_over_size50_box_mean": stats(center15_values),
            "middle_third_over_side_thirds": stats(middle_values),
            "size50_box_rms_percent": stats(rms_values),
        },
        "frames": records,
    }

    outdir.mkdir(parents=True, exist_ok=True)
    json_path = outdir / f"{path.stem}_hollow_sequence_summary.json"
    figure_path = outdir / f"{path.stem}_hollow_sequence_diagnostics.png"
    result["outputs"] = {"summary_json": str(json_path), "diagnostic_png": str(figure_path)}
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    margin_x, margin_y = 8, 7
    zx0, zx1 = max(0, x0 - margin_x), min(average.shape[1], x1 + margin_x)
    zy0, zy1 = max(0, y0 - margin_y), min(average.shape[0], y1 + margin_y)
    core_mean = float(np.mean(average[y0:y1, x0:x1]))
    normalized = average / core_mean
    extent = [zx0 * sx, zx1 * sx, zy1 * sy, zy0 * sy]
    x_axis = np.arange(average.shape[1]) * sx
    y_axis = np.arange(average.shape[0]) * sy
    band = 2
    profile_x = normalized[cy - band : cy + band + 1, :].mean(axis=0)
    profile_y = normalized[:, cx - band : cx + band + 1].mean(axis=1)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    im = axes[0, 0].imshow(
        normalized[zy0:zy1, zx0:zx1],
        origin="upper",
        cmap="turbo",
        extent=extent,
        vmin=0.45,
        vmax=1.45,
        aspect="equal",
    )
    axes[0, 0].set_title(f"Average of {len(records)} frames")
    axes[0, 0].set_xlabel("x / um")
    axes[0, 0].set_ylabel("y / um")
    fig.colorbar(im, ax=axes[0, 0], label="I / mean(size50 box)")

    indices = np.arange(1, len(records) + 1)
    axes[0, 1].plot(indices, center5_values, "o-", label="center 5x5 / box mean")
    axes[0, 1].plot(indices, middle_values, "s-", label="middle third / side thirds")
    axes[0, 1].axhline(1.0, color="gray", ls=":")
    axes[0, 1].set_title("Center-depression stability")
    axes[0, 1].set_xlabel("frame")
    axes[0, 1].set_ylabel("ratio")
    axes[0, 1].set_xticks(indices)
    axes[0, 1].grid(True, alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].plot(x_axis, profile_x, color="black")
    axes[1, 0].axhline(1.0, color="gray", ls=":")
    axes[1, 0].set_xlim(zx0 * sx, zx1 * sx)
    axes[1, 0].set_title("Average X center profile")
    axes[1, 0].set_xlabel("x / um")
    axes[1, 0].set_ylabel("I / mean(size50 box)")
    axes[1, 0].grid(True, alpha=0.25)

    axes[1, 1].plot(y_axis, profile_y, color="black")
    axes[1, 1].axhline(1.0, color="gray", ls=":")
    axes[1, 1].set_xlim(zy0 * sy, zy1 * sy)
    axes[1, 1].set_title("Average Y center profile")
    axes[1, 1].set_xlabel("y / um")
    axes[1, 1].set_ylabel("I / mean(size50 box)")
    axes[1, 1].grid(True, alpha=0.25)
    fig.suptitle(f"Measured center depression across BeamGage sequence: {path.name}")
    fig.savefig(figure_path, dpi=190)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--percent-summary", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    result = run(args.input, args.percent_summary, args.outdir)
    print(json.dumps(result["aggregate"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
