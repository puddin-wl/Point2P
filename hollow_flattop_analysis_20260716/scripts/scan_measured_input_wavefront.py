"""Scan low-order wavefront errors on the measured 2026-07-16 input amplitude."""

from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import matplotlib
import numpy as np
from scipy.ndimage import map_coordinates

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulate_v2_hollow_scan import (
    ANALYSIS_ROOT,
    DEFAULT_BMDATA,
    DEFAULT_CASE_DIR,
    DEFAULT_HOLLOW_SUMMARY,
    DEFAULT_PERCENT_SUMMARY,
    REAL_TEST_ROOT,
    Simulator,
    evaluate_image,
    load_experimental_target,
    normalize_power,
)
from analyze_rect_flattop_size import load_spiricon_frame, robust_background
from src.propagation import forward_fft, intensity


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "04_measured_input_low_order_wavefront"
INPUT_BGDATA = REAL_TEST_ROOT / "20260716" / "G-光斑-1.bgData"
INPUT_SUMMARY = (
    REAL_TEST_ROOT
    / "20260716"
    / "analysis_G_spot_1"
    / "G-光斑-1_beam_size_summary.json"
)


def measured_amplitude_on_doe(simulator: Simulator) -> np.ndarray:
    summary = json.loads(INPUT_SUMMARY.read_text(encoding="utf-8"))
    image, _ = load_spiricon_frame(INPUT_BGDATA)
    background, _ = robust_background(image, corner_px=200)
    signal = np.clip(image - background, 0.0, None).astype(np.float32)
    x0, y0, width, height = (int(value) for value in summary["beam_region_bbox_px"])
    selected = np.zeros_like(signal)
    selected[y0 : y0 + height, x0 : x0 + width] = signal[
        y0 : y0 + height, x0 : x0 + width
    ]
    center_x_px, center_y_px = (
        float(value) for value in summary["moments"]["center_px_global"]
    )
    sx_um = float(summary["pixel_scale_x_um"])
    sy_um = float(summary["pixel_scale_y_um"])
    source_x_px = center_x_px + simulator.X * np.float32(1e6 / sx_um)
    source_y_px = center_y_px + simulator.Y * np.float32(1e6 / sy_um)
    sampled_intensity = map_coordinates(
        selected,
        [source_y_px, source_x_px],
        order=1,
        mode="constant",
        cval=0.0,
        prefilter=False,
    ).astype(np.float32)
    return normalize_power(np.sqrt(sampled_intensity).astype(np.float32))


def reconstruct(
    simulator: Simulator,
    measured_amplitude: np.ndarray,
    defocus: float,
    astigmatism: float,
    spherical: float,
) -> np.ndarray:
    phase_waves = (
        defocus * simulator.normalized_r2
        + astigmatism
        * (
            simulator.normalized_x * simulator.normalized_x
            - simulator.normalized_y * simulator.normalized_y
        )
        + spherical * simulator.normalized_r2 * simulator.normalized_r2
    )
    aberration = np.exp(1j * np.float32(2.0 * np.pi) * phase_waves).astype(
        np.complex64
    )
    field = measured_amplitude * simulator.phase_factor * aberration
    return intensity(forward_fft(field.astype(np.complex64), np), np).astype(
        np.float32
    )


def plot_overview(
    output_path: Path,
    experiment: dict,
    simulator: Simulator,
    measured_amplitude: np.ndarray,
    rows: list[dict],
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), constrained_layout=True)
    panels = [
        (
            "Experiment",
            experiment["normalized"],
            experiment["x_um"],
            experiment["y_um"],
            experiment["metrics"],
        )
    ]
    for row in rows[:5]:
        params = row["params"]
        image = reconstruct(
            simulator,
            measured_amplitude,
            params["defocus_waves"],
            params["astigmatism_waves"],
            params["spherical_waves"],
        )
        normalized, metrics, _ = evaluate_image(image, simulator, experiment)
        panels.append(
            (
                f"{row['name']}\nscore={row['score']:.3f}",
                normalized,
                simulator.base.x_um,
                simulator.base.y_um,
                metrics,
            )
        )
    for ax, (title, image, x_um, y_um, metrics) in zip(axes.ravel(), panels):
        xmask = np.abs(x_um) <= 220
        ymask = np.abs(y_um) <= 100
        roi = image[np.ix_(ymask, xmask)]
        xa = x_um[xmask]
        ya = y_um[ymask]
        im = ax.imshow(
            roi,
            origin="upper",
            extent=[xa[0], xa[-1], ya[-1], ya[0]],
            cmap="turbo",
            vmin=0.35,
            vmax=1.50,
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
    fig.suptitle("Measured input amplitude + low-order wavefront scan, fixed V2 phase")
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        DEFAULT_BMDATA, DEFAULT_PERCENT_SUMMARY, DEFAULT_HOLLOW_SUMMARY
    )
    simulator = Simulator(DEFAULT_CASE_DIR)
    measured_amplitude = measured_amplitude_on_doe(simulator)
    values = np.linspace(-0.6, 0.6, 7)
    rows: list[dict] = []
    started = time.perf_counter()
    total = len(values) ** 3
    index = 0
    for defocus in values:
        for astigmatism in values:
            for spherical in values:
                index += 1
                params = {
                    "defocus_waves": float(defocus),
                    "astigmatism_waves": float(astigmatism),
                    "spherical_waves": float(spherical),
                }
                image = reconstruct(
                    simulator,
                    measured_amplitude,
                    float(defocus),
                    float(astigmatism),
                    float(spherical),
                )
                _, metrics, score_details = evaluate_image(
                    image, simulator, experiment
                )
                name = (
                    f"D{defocus:+.2f}_A{astigmatism:+.2f}_S{spherical:+.2f}waves"
                )
                rows.append(
                    {
                        "family": "measured_input_low_order_wavefront",
                        "name": name,
                        "params": params,
                        "score": score_details["score"],
                        "metrics": metrics,
                        "score_details": score_details,
                    }
                )
                if index == 1 or index % 40 == 0 or index == total:
                    print(
                        f"[{index:3d}/{total}] {name} "
                        f"score={score_details['score']:.3f}",
                        flush=True,
                    )
    elapsed = time.perf_counter() - started
    rows.sort(key=lambda row: row["score"])
    (OUTPUT_DIR / "wavefront_scan_results.json").write_text(
        json.dumps(
            {
                "phase_source": str(DEFAULT_CASE_DIR / "phase_refined.npy"),
                "input_source": str(INPUT_BGDATA),
                "candidate_count": len(rows),
                "elapsed_seconds": elapsed,
                "wave_definition": (
                    "phase/2pi = D*rho^2 + A*(x_normalized^2-y_normalized^2) "
                    "+ S*rho^4; normalized radius is 3.25 mm"
                ),
                "candidates": rows,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    fields = [
        "rank",
        "score",
        "name",
        "defocus_waves",
        "astigmatism_waves",
        "spherical_waves",
        "center_ratio",
        "middle_over_sides",
        "core_rms_fraction",
        "profile_rmse_x",
        "profile_rmse_y",
    ]
    with (OUTPUT_DIR / "wavefront_scan_results.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            writer.writerow(
                {
                    "rank": rank,
                    "score": row["score"],
                    "name": row["name"],
                    **row["params"],
                    "center_ratio": row["metrics"][
                        "center_window_over_core_mean"
                    ],
                    "middle_over_sides": row["metrics"][
                        "middle_third_over_side_thirds"
                    ],
                    "core_rms_fraction": row["metrics"]["core_rms_fraction"],
                    "profile_rmse_x": row["score_details"]["profile_rmse_x"],
                    "profile_rmse_y": row["score_details"]["profile_rmse_y"],
                }
            )
    best = rows[0]
    best_params = best["params"]
    best_image = reconstruct(
        simulator,
        measured_amplitude,
        best_params["defocus_waves"],
        best_params["astigmatism_waves"],
        best_params["spherical_waves"],
    )
    np.save(OUTPUT_DIR / "best_wavefront_focal_intensity.npy", best_image)
    (OUTPUT_DIR / "best_wavefront.json").write_text(
        json.dumps(best, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_overview(
        OUTPUT_DIR / "measured_input_wavefront_overview.png",
        experiment,
        simulator,
        measured_amplitude,
        rows,
    )
    print(
        json.dumps(
            {
                "candidate_count": len(rows),
                "elapsed_seconds": elapsed,
                "best": {
                    "name": best["name"],
                    "params": best["params"],
                    "score": best["score"],
                    "metrics": {
                        key: value
                        for key, value in best["metrics"].items()
                        if not key.startswith("unit_profile")
                    },
                },
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
