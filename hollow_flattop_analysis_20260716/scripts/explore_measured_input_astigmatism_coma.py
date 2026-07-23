"""Exploratory Zernike astigmatism/coma cases for the measured Gaussian input.

This is intentionally separate from the main spherical-aberration conclusion.
It uses the fixed V2 phase, measured 2026-07-16 Gaussian amplitude, and the
15 mm clear pupil. Coefficients are Noll-normalized RMS waves on that pupil.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scan_measured_input_wavefront import measured_amplitude_on_doe
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


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "exploratory_14_astigmatism_coma"
OUTPUT_FIGURE = OUTPUT_DIR / "ASTIGMATISM_COMA_EXPLORATORY_MONTAGE.png"
COEFFICIENT_RMS_WAVES = 0.10


def compact_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in metrics.items()
        if not key.startswith("unit_profile")
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        DEFAULT_BMDATA,
        DEFAULT_PERCENT_SUMMARY,
        DEFAULT_HOLLOW_SUMMARY,
    )
    simulator = Simulator(DEFAULT_CASE_DIR)
    measured_amplitude = measured_amplitude_on_doe(simulator)

    pupil_radius_m = float(simulator.clear_aperture_m / 2.0)
    xn = simulator.X / np.float32(pupil_radius_m)
    yn = simulator.Y / np.float32(pupil_radius_m)
    rho2 = xn * xn + yn * yn

    # Noll-normalized real Zernike modes on the unit disk.
    modes = {
        "astig_0": np.float32(math.sqrt(6.0)) * (xn * xn - yn * yn),
        "astig_45": np.float32(2.0 * math.sqrt(6.0)) * xn * yn,
        "coma_x": np.float32(math.sqrt(8.0)) * (3.0 * rho2 - 2.0) * xn,
        "coma_y": np.float32(math.sqrt(8.0)) * (3.0 * rho2 - 2.0) * yn,
    }
    c = COEFFICIENT_RMS_WAVES
    cases: list[tuple[str, dict[str, float]]] = [
        ("No added aberration", {}),
        (f"Astig 0 deg  {c:+.2f}", {"astig_0": +c}),
        (f"Astig 0 deg  {-c:+.2f}", {"astig_0": -c}),
        (f"Astig 45 deg {c:+.2f}", {"astig_45": +c}),
        (f"Coma X       {c:+.2f}", {"coma_x": +c}),
        (f"Coma X       {-c:+.2f}", {"coma_x": -c}),
        (f"Coma Y       {c:+.2f}", {"coma_y": +c}),
        (
            "Astig 0 +0.10, Coma X +0.06",
            {"astig_0": +0.10, "coma_x": +0.06},
        ),
    ]

    rows: list[dict[str, Any]] = []
    rendered: list[tuple[str, np.ndarray, dict[str, Any]]] = []
    for name, coefficients in cases:
        phase_waves = np.zeros_like(simulator.X, dtype=np.float32)
        for mode_name, coefficient in coefficients.items():
            phase_waves += np.float32(coefficient) * modes[mode_name]
        aberration = np.exp(
            1j * np.float32(2.0 * np.pi) * phase_waves
        ).astype(np.complex64)
        field = measured_amplitude * simulator.phase_factor * aberration
        raw = intensity(
            forward_fft(field.astype(np.complex64), np), np
        ).astype(np.float32)
        normalized, metrics, _ = evaluate_image(raw, simulator, experiment)
        rows.append(
            {
                "name": name,
                "coefficients_rms_waves": coefficients,
                "metrics": compact_metrics(metrics),
            }
        )
        rendered.append((name, normalized, metrics))

    x_um = simulator.base.x_um
    y_um = simulator.base.y_um
    xmask = np.abs(x_um) <= 230.0
    ymask = np.abs(y_um) <= 105.0
    xa = x_um[xmask]
    ya = y_um[ymask]
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
    for ax, (name, image, metrics) in zip(axes.ravel(), rendered):
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
        ax.set_title(
            f"{name} RMS waves\n"
            f"center={metrics['center_window_over_core_mean']:.3f}, "
            f"mid/sides={metrics['middle_third_over_side_thirds']:.3f}, "
            f"RMS={100.0 * metrics['core_rms_fraction']:.1f}%",
            fontsize=10,
        )
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(
        "Exploratory Noll-Zernike astigmatism and coma\n"
        "Measured Gaussian amplitude + fixed V2 phase; 15 mm pupil",
        fontsize=15,
    )
    fig.savefig(OUTPUT_FIGURE, dpi=180)
    plt.close(fig)

    payload = {
        "status": "exploratory_not_mainline",
        "phase_source": str(DEFAULT_CASE_DIR / "phase_refined.npy"),
        "pupil_diameter_mm": simulator.clear_aperture_m * 1e3,
        "definitions": {
            "astig_0": "sqrt(6)*(x_n^2-y_n^2)",
            "astig_45": "2*sqrt(6)*x_n*y_n",
            "coma_x": "sqrt(8)*(3*rho^2-2)*x_n",
            "coma_y": "sqrt(8)*(3*rho^2-2)*y_n",
        },
        "cases": rows,
    }
    (OUTPUT_DIR / "astigmatism_coma_exploratory_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    print(str(OUTPUT_FIGURE))


if __name__ == "__main__":
    main()
