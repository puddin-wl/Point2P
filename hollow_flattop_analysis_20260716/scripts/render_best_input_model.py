"""Visualize the input-amplitude model that best reproduced the hollow flat top."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulate_v2_hollow_scan import (
    ANALYSIS_ROOT,
    DEFAULT_CASE_DIR,
    Simulator,
    normalize_power,
)


BEST_JSON = (
    ANALYSIS_ROOT
    / "results"
    / "02_elliptical_input_counterexample"
    / "best_elliptical_dip.json"
)
OUTPUT_DIR = BEST_JSON.parent


def main() -> None:
    best = json.loads(BEST_JSON.read_text(encoding="utf-8"))
    params = best["params"]
    simulator = Simulator(DEFAULT_CASE_DIR)
    depth = float(params["amplitude_dip_depth"])
    sigma_x_m = float(params["amplitude_dip_sigma_x_mm"]) * 1e-3
    sigma_y_m = float(params["amplitude_dip_sigma_y_mm"]) * 1e-3
    offset_x_m = float(params.get("amplitude_dip_offset_x_mm", 0.0)) * 1e-3
    offset_y_m = float(params.get("amplitude_dip_offset_y_mm", 0.0)) * 1e-3
    transmission_amplitude = 1.0 - depth * np.exp(
        -0.5
        * (
            (simulator.X - offset_x_m) ** 2 / sigma_x_m**2
            + (simulator.Y - offset_y_m) ** 2 / sigma_y_m**2
        )
    )
    modified_amplitude = normalize_power(
        (simulator.base_amplitude * transmission_amplitude).astype(np.float32)
    )
    base_intensity = simulator.base_amplitude.astype(np.float64) ** 2
    modified_intensity = modified_amplitude.astype(np.float64) ** 2
    base_intensity /= float(np.max(base_intensity))
    modified_intensity /= float(np.max(modified_intensity))
    x_mm = simulator.X[0, :].astype(np.float64) * 1e3
    y_mm = simulator.Y[:, 0].astype(np.float64) * 1e3
    cx = int(np.argmin(np.abs(x_mm)))
    cy = int(np.argmin(np.abs(y_mm)))
    view = 7.5
    xmask = np.abs(x_mm) <= view
    ymask = np.abs(y_mm) <= view

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for ax, image, title in (
        (axes[0, 0], base_intensity, "Ideal 6.5 mm Gaussian intensity"),
        (
            axes[0, 1],
            modified_intensity,
            "Modeled input intensity required by best hollow fit",
        ),
    ):
        roi = image[np.ix_(ymask, xmask)]
        im = ax.imshow(
            roi,
            origin="upper",
            extent=[-view, view, view, -view],
            cmap="turbo",
            vmin=0,
            vmax=1,
            aspect="equal",
        )
        ax.set_title(title)
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
        fig.colorbar(im, ax=ax, label="I / peak")

    axes[1, 0].plot(x_mm, base_intensity[cy, :], label="ideal Gaussian")
    axes[1, 0].plot(x_mm, modified_intensity[cy, :], label="modeled input")
    axes[1, 0].set_xlim(-7.5, 7.5)
    axes[1, 0].set_title("Input X center profile")
    axes[1, 0].set_xlabel("x / mm")
    axes[1, 0].set_ylabel("I / peak")
    axes[1, 0].grid(True, alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].plot(y_mm, base_intensity[:, cx], label="ideal Gaussian")
    axes[1, 1].plot(y_mm, modified_intensity[:, cx], label="modeled input")
    axes[1, 1].set_xlim(-7.5, 7.5)
    axes[1, 1].set_title("Input Y center profile")
    axes[1, 1].set_xlabel("y / mm")
    axes[1, 1].set_ylabel("I / peak")
    axes[1, 1].grid(True, alpha=0.25)
    axes[1, 1].legend()
    fig.suptitle(
        "Best phenomenological input model\n"
        f"amplitude dip depth={depth:.3f}, sigma_x={sigma_x_m*1e3:.3f} mm, "
        f"sigma_y={sigma_y_m*1e3:.3f} mm"
    )
    figure_path = OUTPUT_DIR / "best_input_model.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    center_transmission_amplitude = 1.0 - depth
    result = {
        "source_best_candidate": str(BEST_JSON),
        "phase_source": str(DEFAULT_CASE_DIR / "phase_refined.npy"),
        "model": params,
        "center_local_amplitude_transmission_before_power_renormalization": center_transmission_amplitude,
        "center_local_intensity_transmission_before_power_renormalization": center_transmission_amplitude**2,
        "normalized_input_center_intensity_over_peak": float(modified_intensity[cy, cx]),
        "note": (
            "This is a phenomenological input-amplitude model that reproduces the "
            "focal hollow with the fixed V2 phase. It is not yet proof that the real "
            "laser input has this exact centered elliptical depression."
        ),
        "output_figure": str(figure_path),
    }
    (OUTPUT_DIR / "best_input_model.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
