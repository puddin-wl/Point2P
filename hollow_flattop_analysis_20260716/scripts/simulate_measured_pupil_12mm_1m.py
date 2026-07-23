"""Compare no pupil with a 12 mm field-lens pupil at a 1 m distance.

The input is the full-frame measured 2026-07-16 Gaussian amplitude carrying
the fixed V2 phase. No blaze, Zernike aberration, or other wavefront term is
added. After 1 m angular-spectrum propagation, the field is optionally clipped
by a centered 12 mm circular pupil and Fourier transformed by the 429 mm field
lens.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scan_slm_to_field_lens_free_space import (
    ANALYSIS_ROOT,
    PAD_N,
    SLM_PITCH_M,
    WAVELENGTH_M,
    embed_slm,
    load_phase_variants,
    measured_amplitude_on_slm,
    slm_coordinates,
)
from simulate_field_lens_pupil_fourier import (
    FIELD_LENS_FOCAL_LENGTH_M,
    forward_fourier,
    normalize_to_core,
)
from simulate_measured_pupil_13_14mm import evaluate_pupil


OUTPUT_DIR = ANALYSIS_ROOT / "results" / "12_measured_pupil_12mm_1m"
OUTPUT_FIGURE = OUTPUT_DIR / "measured_pupil_12mm_z1m_comparison.png"
DISTANCE_M = 1.0
PUPIL_DIAMETER_MM = 12.0


def propagate_to_field_lens(field0: np.ndarray, z_m: float) -> np.ndarray:
    """Angular-spectrum propagation with the same sampling as prior scans."""
    frequency = np.fft.fftfreq(PAD_N, d=SLM_PITCH_M)
    fx, fy = np.meshgrid(frequency, frequency)
    k = 2.0 * np.pi / WAVELENGTH_M
    angular_argument = np.maximum(
        0.0,
        1.0 - (WAVELENGTH_M * fx) ** 2 - (WAVELENGTH_M * fy) ** 2,
    )
    delta_kz = (k * (np.sqrt(angular_argument) - 1.0)).astype(np.float32)
    spectrum = np.fft.fft2(field0, norm="ortho").astype(np.complex64)
    transfer = np.exp(1j * delta_kz * z_m).astype(np.complex64)
    return np.fft.ifft2(spectrum * transfer, norm="ortho").astype(
        np.complex64
    )


def plot_comparison(
    focal_no_pupil: np.ndarray,
    focal_clipped: np.ndarray,
    x_focal_um: np.ndarray,
    y_focal_um: np.ndarray,
    metrics: dict[str, float],
) -> None:
    """Render one comparison image with maps and central profiles."""
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 8.2),
        constrained_layout=True,
    )
    xmask = np.abs(x_focal_um) <= 250.0
    ymask = np.abs(y_focal_um) <= 130.0
    xa = x_focal_um[xmask]
    ya = y_focal_um[ymask]

    for ax, image, title in (
        (
            axes[0, 0],
            focal_no_pupil,
            (
                "No pupil\n"
                f"center={metrics['no_pupil_center_ratio']:.4f}, "
                f"RMS={100.0 * metrics['no_pupil_core_rms_fraction']:.2f}%"
            ),
        ),
        (
            axes[0, 1],
            focal_clipped,
            (
                "12 mm circular pupil\n"
                f"center={metrics['clipped_center_ratio']:.4f}, "
                f"RMS={100.0 * metrics['clipped_core_rms_fraction']:.2f}%"
            ),
        ),
    ):
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
        ax.set_title(title)
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

    x_center = int(np.argmin(np.abs(x_focal_um)))
    y_center = int(np.argmin(np.abs(y_focal_um)))
    axes[1, 0].plot(
        x_focal_um[xmask],
        focal_no_pupil[y_center, xmask],
        label="no pupil",
    )
    axes[1, 0].plot(
        x_focal_um[xmask],
        focal_clipped[y_center, xmask],
        label="12 mm pupil",
    )
    axes[1, 0].set_title("Horizontal center profile")
    axes[1, 0].set_xlabel("x / um")
    axes[1, 0].set_ylabel("normalized intensity")
    axes[1, 0].set_ylim(0.0, 1.6)
    axes[1, 0].grid(True, alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].plot(
        y_focal_um[ymask],
        focal_no_pupil[ymask, x_center],
        label="no pupil",
    )
    axes[1, 1].plot(
        y_focal_um[ymask],
        focal_clipped[ymask, x_center],
        label="12 mm pupil",
    )
    axes[1, 1].set_title("Vertical center profile")
    axes[1, 1].set_xlabel("y / um")
    axes[1, 1].set_ylabel("normalized intensity")
    axes[1, 1].set_ylim(0.0, 1.6)
    axes[1, 1].grid(True, alpha=0.25)
    axes[1, 1].legend()

    fig.suptitle(
        "Measured Gaussian + V2, 1.0 m propagation: 12 mm pupil test\n"
        f"power removed={metrics['pupil_power_loss_percent']:.4f}%, "
        f"center change={metrics['center_ratio_change']:+.4f}, "
        f"middle/sides "
        f"{metrics['no_pupil_middle_over_sides']:.4f}"
        f" -> {metrics['clipped_middle_over_sides']:.4f}, "
        f"50% IoU={metrics['footprint_iou_clipped_vs_no_pupil']:.4f}"
    )
    fig.savefig(OUTPUT_FIGURE, dpi=180)
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    _, _, x_slm, y_slm = slm_coordinates()
    measured_amplitude = measured_amplitude_on_slm(x_slm, y_slm)
    phase_v2 = load_phase_variants()["v2_no_blaze"]
    field_slm = measured_amplitude * np.exp(1j * phase_v2).astype(
        np.complex64
    )
    field0 = embed_slm(field_slm.astype(np.complex64))
    field_lens = propagate_to_field_lens(field0, DISTANCE_M)

    axis_pad = (
        np.arange(PAD_N, dtype=np.float32) - PAD_N // 2
    ) * np.float32(SLM_PITCH_M)
    x_pad, y_pad = np.meshgrid(axis_pad, axis_pad)
    pupil_radius_m = np.float32(PUPIL_DIAMETER_MM * 0.5e-3)
    pupil_mask = x_pad * x_pad + y_pad * y_pad <= pupil_radius_m**2

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
    focal_no_pupil = normalize_to_core(
        np.abs(forward_fourier(field_lens)) ** 2,
        x_focal_um,
        y_focal_um,
    )
    metrics, focal_clipped = evaluate_pupil(
        field_lens,
        pupil_mask,
        PUPIL_DIAMETER_MM,
        DISTANCE_M,
        focal_no_pupil,
        x_focal_um,
        y_focal_um,
    )
    plot_comparison(
        focal_no_pupil,
        focal_clipped,
        x_focal_um,
        y_focal_um,
        metrics,
    )
    print(json.dumps(metrics, ensure_ascii=False, indent=2))
    print(str(OUTPUT_FIGURE))


if __name__ == "__main__":
    main()
