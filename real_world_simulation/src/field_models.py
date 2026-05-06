"""DOE-plane real input field models."""

from __future__ import annotations

from typing import Any

import numpy as np


def doe_axes(shape: tuple[int, int], dx_doe_m: float) -> tuple[np.ndarray, np.ndarray]:
    """Return DOE-plane x/y coordinate arrays in meters."""
    ny, nx = int(shape[0]), int(shape[1])
    dtype = np.float32
    x = (np.arange(nx, dtype=dtype) - nx // 2) * dtype(dx_doe_m)
    y = (np.arange(ny, dtype=dtype) - ny // 2) * dtype(dx_doe_m)
    return np.meshgrid(x, y)


def make_real_input_field(
    shape: tuple[int, int],
    dx_doe_m: float,
    wavelength_m: float,
    focal_length_m: float,
    params: dict[str, Any],
) -> tuple[np.ndarray, dict[str, float]]:
    """Build a normalized real input field in the DOE plane.

    The unclipped Gaussian power is normalized to 1 first. The circular clear
    aperture is then applied and reported as ``aperture_throughput_percent``.
    The transmitted field is normalized again before propagation so diffraction
    efficiency remains relative to the light that actually passed the aperture.
    """
    x, y = doe_axes(shape, dx_doe_m)
    offset_x_m = float(params.get("offset_x_mm", 0.0)) * 1e-3
    offset_y_m = float(params.get("offset_y_mm", 0.0)) * 1e-3
    diameter_x_m = float(params.get("diameter_1e2_x_mm", 5.0)) * 1e-3
    diameter_y_m = float(params.get("diameter_1e2_y_mm", 5.0)) * 1e-3
    aperture_radius_m = float(params.get("clear_aperture_mm", 15.0)) * 0.5e-3
    if diameter_x_m <= 0 or diameter_y_m <= 0:
        raise ValueError("Gaussian 1/e^2 diameters must be positive.")
    if aperture_radius_m <= 0:
        raise ValueError("Clear aperture diameter must be positive.")

    wx = np.float32(0.5 * diameter_x_m)
    wy = np.float32(0.5 * diameter_y_m)
    xr = x - np.float32(offset_x_m)
    yr = y - np.float32(offset_y_m)
    amp = np.exp(-((xr * xr) / (wx * wx) + (yr * yr) / (wy * wy))).astype(np.float32)
    pre_power = float(np.sum(amp * amp, dtype=np.float64))
    if pre_power <= 0:
        raise ValueError("Unclipped Gaussian has zero power.")
    amp *= np.float32(1.0 / np.sqrt(pre_power))

    aperture = (x * x + y * y) <= np.float32(aperture_radius_m * aperture_radius_m)
    amp = np.where(aperture, amp, np.float32(0.0)).astype(np.float32, copy=False)
    transmitted_power = float(np.sum(amp * amp, dtype=np.float64))
    if transmitted_power <= 0:
        raise ValueError("No input power passed the aperture.")
    aperture_throughput_percent = 100.0 * transmitted_power
    amp *= np.float32(1.0 / np.sqrt(transmitted_power))

    k = np.float32(2.0 * np.pi / wavelength_m)
    phase = np.zeros(shape, dtype=np.float32)

    divergence_edge_mrad = float(params.get("divergence_edge_mrad", 0.0))
    if divergence_edge_mrad:
        theta_edge = np.float32(divergence_edge_mrad * 1e-3)
        reference_radius = np.float32(0.25 * (diameter_x_m + diameter_y_m))
        if reference_radius <= 0:
            raise ValueError("Divergence reference radius must be positive.")
        phase += k * theta_edge * (xr * xr + yr * yr) / np.float32(2.0 * reference_radius)

    shift_x_m = float(params.get("pointing_shift_x_um", 0.0)) * 1e-6
    shift_y_m = float(params.get("pointing_shift_y_um", 0.0)) * 1e-6
    if shift_x_m or shift_y_m:
        theta_x = np.float32(shift_x_m / focal_length_m)
        theta_y = np.float32(shift_y_m / focal_length_m)
        phase += k * (theta_x * x + theta_y * y)

    field = amp.astype(np.complex64) * np.exp(1j * phase).astype(np.complex64)
    meta = {
        "aperture_throughput_percent": float(aperture_throughput_percent),
        "transmitted_power_before_renormalization": float(transmitted_power),
        "diameter_1e2_x_mm": float(params.get("diameter_1e2_x_mm", 5.0)),
        "diameter_1e2_y_mm": float(params.get("diameter_1e2_y_mm", 5.0)),
        "clear_aperture_mm": float(params.get("clear_aperture_mm", 15.0)),
    }
    return field, meta

