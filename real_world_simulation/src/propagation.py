"""Propagation utilities for fixed DOE simulation."""

from __future__ import annotations

import numpy as np


def forward_fft(field_in: np.ndarray) -> np.ndarray:
    """Centered orthonormal Fourier-lens propagation."""
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(field_in), norm="ortho"))


def intensity(field: np.ndarray) -> np.ndarray:
    """Return optical intensity as ``abs(field)**2``."""
    amp = np.abs(field)
    return (amp * amp).astype(np.float32, copy=False)


def angular_spectrum_defocus(
    field_at_focus: np.ndarray,
    focal_dx_m: float,
    focal_dy_m: float,
    wavelength_m: float,
    defocus_m: float,
) -> np.ndarray:
    """Propagate a focal-plane field by ``defocus_m`` using angular spectrum."""
    if abs(defocus_m) < 1e-15:
        return field_at_focus
    ny, nx = field_at_focus.shape
    fx = np.fft.fftfreq(nx, d=float(focal_dx_m)).astype(np.float32)
    fy = np.fft.fftfreq(ny, d=float(focal_dy_m)).astype(np.float32)
    fx2, fy2 = np.meshgrid(fx, fy)
    lam = np.float32(wavelength_m)
    arg = np.float32(1.0) - (lam * fx2) ** 2 - (lam * fy2) ** 2
    root = np.sqrt(np.maximum(arg, np.float32(0.0))).astype(np.float32)
    k = np.float32(2.0 * np.pi / wavelength_m)
    transfer = np.exp(1j * k * np.float32(defocus_m) * root).astype(np.complex64)
    spectrum = np.fft.fft2(np.fft.ifftshift(field_at_focus), norm="ortho")
    out = np.fft.fftshift(np.fft.ifft2(spectrum * transfer, norm="ortho"))
    return out.astype(np.complex64, copy=False)


def propagate_after_doe(
    field_after_doe: np.ndarray,
    focal_dx_m: float,
    focal_dy_m: float,
    wavelength_m: float,
    defocus_m: float = 0.0,
) -> np.ndarray:
    """Propagate from the DOE plane to the requested observation plane."""
    focal_field = forward_fft(field_after_doe)
    return angular_spectrum_defocus(
        focal_field,
        focal_dx_m=focal_dx_m,
        focal_dy_m=focal_dy_m,
        wavelength_m=wavelength_m,
        defocus_m=defocus_m,
    )

