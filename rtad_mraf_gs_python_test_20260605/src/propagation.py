"""Fourier-lens propagation utilities.

Matrix convention throughout this project:
    - axis 0 / rows are y
    - axis 1 / columns are x
"""

from __future__ import annotations

from typing import Any


def forward_fft(field_in: Any, xp: Any) -> Any:
    """Fraunhofer forward propagation through a Fourier lens.

    Uses an orthonormal FFT so total discrete power is preserved by the
    transform. The output is shifted so the DC/zero spatial frequency is at
    the center pixel.
    """
    return xp.fft.fftshift(xp.fft.fft2(xp.fft.ifftshift(field_in), norm="ortho"))


def backward_fft(field_out: Any, xp: Any) -> Any:
    """Inverse Fourier-lens propagation back to the DOE/input plane."""
    return xp.fft.fftshift(xp.fft.ifft2(xp.fft.ifftshift(field_out), norm="ortho"))


def intensity(field: Any, xp: Any) -> Any:
    """Return optical intensity ``abs(field)**2``."""
    amp = xp.abs(field)
    return amp * amp


def l2_norm(array: Any, xp: Any, eps: float = 1e-20) -> Any:
    """Return ``sqrt(sum(abs(array)**2))`` with a small floor."""
    value = xp.sqrt(xp.sum(xp.abs(array) ** 2))
    return xp.maximum(value, eps)


def normalize_power(array: Any, xp: Any, target_power: float = 1.0, eps: float = 1e-20) -> Any:
    """Scale an amplitude or field array to a requested total power."""
    norm = l2_norm(array, xp, eps=eps)
    return array * (target_power ** 0.5 / norm)


def normalize_mean_in_mask(values: Any, mask: Any, xp: Any, target_mean: float = 1.0, eps: float = 1e-20) -> Any:
    """Scale ``values`` so the mean inside ``mask`` equals ``target_mean``."""
    mean = xp.mean(values[mask])
    return values * (target_mean / xp.maximum(mean, eps))


def make_input_gaussian(
    shape: tuple[int, int],
    dx_doe_m: float,
    gaussian_1e2_diameter_m: float | None,
    clear_aperture_m: float,
    xp: Any,
    dtype: Any,
    gaussian_1e2_diameter_x_m: float | None = None,
    gaussian_1e2_diameter_y_m: float | None = None,
) -> Any:
    """Build the normalized input amplitude in the DOE plane.

    Diameters are 1/e^2 intensity diameters. The optional x/y values support
    an elliptical Gaussian. Each missing axis falls back to the legacy scalar
    ``gaussian_1e2_diameter_m``. For radii ``wx`` and ``wy``, the amplitude is
    ``A = exp(-(x^2/wx^2 + y^2/wy^2))``. A hard circular clear aperture is
    then applied.
    """
    diameter_x_m = gaussian_1e2_diameter_x_m
    diameter_y_m = gaussian_1e2_diameter_y_m
    if diameter_x_m is None:
        diameter_x_m = gaussian_1e2_diameter_m
    if diameter_y_m is None:
        diameter_y_m = gaussian_1e2_diameter_m
    if diameter_x_m is None or diameter_y_m is None:
        raise ValueError("Gaussian 1/e^2 diameter must be supplied for both x and y axes")
    if diameter_x_m <= 0 or diameter_y_m <= 0:
        raise ValueError("Gaussian 1/e^2 diameters must be positive")

    Ny, Nx = int(shape[0]), int(shape[1])
    x = (xp.arange(Nx, dtype=dtype) - (Nx // 2)) * dtype(dx_doe_m)
    y = (xp.arange(Ny, dtype=dtype) - (Ny // 2)) * dtype(dx_doe_m)
    X, Y = xp.meshgrid(x, y)
    r2 = X * X + Y * Y
    wx = dtype(diameter_x_m / 2.0)
    wy = dtype(diameter_y_m / 2.0)
    aperture_radius = dtype(clear_aperture_m / 2.0)
    amp = xp.exp(-(X * X / (wx * wx) + Y * Y / (wy * wy))).astype(dtype, copy=False)
    amp = xp.where(r2 <= aperture_radius * aperture_radius, amp, dtype(0.0))
    return normalize_power(amp, xp, target_power=1.0).astype(dtype, copy=False)
