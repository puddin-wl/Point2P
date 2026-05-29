"""Generate Romero-Dickey initial phase in Python (replicates MATLAB generate_initial_phase.m).

Usage:
  python make_phase0.py --beam 7 --out make_phase0_output/
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.special import erf

THIS_DIR = Path(__file__).resolve().parent


def make_grid(N: int, dx_doe_m: float):
    index = np.arange(-N / 2, N / 2)
    x_m = index * dx_doe_m
    y_m = index * dx_doe_m
    X_m, Y_m = np.meshgrid(x_m, y_m)
    grid = dict(N=N, dx_doe_m=dx_doe_m, extent_m=N * dx_doe_m,
                x_min_m=float(x_m[0]), x_max_m=float(x_m[-1]),
                y_min_m=float(y_m[0]), y_max_m=float(y_m[-1]))
    return X_m, Y_m, x_m, y_m, grid


def gaussian_input_field(X_m, Y_m, input_1e2_radius_m, aperture_radius_m):
    r2 = X_m**2 + Y_m**2
    w = input_1e2_radius_m
    amplitude = np.exp(-r2 / w**2)
    aperture_mask = r2 <= aperture_radius_m**2
    amplitude *= aperture_mask
    return amplitude.astype(np.complex128), aperture_mask


def romero_dickey_phase_1d(x_m, ri_m, Ro_m, lambda_m, f_m):
    xi = x_m / ri_m
    phi_dimless = xi * np.sqrt(np.pi) / 2.0 * erf(xi) + 0.5 * np.exp(-xi**2) - 0.5
    beta = 2.0 * np.pi * ri_m * Ro_m / (lambda_m * f_m)
    phase_rad = beta * phi_dimless
    return phase_rad, dict(ri_m=ri_m, Ro_m=Ro_m, beta=beta,
                           phi_dimless_min=float(phi_dimless.min()),
                           phi_dimless_max=float(phi_dimless.max()))


def fresnel_to_focal_fft(U_doe, dx_doe_m, lambda_m, f_m):
    N = U_doe.shape[0]
    df = 1.0 / (N * dx_doe_m)
    freq = (np.arange(N) - N / 2) * df
    focal_x_m = lambda_m * f_m * freq
    focal_y_m = focal_x_m.copy()
    U_focal = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(U_doe))) * dx_doe_m**2 / (1j * lambda_m * f_m)
    intensity = np.abs(U_focal)**2
    max_intensity = intensity.max()
    intensity_norm = intensity / max_intensity if max_intensity > 0 else intensity
    return focal_x_m, focal_y_m, intensity_norm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--beam", type=float, required=True, help="Beam 1/e^2 diameter in mm")
    parser.add_argument("--out", default=None, help="Output directory")
    parser.add_argument("--f", type=float, default=200e-3, help="Focal length in m")
    parser.add_argument("--N", type=int, default=2048)
    args = parser.parse_args()

    lambda_m = 532e-9
    f_m = args.f
    N = int(args.N)
    aperture_diameter_m = 15e-3
    aperture_radius_m = aperture_diameter_m / 2
    input_1e2_diameter_m = args.beam * 1e-3
    input_1e2_radius_m = input_1e2_diameter_m / 2
    input_1e_radius_m = input_1e2_radius_m / np.sqrt(2)
    target_size_x_m = 330e-6
    target_size_y_m = 120e-6
    requested_focal_dx_m = 2.5e-6

    doe_grid_extent_m = lambda_m * f_m / requested_focal_dx_m
    dx_doe_m = doe_grid_extent_m / N
    focal_dx_m = lambda_m * f_m / (N * dx_doe_m)
    Ro_x_m = target_size_x_m / np.sqrt(np.pi)
    Ro_y_m = target_size_y_m / np.sqrt(np.pi)
    beta_x = 2.0 * np.pi * input_1e_radius_m * Ro_x_m / (lambda_m * f_m)
    beta_y = 2.0 * np.pi * input_1e_radius_m * Ro_y_m / (lambda_m * f_m)

    print(f"f = {f_m*1e3:.0f} mm, beam = {args.beam} mm, N = {N}")
    print(f"dx_doe = {dx_doe_m*1e6:.4f} um, DOE extent = {doe_grid_extent_m*1e3:.2f} mm")
    print(f"focal_dx = {focal_dx_m*1e6:.4f} um")
    print(f"beta_x = {beta_x:.6f}, beta_y = {beta_y:.6f}")

    X_m, Y_m, x_m, y_m, grid = make_grid(N, dx_doe_m)
    input_field, aperture_mask = gaussian_input_field(X_m, Y_m, input_1e2_radius_m, aperture_radius_m)

    phase_x, info_x = romero_dickey_phase_1d(X_m, input_1e_radius_m, Ro_x_m, lambda_m, f_m)
    phase_y, info_y = romero_dickey_phase_1d(Y_m, input_1e_radius_m, Ro_y_m, lambda_m, f_m)
    phase_unwrapped = phase_x + phase_y
    phase_unwrapped[~aperture_mask] = 0.0
    phase_wrapped = np.mod(phase_unwrapped, 2.0 * np.pi)
    phase_wrapped[~aperture_mask] = np.nan

    unwrapped_min = phase_unwrapped[aperture_mask].min()
    unwrapped_max = phase_unwrapped[aperture_mask].max()
    print(f"phase unwrapped range in aperture: [{unwrapped_min:.6g}, {unwrapped_max:.6g}] rad")

    doe_field = input_field * np.exp(1j * phase_unwrapped)
    focal_x_m, focal_y_m, I_norm = fresnel_to_focal_fft(doe_field, dx_doe_m, lambda_m, f_m)

    out_dir = Path(args.out) if args.out else THIS_DIR / "make_phase0_output"
    out_dir.mkdir(parents=True, exist_ok=True)

    savemat(
        str(out_dir / "phase0.mat"),
        {
            "phase0_unwrapped_rad": phase_unwrapped.astype(np.float32),
            "phase0_wrapped_rad": phase_wrapped.astype(np.float32),
            "phase_info": {
                "method": "romero_dickey_separable",
                "phase_sign": 1,
                "phase_scale_x": 1,
                "phase_scale_y": 1,
                "unwrapped_min_rad": float(unwrapped_min),
                "unwrapped_max_rad": float(unwrapped_max),
                "x": info_x, "y": info_y,
            },
            "x_m": x_m.astype(np.float64),
            "y_m": y_m.astype(np.float64),
            "focal_x_m": focal_x_m.astype(np.float64),
            "focal_y_m": focal_y_m.astype(np.float64),
        },
        do_compression=True,
    )

    # Quick center profile check
    cy = N // 2
    prof_x = I_norm[cy, :]
    prof_y = I_norm[:, N // 2]
    def crossing_width(axis_um, profile, level):
        mid = len(profile) // 2
        if profile[mid] < level:
            return np.nan
        left = np.where(profile[:mid] >= level)[0]
        right = np.where(profile[mid:] < level)[0]
        if len(left) == 0 or len(right) == 0:
            return np.nan
        return float(axis_um[mid + right[0]] - axis_um[left[-1]])

    axis_um = focal_x_m * 1e6
    size50_x = crossing_width(axis_um, prof_x, 0.5)
    size50_y = crossing_width(axis_um, prof_y, 0.5)
    print(f"Initial RD phase size50_x: {size50_x:.1f} um, size50_y: {size50_y:.1f} um")
    print(f"Saved: {out_dir / 'phase0.mat'}")

    config_snapshot = dict(
        lambda_m=lambda_m, f_m=f_m, N=N,
        aperture_diameter_m=aperture_diameter_m,
        input_1e2_diameter_m=input_1e2_diameter_m,
        target_size_x_m=target_size_x_m, target_size_y_m=target_size_y_m,
        dx_doe_m=dx_doe_m, focal_dx_m=focal_dx_m,
        beta_x=beta_x, beta_y=beta_y,
        Ro_x_m=Ro_x_m, Ro_y_m=Ro_y_m,
    )
    savemat(str(out_dir / "config_snapshot.mat"), {"cfg": config_snapshot, "grid": grid}, do_compression=True)
    print(f"Saved: {out_dir / 'config_snapshot.mat'}")


if __name__ == "__main__":
    main()
