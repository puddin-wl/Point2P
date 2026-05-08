"""Numerical initial phase for a truncated (non-Gaussian) input beam.

Uses stationary-phase energy conservation to compute a separable 1D phase
that maps an arbitrary input intensity profile to a raised-cosine flat-top.

Physical model:
  Expanded Gaussian (D_exp=25mm)
  → 5mm hard aperture
  → angular spectrum propagation (z=100mm)
  → DOE plane: extract 1D center-line intensity
  → numerical CDF matching → u(x) mapping
  → integrate → φ_x(x), φ_y(y)
  → build separable 2D phase φ(x,y) = φ_x(x) + φ_y(y)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = THIS_DIR.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "lab_test_f300mm"))

from src.backend import get_backend
from src.utils import ensure_unique_dir, timestamp, save_json


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
PHYSICAL = {
    "wavelength_m": 532e-9,
    "focal_length_m": 300e-3,
}
GRID = {
    "N": 2048,
    "focal_dx_um": 2.5,
    "focal_dy_um": 2.5,
    "dx_doe_m": 3.1171875e-05,
}
TARGET_CONFIG = {
    "W50_um": 330.0,
    "H50_um": 120.0,
    "delta_x_um": 15.0,
    "delta_y_um": 8.0,
    "guard_x_um": 20.0,
    "guard_y_um": 12.0,
}

D_EXP_MM = 25.0
Z_MM = 100.0
APERTURE_DIAMETER_MM = 5.0
CLEAR_APERTURE_MM = 15.0


# ---------------------------------------------------------------------------
# Copy of angular spectrum / truncated Gaussian (self-contained)
# ---------------------------------------------------------------------------
def angular_spectrum_propagate(field_in, dx_m, dy_m, wavelength_m, distance_m, xp):
    """Angular spectrum propagation between two parallel planes."""
    if abs(distance_m) < 1e-15:
        return field_in
    ny, nx = field_in.shape
    fx = xp.fft.fftfreq(nx, d=float(dx_m)).astype(xp.float32)
    fy = xp.fft.fftfreq(ny, d=float(dy_m)).astype(xp.float32)
    FX, FY = xp.meshgrid(fx, fy)
    lam = xp.float32(wavelength_m)
    arg = xp.float32(1.0) - (lam * FX) ** 2 - (lam * FY) ** 2
    root = xp.sqrt(xp.maximum(arg, xp.float32(0.0)))
    k = xp.float32(2.0 * xp.pi / wavelength_m)
    transfer = xp.exp(1j * k * xp.float32(distance_m) * root).astype(xp.complex64)
    spectrum = xp.fft.fft2(xp.fft.ifftshift(field_in), norm="ortho")
    out = xp.fft.fftshift(xp.fft.ifft2(spectrum * transfer, norm="ortho"))
    return out.astype(xp.complex64, copy=False)


def make_truncated_gaussian_at_aperture(shape, dx_m, D_exp_mm, aperture_diameter_mm, xp, dtype):
    """Gaussian truncated by hard circular aperture at aperture plane."""
    Ny, Nx = int(shape[0]), int(shape[1])
    x = (xp.arange(Nx, dtype=dtype) - (Nx // 2)) * dtype(dx_m)
    y = (xp.arange(Ny, dtype=dtype) - (Ny // 2)) * dtype(dx_m)
    X, Y = xp.meshgrid(x, y)
    r2 = X * X + Y * Y

    w = dtype(D_exp_mm * 1e-3 / 2.0)
    amp = xp.exp(-r2 / (w * w)).astype(dtype, copy=False)

    pre_power = dtype(xp.sum(amp * amp))
    amp = xp.where(pre_power > 0, amp / xp.sqrt(pre_power), amp)

    r_ap = dtype(aperture_diameter_mm * 1e-3 / 2.0)
    mask = r2 <= (r_ap * r_ap)
    amp = xp.where(mask, amp, dtype(0.0)).astype(dtype, copy=False)

    transmitted = dtype(xp.sum(amp * amp))
    throughput = 100.0 * float(transmitted) if float(transmitted) > 0 else 0.0
    amp = xp.where(transmitted > 0, amp / xp.sqrt(transmitted), amp)

    return amp.astype(xp.complex64, copy=False), throughput


# ---------------------------------------------------------------------------
# 1D raised-cosine edge (copied from rtad_target.py for self-containment)
# ---------------------------------------------------------------------------
def raised_cosine_edge(u, u0, u1):
    """1D raised-cosine falling edge: 1 for |u|≤u0, 0 for |u|≥u1."""
    u = np.abs(np.asarray(u, dtype=np.float64))
    C = np.zeros_like(u, dtype=np.float64)
    C[u <= u0] = 1.0
    idx = (u > u0) & (u < u1)
    t = (u[idx] - u0) / (u1 - u0)
    C[idx] = 0.5 * (1.0 + np.cos(np.pi * t))
    return C


# ---------------------------------------------------------------------------
# Numerical 1D phase via energy conservation
# ---------------------------------------------------------------------------
def compute_1d_phase_from_profiles(
    x_doe_m: np.ndarray,
    I_in_1d: np.ndarray,
    u_focal_um: np.ndarray,
    I_out_1d: np.ndarray,
    wavelength_m: float,
    focal_length_m: float,
    target_edge_um: float,
) -> np.ndarray:
    """Compute 1D DOE phase from input/output intensity profiles.

    Uses stationary-phase energy conservation:
      dφ/dx = (2π/(λf)) · u(x)
    where CDF_in(x) = CDF_out(u).

    ``target_edge_um`` clamps the mapping to the physical target support
    (e.g. ~200 um = a1 + guard_x).

    Returns φ(x) in radians on the same grid as x_doe_m.
    """
    x = np.asarray(x_doe_m, dtype=np.float64).ravel()
    I_in = np.asarray(I_in_1d, dtype=np.float64).ravel()
    u = np.asarray(u_focal_um, dtype=np.float64).ravel()
    I_out = np.asarray(I_out_1d, dtype=np.float64).ravel()

    eps = 1e-30
    cdf_in = np.cumsum(np.maximum(I_in, 0))
    cdf_in /= max(cdf_in[-1], eps)

    cdf_out = np.cumsum(np.maximum(I_out, 0))
    cdf_out /= max(cdf_out[-1], eps)

    # Restrict interpolation to the non-flat region of cdf_out
    active = (cdf_out > 1e-12) & (cdf_out < 1.0 - 1e-12)
    if np.any(active):
        u_active = u[active]
        cdf_active = cdf_out[active]
        # Ensure strictly increasing by adding a tiny ramp
        cdf_active = np.maximum.accumulate(cdf_active)
    else:
        u_active = u
        cdf_active = cdf_out

    u_map = np.interp(cdf_in, cdf_active, u_active)

    # Clamp to target support; outside the aperture there is no light,
    # so the mapping is pinned to the target edge.
    edge = abs(float(target_edge_um))
    u_map = np.clip(u_map, -edge, edge)

    dphi_dx = (2.0 * np.pi / (wavelength_m * focal_length_m)) * (u_map * 1e-6)
    phi = np.zeros_like(x)
    dx = x[1] - x[0]
    mid = len(x) // 2
    for i in range(mid + 1, len(x)):
        phi[i] = phi[i - 1] + 0.5 * dx * (dphi_dx[i] + dphi_dx[i - 1])
    for i in range(mid - 1, -1, -1):
        phi[i] = phi[i + 1] - 0.5 * dx * (dphi_dx[i] + dphi_dx[i + 1])

    return phi.astype(np.float32)


# ---------------------------------------------------------------------------
# Forward FFT (copy from lab_test_f300mm/src/propagation.py for self-containment)
# ---------------------------------------------------------------------------
def forward_fft(field_in, xp):
    return xp.fft.fftshift(xp.fft.fft2(xp.fft.ifftshift(field_in), norm="ortho"))


def intensity(field, xp):
    amp = xp.abs(field)
    return amp * amp


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Numerical initial phase for truncated beam")
    parser.add_argument("--cpu", action="store_true", help="Force CPU backend")
    parser.add_argument("--D-exp", type=float, default=D_EXP_MM)
    parser.add_argument("--z-mm", type=float, default=Z_MM)
    args = parser.parse_args()

    D_exp_mm = args.D_exp
    z_mm = args.z_mm

    backend = get_backend(use_cupy=not args.cpu, verbose=True)
    xp = backend.xp
    dtype = backend.float_dtype

    # ------------------------------------------------------------------
    # 1. Compute truncated beam field at DOE plane
    # ------------------------------------------------------------------
    print(f"\nComputing truncated beam field (D_exp={D_exp_mm}mm, z={z_mm}mm)...")
    N = int(GRID["N"])
    dx_doe_m = float(GRID["dx_doe_m"])
    shape = (N, N)
    wavelength_m = float(PHYSICAL["wavelength_m"])

    field_ap, throughput = make_truncated_gaussian_at_aperture(
        shape, dx_doe_m, D_exp_mm, APERTURE_DIAMETER_MM, xp, dtype
    )
    print(f"  Aperture throughput: {throughput:.2f}%")

    if z_mm > 1e-6:
        field_doe = angular_spectrum_propagate(
            field_ap, dx_doe_m, dx_doe_m, wavelength_m, z_mm * 1e-3, xp
        )
    else:
        field_doe = field_ap

    # Apply 15mm clear aperture
    x_doe = (xp.arange(N, dtype=dtype) - (N // 2)) * dtype(dx_doe_m)
    y_doe = (xp.arange(N, dtype=dtype) - (N // 2)) * dtype(dx_doe_m)
    Xd, Yd = xp.meshgrid(x_doe, y_doe)
    r2_doe = Xd * Xd + Yd * Yd
    r_clear = dtype(CLEAR_APERTURE_MM * 1e-3 / 2.0)
    clear_mask = r2_doe <= (r_clear * r_clear)
    field_doe = xp.where(clear_mask, field_doe, dtype(0.0) + 0j * dtype(0.0))
    field_doe = field_doe.astype(xp.complex64)

    # Renormalise power after clear aperture
    p_doe = dtype(xp.sum(xp.abs(field_doe) ** 2))
    if float(p_doe) > 0:
        field_doe = field_doe / xp.sqrt(p_doe)

    # Extract amplitude and intensity
    amp_doe = xp.abs(field_doe)
    I_doe = amp_doe * amp_doe

    # Move to NumPy
    if hasattr(I_doe, "get"):
        I_doe_np = I_doe.get()
        amp_doe_np = amp_doe.get()
        field_doe_np = field_doe.get()
    else:
        I_doe_np = np.asarray(I_doe)
        amp_doe_np = np.asarray(amp_doe)
        field_doe_np = np.asarray(field_doe)

    # DOE-plane coordinate axes in meters
    x_doe_np = (np.arange(N, dtype=np.float64) - N // 2) * dx_doe_m
    y_doe_np = (np.arange(N, dtype=np.float64) - N // 2) * dx_doe_m

    # ------------------------------------------------------------------
    # 2. Extract 1D marginal intensity profiles
    #    (integrate over the orthogonal direction for energy conservation)
    # ------------------------------------------------------------------
    I_x = np.sum(I_doe_np, axis=0)   # ∫ I(x,y) dy  → marginal along x
    I_y = np.sum(I_doe_np, axis=1)   # ∫ I(x,y) dx  → marginal along y

    print(f"  DOE-plane I_x (marginal): max={I_x.max():.6f}, sum={I_x.sum():.6f}")
    print(f"  DOE-plane I_y (marginal): max={I_y.max():.6f}, sum={I_y.sum():.6f}")

    # ------------------------------------------------------------------
    # 3. Define 1D RTAD output target profiles
    # ------------------------------------------------------------------
    W50 = float(TARGET_CONFIG["W50_um"])
    H50 = float(TARGET_CONFIG["H50_um"])
    delta_x = float(TARGET_CONFIG["delta_x_um"])
    delta_y = float(TARGET_CONFIG["delta_y_um"])
    guard_x = float(TARGET_CONFIG["guard_x_um"])
    guard_y = float(TARGET_CONFIG["guard_y_um"])

    a50 = W50 / 2.0  # 165 um
    b50 = H50 / 2.0  # 60 um
    a0 = a50 - delta_x  # 150 um
    a1 = a50 + delta_x  # 180 um
    b0 = b50 - delta_y  # 52 um
    b1 = b50 + delta_y  # 68 um
    a2 = a1 + guard_x   # 200 um
    b2 = b1 + guard_y   # 80 um

    focal_dx = float(GRID["focal_dx_um"])
    u_axis = (np.arange(N, dtype=np.float64) - N // 2) * focal_dx  # focal-plane axis in um
    v_axis = (np.arange(N, dtype=np.float64) - N // 2) * focal_dx

    I_out_x = raised_cosine_edge(u_axis, a0, a1)  # 1D x-target at y=0 (C_y(|0|)=1)
    I_out_y = raised_cosine_edge(v_axis, b0, b1)  # 1D y-target at x=0 (C_x(|0|)=1)

    print(f"  Target x: flat=±{a0}um, edge=±{a1}um, guard=±{a2}um")
    print(f"  Target y: flat=±{b0}um, edge=±{b1}um, guard=±{b2}um")

    # ------------------------------------------------------------------
    # 4. Compute numerical 1D phases
    # ------------------------------------------------------------------
    print("\nComputing numerical 1D phases via energy conservation...")
    phi_x = compute_1d_phase_from_profiles(
        x_doe_np, I_x, u_axis, I_out_x, wavelength_m, float(PHYSICAL["focal_length_m"]),
        target_edge_um=a2,
    )
    phi_y = compute_1d_phase_from_profiles(
        y_doe_np, I_y, v_axis, I_out_y, wavelength_m, float(PHYSICAL["focal_length_m"]),
        target_edge_um=b2,
    )

    phi_x_range = float(phi_x.max() - phi_x.min())
    phi_y_range = float(phi_y.max() - phi_y.min())
    print(f"  φ_x: min={phi_x.min():.2f}, max={phi_x.max():.2f} rad")
    print(f"  φ_y: min={phi_y.min():.2f}, max={phi_y.max():.2f} rad")
    print(f"  φ_x range: {phi_x_range:.2f} rad ({phi_x_range/(2*np.pi):.1f} × 2π)")
    print(f"  φ_y range: {phi_y_range:.2f} rad ({phi_y_range/(2*np.pi):.1f} × 2π)")

    # ------------------------------------------------------------------
    # 5. Build separable 2D phase
    # ------------------------------------------------------------------
    X_phi, Y_phi = np.meshgrid(phi_x, phi_y)  # rows=y, cols=x
    phase_2d = (X_phi + Y_phi).astype(np.float32)

    # Wrap to [0, 2π)
    phase_wrapped = np.mod(phase_2d, np.float32(2.0 * np.pi))

    print(f"\n  Phase 2D: shape={phase_wrapped.shape}, "
          f"min={phase_wrapped.min():.4f}, max={phase_wrapped.max():.4f}")

    # ------------------------------------------------------------------
    # 6. Quick forward propagation to check RMS
    # ------------------------------------------------------------------
    print("\n=== Quick check: forward propagation with numerical phase ===")

    # Build DOE-plane field: truncated amplitude × exp(i * phase_2d)
    # Note: diffraction phase from angular spectrum is included in amp_doe_np (complex)
    # But for a pure phase DOE, we use a real amplitude + our computed phase.
    # The diffraction phase between aperture and DOE is small (near-field propagation
    # with Fresnel number ~117) and we absorb it into the phase.
    amp_input = np.abs(field_doe_np).astype(np.float32)
    phase_dev = xp.asarray(phase_wrapped, dtype=xp.float32)
    amp_dev = xp.asarray(amp_input, dtype=xp.float32)

    doe_field = amp_dev.astype(xp.complex64) * xp.exp(1j * phase_dev)
    focal_field = forward_fft(doe_field, xp)
    I_focal = intensity(focal_field, xp)

    if hasattr(I_focal, "get"):
        I_focal_np = I_focal.get()
    else:
        I_focal_np = np.asarray(I_focal)

    # Compute quick RMS over flat region
    # mask_flat: 2D array, rows=y, cols=x
    XX, YY = np.meshgrid(u_axis, v_axis)
    mask_flat = (np.abs(XX) <= a0) & (np.abs(YY) <= b0)

    flat_vals = I_focal_np[mask_flat]
    if flat_vals.size > 0:
        flat_mean = np.mean(flat_vals)
        flat_std = np.std(flat_vals)
        rms_pct = 100.0 * flat_std / flat_mean if flat_mean > 0 else np.nan
    else:
        rms_pct = np.nan

    # Center profiles
    I_norm = I_focal_np / flat_mean if flat_mean > 0 else I_focal_np
    cy_f = I_focal_np.shape[0] // 2
    cx_f = I_focal_np.shape[1] // 2
    prof_x = I_norm[cy_f, :]
    prof_y = I_norm[:, cx_f]

    # size50 via threshold crossings
    def crossing_width(axis_um, profile, level):
        """Find width at a given level via linear interpolation."""
        mid = len(profile) // 2
        if profile[mid] < level:
            return np.nan
        # left crossing
        left = np.nan
        for i in range(mid, 0, -1):
            if (profile[i-1] <= level <= profile[i]) or (profile[i] <= level <= profile[i-1]):
                t = (level - profile[i-1]) / (profile[i] - profile[i-1]) if abs(profile[i] - profile[i-1]) > 1e-30 else 0.5
                left = axis_um[i-1] + t * (axis_um[i] - axis_um[i-1])
                break
        # right crossing
        right = np.nan
        for i in range(mid, len(profile) - 1):
            if (profile[i] >= level >= profile[i+1]) or (profile[i+1] >= level >= profile[i]):
                t = (level - profile[i]) / (profile[i+1] - profile[i]) if abs(profile[i+1] - profile[i]) > 1e-30 else 0.5
                right = axis_um[i] + t * (axis_um[i+1] - axis_um[i])
                break
        if np.isfinite(left) and np.isfinite(right):
            return float(right - left)
        return np.nan

    size50_x = crossing_width(u_axis, prof_x, 0.5)
    size50_y = crossing_width(v_axis, prof_y, 0.5)

    # e⁻² efficiency
    size13x = crossing_width(u_axis, prof_x, np.exp(-2))
    size13y = crossing_width(v_axis, prof_y, np.exp(-2))
    total_power = float(np.sum(I_focal_np))
    if np.isfinite(size13x) and np.isfinite(size13y) and total_power > 0:
        x_roi = np.abs(XX) <= size13x / 2.0
        y_roi = np.abs(YY) <= size13y / 2.0
        eff = 100.0 * float(np.sum(I_focal_np[y_roi & x_roi])) / total_power
    else:
        eff = np.nan

    print(f"  RMS nonuniformity: {rms_pct:.4f}%")
    print(f"  size50_x: {size50_x:.1f} um  (target: {W50} um)")
    print(f"  size50_y: {size50_y:.1f} um  (target: {H50} um)")
    print(f"  e^-2 efficiency: {eff:.2f}%")

    # ------------------------------------------------------------------
    # 7. Save
    # ------------------------------------------------------------------
    stamp = timestamp()
    outdir = ensure_unique_dir(THIS_DIR / "artifacts" / f"{stamp}_numerical_phase")
    print(f"\nOutput: {outdir}")

    np.save(outdir / "phase_numerical.npy", phase_wrapped)
    np.save(outdir / "phase_numerical_unwrapped.npy", phase_2d)
    np.save(outdir / "amplitude_at_doe.npy", amp_input)
    np.save(outdir / "field_at_doe.npy", field_doe_np)
    np.save(outdir / "phi_x_1d.npy", phi_x)
    np.save(outdir / "phi_y_1d.npy", phi_y)
    np.savez(outdir / "focal_check.npz",
             I_focal=I_focal_np, I_norm=I_norm,
             prof_x=prof_x, prof_y=prof_y,
             u_axis=u_axis, v_axis=v_axis)

    # Save as MAT for MATLAB compatibility
    try:
        from scipy.io import savemat
        savemat(outdir / "phase_numerical.mat",
                {"phase_numerical_wrapped_rad": phase_wrapped,
                 "phase_numerical_unwrapped_rad": phase_2d},
                do_compression=True)
    except Exception:
        pass

    config = {
        "physical": PHYSICAL,
        "grid": GRID,
        "target": TARGET_CONFIG,
        "D_exp_mm": D_exp_mm,
        "z_mm": z_mm,
        "aperture_diameter_mm": APERTURE_DIAMETER_MM,
        "clear_aperture_mm": CLEAR_APERTURE_MM,
        "rms_numerical_phase_pct": float(rms_pct) if np.isfinite(rms_pct) else None,
        "size50_x_um": size50_x,
        "size50_y_um": size50_y,
        "efficiency_e2_pct": eff,
        "phi_x_range_rad": phi_x_range,
        "phi_y_range_rad": phi_y_range,
    }
    save_json(outdir / "config.json", config)

    print("Done.")


if __name__ == "__main__":
    main()
