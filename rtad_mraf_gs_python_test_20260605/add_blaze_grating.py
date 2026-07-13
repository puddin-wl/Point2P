"""Add blazed grating to DOE phase for focal-spot shift.

A linear phase ramp (blazed grating) in the DOE plane translates the focal spot:
    Δx_focal = λf / T_x,   φ_blaze(x) = 2π · x · Δx / (λf)

Usage:
  python add_blaze_grating.py <phase.npy> --shift-x 200 --shift-y 200 --out <output_dir>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat


def add_blaze_grating(
    phase: np.ndarray,
    dx_doe_m: float,
    wavelength_m: float,
    focal_length_m: float,
    shift_x_um: float,
    shift_y_um: float,
) -> np.ndarray:
    """Add blazed grating to shift focal spot by (shift_x_um, shift_y_um).

    Parameters
    ----------
    phase : np.ndarray
        Input DOE phase in radians, shape (Ny, Nx).
    dx_doe_m : float
        DOE pixel size in meters.
    wavelength_m : float
        Laser wavelength in meters.
    focal_length_m : float
        Lens focal length in meters.
    shift_x_um : float
        Desired X shift in focal plane (μm).
    shift_y_um : float
        Desired Y shift in focal plane (μm).

    Returns
    -------
    np.ndarray
        New phase with blazed grating added, wrapped to [0, 2π).
    """
    Ny, Nx = phase.shape
    shift_x_m = shift_x_um * 1e-6
    shift_y_m = shift_y_um * 1e-6

    # Coordinate axes (centered, in meters)
    x_m = (np.arange(Nx, dtype=np.float64) - Nx // 2) * dx_doe_m
    y_m = (np.arange(Ny, dtype=np.float64) - Ny // 2) * dx_doe_m
    X_m, Y_m = np.meshgrid(x_m, y_m)

    # Blazed grating phase: φ = 2π/(λf) · (Δx·x + Δy·y)
    k_factor = 2.0 * np.pi / (wavelength_m * focal_length_m)
    blaze = k_factor * (shift_x_m * X_m + shift_y_m * Y_m)

    phase_new = phase + blaze.astype(phase.dtype)
    phase_new = np.mod(phase_new, 2.0 * np.pi)

    return phase_new


def main():
    parser = argparse.ArgumentParser(
        description="Add blazed grating to shift focal spot position."
    )
    parser.add_argument("phase", help="Path to refined phase .npy file")
    parser.add_argument("--shift-x", type=float, required=True,
                        help="X shift in focal plane (μm)")
    parser.add_argument("--shift-y", type=float, required=True,
                        help="Y shift in focal plane (μm)")
    parser.add_argument("--f", type=float, default=0.429,
                        help="Focal length in m (default: 0.429)")
    parser.add_argument("--wavelength", type=float, default=532e-9,
                        help="Wavelength in m (default: 532e-9)")
    parser.add_argument("--focal-dx", type=float, default=2.5,
                        help="Focal-plane sampling in μm (default: 2.5)")
    parser.add_argument("--N", type=int, default=2048,
                        help="Grid size (default: 2048)")
    parser.add_argument("--out", required=True,
                        help="Output directory")
    parser.add_argument("--verify", action="store_true",
                        help="Verify shift via FFT propagation")
    args = parser.parse_args()

    phase_path = Path(args.phase)
    if not phase_path.exists():
        print(f"ERROR: phase file not found: {phase_path}")
        sys.exit(1)

    phase = np.load(str(phase_path)).astype(np.float64)
    print(f"Loaded phase: {phase.shape}, "
          f"[{phase.min():.4f}, {phase.max():.4f}]")

    # Compute dx_doe from focal sampling
    lambda_m = args.wavelength
    f_m = args.f
    N = args.N
    focal_dx_m = args.focal_dx * 1e-6
    doe_extent_m = lambda_m * f_m / focal_dx_m
    dx_doe_m = doe_extent_m / N

    print(f"λ = {lambda_m*1e9:.0f} nm, f = {f_m*1e3:.0f} mm")
    print(f"dx_doe = {dx_doe_m*1e6:.4f} μm, DOE extent = {doe_extent_m*1e3:.2f} mm")
    print(f"Shift: Δx = {args.shift_x} μm, Δy = {args.shift_y} μm")

    # Grating period in pixels
    T_x_px = lambda_m * f_m / (dx_doe_m * args.shift_x * 1e-6)
    T_y_px = lambda_m * f_m / (dx_doe_m * args.shift_y * 1e-6)
    print(f"Grating period: T_x = {T_x_px:.1f} px, T_y = {T_y_px:.1f} px")
    print(f"Phase ramp: {360.0/T_x_px:.2f} °/px (x), "
          f"{360.0/T_y_px:.2f} °/px (y)")

    phase_blazed = add_blaze_grating(
        phase, dx_doe_m, lambda_m, f_m,
        args.shift_x, args.shift_y,
    )

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save outputs
    np.save(str(out_dir / "phase_blazed.npy"), phase_blazed.astype(np.float32))
    savemat(
        str(out_dir / "phase_blazed.mat"),
        {"phase_blazed_rad": phase_blazed.astype(np.float64)},
        do_compression=True,
    )
    print(f"Saved: {out_dir / 'phase_blazed.npy'}")
    print(f"Saved: {out_dir / 'phase_blazed.mat'}")

    # ---- verification via FFT ----
    if args.verify:
        print("\n--- Verification: FFT propagation ---")
        from src.propagation import forward_fft, intensity, make_input_gaussian

        # Build input amplitude
        amp = make_input_gaussian(
            shape=phase.shape,
            dx_doe_m=dx_doe_m,
            gaussian_1e2_diameter_m=6.5e-3,  # from the run
            clear_aperture_m=15e-3,
            xp=np,
            dtype=np.float32,
        )

        def propagate(ph):
            field = amp * np.exp(1j * ph).astype(np.complex64)
            return intensity(forward_fft(field, np), np)

        I_orig = propagate(phase)
        I_blazed = propagate(phase_blazed)

        # Focal-plane axes
        df = 1.0 / (N * dx_doe_m)
        freq = (np.arange(N) - N / 2) * df
        focal_x = lambda_m * f_m * freq * 1e6  # μm
        focal_y = focal_x.copy()

        # Find centroids via center of mass
        def centroid_x(I):
            total = I.sum()
            if total == 0:
                return 0.0
            col_sum = I.sum(axis=0)
            return float(np.sum(focal_x * col_sum) / total)

        def centroid_y(I):
            total = I.sum()
            if total == 0:
                return 0.0
            row_sum = I.sum(axis=1)
            return float(np.sum(focal_y * row_sum) / total)

        cx_orig = centroid_x(I_orig)
        cy_orig = centroid_y(I_orig)
        cx_blazed = centroid_x(I_blazed)
        cy_blazed = centroid_y(I_blazed)

        print(f"Original  centroid: x={cx_orig:.2f} μm, y={cy_orig:.2f} μm")
        print(f"Blazed    centroid: x={cx_blazed:.2f} μm, y={cy_blazed:.2f} μm")
        print(f"Measured shift:     Δx={cx_blazed - cx_orig:.2f} μm, "
              f"Δy={cy_blazed - cy_orig:.2f} μm")
        print(f"Requested shift:    Δx={args.shift_x:.2f} μm, "
              f"Δy={args.shift_y:.2f} μm")

        # Also find peak position
        peak_orig = np.unravel_index(np.argmax(I_orig), I_orig.shape)
        peak_blazed = np.unravel_index(np.argmax(I_blazed), I_blazed.shape)
        print(f"Original  peak: x={focal_x[peak_orig[1]]:.2f} μm, "
              f"y={focal_y[peak_orig[0]]:.2f} μm")
        print(f"Blazed    peak: x={focal_x[peak_blazed[1]]:.2f} μm, "
              f"y={focal_y[peak_blazed[0]]:.2f} μm")

        # Save NPZ for quick comparison
        np.savez_compressed(
            str(out_dir / "verification.npz"),
            I_orig=I_orig.astype(np.float32),
            I_blazed=I_blazed.astype(np.float32),
            focal_x_um=focal_x.astype(np.float64),
            focal_y_um=focal_y.astype(np.float64),
            centroid_orig_x=cx_orig,
            centroid_orig_y=cy_orig,
            centroid_blazed_x=cx_blazed,
            centroid_blazed_y=cy_blazed,
            shift_x_requested=args.shift_x,
            shift_y_requested=args.shift_y,
        )
        print(f"Saved: {out_dir / 'verification.npz'}")

    print("Done.")


if __name__ == "__main__":
    main()
