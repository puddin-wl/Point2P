"""Phase-shift sweep → SLM BMP generation.

Simulates SLM misalignment compensation: shifts the DOE phase pattern
relative to the fixed Gaussian beam, adds 200μm X+Y blazed grating,
then converts each to SLM 1024×1024 BMP.

Naming convention: shift_Y{+/-}NN_X{+/-}NN

Usage:
  python shift_sweep_slm.py <phase.npy> --out <output_root>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

# ── reuse existing SLM conversion ──────────────────────────────────────────
from convert_to_slm import convert_to_slm, save_outputs


def load_phase(path: str) -> tuple[np.ndarray, str]:
    """Load phase from .npy or .mat."""
    p = Path(path)
    if p.suffix == ".npy":
        return np.load(path).astype(np.float64), "npy"
    if p.suffix == ".mat":
        try:
            import h5py
            with h5py.File(path, "r") as f:
                for k in f.keys():
                    if "phase" in k.lower() and hasattr(f[k], "shape"):
                        return np.array(f[k]).T.astype(np.float64), "h5py"
        except Exception:
            pass
        from scipy.io import loadmat
        mat = loadmat(path)
        for k in mat:
            if "phase" in k.lower():
                return mat[k].astype(np.float64), "loadmat"
        raise KeyError(f"No phase variable found in {path}")
    raise ValueError(f"Unsupported format: {p.suffix}")


def shift_phase(
    phase: np.ndarray,
    shift_y_px: int,
    shift_x_px: int,
) -> np.ndarray:
    """Roll the phase array by integer pixel shifts.

    Positive shift_y = move phase DOWN (beam effectively moves UP relative
    to phase centre), which compensates for "top bright, bottom dark" when
    the SLM is physically too low.
    """
    return np.roll(np.roll(phase, shift_y_px, axis=0), shift_x_px, axis=1)


def add_blaze(
    phase: np.ndarray,
    dx_doe_m: float,
    wavelength_m: float,
    focal_length_m: float,
    shift_x_um: float,
    shift_y_um: float,
) -> np.ndarray:
    """Add blazed grating for focal-spot shift."""
    Ny, Nx = phase.shape
    x_m = (np.arange(Nx, dtype=np.float64) - Nx // 2) * dx_doe_m
    y_m = (np.arange(Ny, dtype=np.float64) - Ny // 2) * dx_doe_m
    X_m, Y_m = np.meshgrid(x_m, y_m)
    k = 2.0 * np.pi / (wavelength_m * focal_length_m)
    blaze = k * (shift_x_um * 1e-6 * X_m + shift_y_um * 1e-6 * Y_m)
    return np.mod(phase + blaze.astype(phase.dtype), 2.0 * np.pi)


def parse_args():
    p = argparse.ArgumentParser(
        description="Phase shift sweep → SLM BMP (with blazed grating)"
    )
    p.add_argument("phase", help="Refined phase .npy path")
    p.add_argument("--out", required=True, help="Output root directory")
    p.add_argument("--f", type=float, default=0.429, help="Focal length (m)")
    p.add_argument("--wavelength", type=float, default=532e-9,
                   help="Wavelength (m)")
    p.add_argument("--focal-dx", type=float, default=2.5,
                   help="Focal-plane sampling (μm)")
    p.add_argument("--N", type=int, default=2048, help="Grid size")
    p.add_argument("--blaze-x", type=float, default=200.0,
                   help="Blazed grating X shift (μm)")
    p.add_argument("--blaze-y", type=float, default=200.0,
                   help="Blazed grating Y shift (μm)")
    p.add_argument("--shifts-y", type=str, default="-20,-15,-10,-5,0,5,10,15,20",
                   help="Comma-separated Y pixel shifts")
    p.add_argument("--shifts-x", type=str, default="0",
                   help="Comma-separated X pixel shifts (applied in a grid with Y)")
    p.add_argument("--slm-res", type=int, default=1024)
    p.add_argument("--slm-pitch", type=float, default=17.0)
    return p.parse_args()


def main():
    args = parse_args()

    phase, loader = load_phase(args.phase)
    print(f"Loaded phase: {phase.shape} [{phase.min():.4f}, {phase.max():.4f}] "
          f"({loader})")

    # physical params
    lam = args.wavelength
    f_m = args.f
    focal_dx_m = args.focal_dx * 1e-6
    doe_extent = lam * f_m / focal_dx_m
    dx_doe = doe_extent / args.N
    print(f"dx_doe = {dx_doe*1e6:.4f} μm  "
          f"(1 px shift = {dx_doe*1e6:.1f} μm on DOE)")

    shifts_y = [int(s.strip()) for s in args.shifts_y.split(",")]
    shifts_x = [int(s.strip()) for s in args.shifts_x.split(",")]

    # aperture info for sanity checks
    aperture_diam_m = 15e-3
    aperture_px = aperture_diam_m / dx_doe
    slm_phys = args.slm_res * args.slm_pitch * 1e-6
    slm_px = slm_phys / dx_doe
    margin_px = (slm_px - aperture_px) / 2.0
    print(f"Aperture: {aperture_px:.1f} px, SLM window: {slm_px:.1f} px, "
          f"margin: {margin_px:.1f} px each side")
    max_abs_shift = max(abs(s) for s in shifts_y + shifts_x)
    if max_abs_shift > margin_px:
        print(f"⚠  WARNING: max shift |{max_abs_shift}| px exceeds margin "
              f"({margin_px:.1f} px) — phase edge may enter aperture!")

    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)

    summary: list[dict] = []

    for sy in shifts_y:
        for sx in shifts_x:
            label = f"shift_Y{sy:+d}_X{sx:+d}"
            print(f"\n{'─'*50}\n  {label}\n{'─'*50}")

            # 1) shift phase
            ph_shifted = shift_phase(phase, sy, sx)

            # 2) add blazed grating
            ph_blazed = add_blaze(ph_shifted, dx_doe, lam, f_m,
                                  args.blaze_x, args.blaze_y)

            # 3) convert to SLM
            ph_slm = convert_to_slm(
                ph_blazed,
                dx_doe_um=dx_doe * 1e6,
                N=args.N,
                lambda_m=lam,
                f_m=f_m,
                focal_dx_um=args.focal_dx,
                slm_res=args.slm_res,
                slm_pitch_um=args.slm_pitch,
            )

            # 4) save
            out_dir = root / label / "slm"
            save_outputs(ph_slm, str(out_dir),
                         label=f"phase_slm_{args.slm_res}x{args.slm_res}")

            # also save the full 2048 blazed phase for record
            full_dir = root / label
            np.save(str(full_dir / "phase_blazed_2048.npy"),
                    ph_blazed.astype(np.float32))
            savemat(
                str(full_dir / "phase_blazed_2048.mat"),
                {"phase_blazed_rad": ph_blazed.astype(np.float64)},
                do_compression=True,
            )

            summary.append({
                "label": label,
                "shift_y_px": sy,
                "shift_x_px": sx,
                "shift_y_um": round(sy * dx_doe * 1e6, 2),
                "shift_x_um": round(sx * dx_doe * 1e6, 2),
                "blaze_x_um": args.blaze_x,
                "blaze_y_um": args.blaze_y,
                "slm_bmp": str(out_dir / f"phase_slm_{args.slm_res}x{args.slm_res}.bmp"),
            })

    # ── summary table ──
    print(f"\n{'='*70}")
    print(f"{'Label':<22s} {'shift_y':>8s} {'shift_x':>8s}  "
          f"  {'Y μm':>8s} {'X μm':>8s}  BMP")
    print(f"{'─'*22} {'─'*8} {'─'*8}  {'─'*8} {'─'*8}  {'─'*40}")
    for s in summary:
        bmp_ok = Path(s["slm_bmp"]).exists()
        print(f"{s['label']:<22s} {s['shift_y_px']:>+4d} px "
              f"{s['shift_x_px']:>+4d} px  "
              f"{s['shift_y_um']:>+8.1f} {s['shift_x_um']:>+8.1f}  "
              f"{'✓' if bmp_ok else '✗'} {s['slm_bmp']}")
    print(f"{'='*70}")
    print(f"Done. {len(summary)} SLM BMPs generated → {root}")

    # quick reference: which direction is which
    print(f"""
┌─────────────────────────────────────────┐
│  Shift convention (for the experimenter) │
├─────────────────────────────────────────┤
│  +Y shift → phase moves DOWN on SLM     │
│          → fixes "TOP bright" problem   │
│  -Y shift → phase moves UP on SLM       │
│          → fixes "BOTTOM bright" problem│
│                                         │
│  +X shift → phase moves RIGHT on SLM    │
│  -X shift → phase moves LEFT on SLM     │
│                                         │
│  A positive blaze shifts focal spot     │
│  toward +X / +Y in the focal plane.     │
│                                         │
│  1 px shift = {dx_doe*1e3:.1f} μm on DOE │
└─────────────────────────────────────────┘
""")


if __name__ == "__main__":
    main()
