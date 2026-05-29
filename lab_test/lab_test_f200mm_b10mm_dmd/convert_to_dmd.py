"""Convert computational DOE phase (2048×2048) to DMD superpixel format (192×192).

DMD parameters:
  - Pixel pitch: 13.7 µm
  - Native resolution: 1024×768
  - Square active area: 768×768
  - Superpixel: 4×4 mirrors → effective resolution 192×192
  - Superpixel pitch: 54.8 µm
  - Physical active area: 10.52 × 10.52 mm

Usage:
  python convert_to_dmd.py <phase.npy> --out <output_dir>
  python convert_to_dmd.py <phase.mat> --phase-var phase_refined --out <output_dir>
  python convert_to_dmd.py <phase.npy> --f 200e-3 --beam 5 --out <output_dir>
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import zoom
from scipy.io import savemat


# DMD constants
DMD_NATIVE_W = 1024
DMD_NATIVE_H = 768
DMD_PIXEL_UM = 13.7
DMD_SQUARE = 768  # square region
DMD_SUPER = 4      # 4×4 mirrors per superpixel
DMD_RES = DMD_SQUARE // DMD_SUPER  # 192
DMD_PITCH_UM = DMD_PIXEL_UM * DMD_SUPER  # 54.8 µm
DMD_PHYS_MM = DMD_RES * DMD_PITCH_UM * 1e-3  # 10.52 mm


def compute_dx_doe(lambda_m: float, f_m: float, N: int, focal_dx_um: float) -> float:
    doe_extent = lambda_m * f_m / (focal_dx_um * 1e-6)
    return doe_extent / N


def load_phase(path: str, var: str | None = None) -> tuple[np.ndarray, str]:
    """Load phase from .npy or .mat. Returns (phase, label)."""
    p = Path(path)
    if p.suffix == '.npy':
        phase = np.load(path).astype(np.float64)
        return phase, "npy"
    elif p.suffix == '.mat':
        try:
            import h5py
            with h5py.File(path, 'r') as f:
                if var and var in f:
                    phase = np.array(f[var]).T
                elif var is None:
                    # find first phase-like dataset
                    for k in f.keys():
                        if 'phase' in k.lower() and hasattr(f[k], 'shape'):
                            phase = np.array(f[k]).T
                            break
                    else:
                        raise KeyError(f"No phase variable found in {path}")
                else:
                    raise KeyError(f"Variable '{var}' not found in {path}")
            return phase, "h5py"
        except Exception:
            from scipy.io import loadmat
            mat = loadmat(path)
            if var and var in mat:
                phase = mat[var]
            elif var is None:
                for k in mat:
                    if 'phase' in k.lower():
                        phase = mat[k]
                        break
                else:
                    raise KeyError(f"No phase variable found in {path}")
            else:
                raise KeyError(f"Variable '{var}' not found in {path}")
            return phase.astype(np.float64), "loadmat"
    else:
        raise ValueError(f"Unsupported file format: {p.suffix}")


def convert_to_dmd(
    phase: np.ndarray,
    dx_doe_um: float | None = None,
    N: int = 2048,
    lambda_m: float = 532e-9,
    f_m: float = 200e-3,
    focal_dx_um: float = 2.5,
    dmd_res: int = DMD_RES,
    dmd_pitch_um: float = DMD_PITCH_UM,
) -> np.ndarray:
    """Convert computational phase to DMD superpixel format.

    Parameters
    ----------
    phase : (N, N) float64
        Computational DOE phase in radians [0, 2π).
    dx_doe_um : float or None
        Computational pixel pitch in µm. Auto-computed if None.
    N : int
        Computational grid size (default 2048).
    lambda_m, f_m, focal_dx_um : float
        Used to compute dx_doe_um if not explicitly provided.
    dmd_res : int
        DMD superpixel resolution (default 192).
    dmd_pitch_um : float
        DMD superpixel pitch in µm (default 54.8).

    Returns
    -------
    phase_dmd : (dmd_res, dmd_res) float64
        DMD phase in radians [0, 2π).
    """
    if dx_doe_um is None:
        dx_doe_um = compute_dx_doe(lambda_m, f_m, N, focal_dx_um) * 1e6

    phys_mm = dmd_res * dmd_pitch_um * 1e-3
    crop_px = int(phys_mm / (dx_doe_um * 1e-3))
    cx, cy = N // 2, N // 2

    x0 = cx - crop_px // 2
    y0 = cy - crop_px // 2
    phase_crop = phase[y0:y0 + crop_px, x0:x0 + crop_px]

    zoom_ratio = dmd_res / crop_px
    print(f"dx_doe = {dx_doe_um:.4f} um")
    print(f"DMD physical: {phys_mm:.2f} × {phys_mm:.2f} mm")
    print(f"Crop: {crop_px}×{crop_px} px → {dmd_res}×{dmd_res} (zoom {zoom_ratio:.4f}x)")

    # Complex-field interpolation
    cfield = np.exp(1j * phase_crop)
    real_z = zoom(cfield.real, zoom_ratio, order=3)
    imag_z = zoom(cfield.imag, zoom_ratio, order=3)
    phase_dmd = np.arctan2(imag_z, real_z)
    phase_dmd = np.mod(phase_dmd, 2.0 * np.pi)

    return phase_dmd


def save_outputs(phase_dmd: np.ndarray, out_dir: str, label: str = "phase_dmd") -> None:
    """Save DMD phase as .npy, .mat, and 8-bit .png."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # .npy
    npy_path = out / f"{label}.npy"
    np.save(npy_path, phase_dmd.astype(np.float32))
    print(f"Saved: {npy_path}")

    # .mat
    mat_path = out / f"{label}.mat"
    savemat(mat_path, {"phase_dmd_rad": phase_dmd.astype(np.float64)}, do_compression=True)
    print(f"Saved: {mat_path}")

    # .png (8-bit grayscale)
    from PIL import Image
    gray = (phase_dmd / (2.0 * np.pi) * 255).clip(0, 255).astype(np.uint8)
    png_path = out / f"{label}.png"
    Image.fromarray(gray, mode="L").save(png_path)
    print(f"Saved: {png_path}  ({gray.min()}-{gray.max()}, mean={gray.mean():.1f})")


def main():
    parser = argparse.ArgumentParser(
        description="Convert computational DOE phase to DMD superpixel format."
    )
    parser.add_argument("phase", help="Path to phase .npy or .mat (2048×2048)")
    parser.add_argument("--phase-var", default=None, help="Variable name if .mat file")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--label", default="phase_dmd_192x192", help="Output filename prefix")
    parser.add_argument("--f", type=float, default=200e-3, help="Focal length in m")
    parser.add_argument("--wavelength", type=float, default=532e-9, help="Wavelength in m")
    parser.add_argument("--focal-dx", type=float, default=2.5, help="Focal-plane sampling in µm")
    parser.add_argument("--dx-doe", type=float, default=None, help="Computational dx_doe in µm (auto)")
    parser.add_argument("--dmd-res", type=int, default=192, help="DMD superpixel resolution")
    parser.add_argument("--dmd-pitch", type=float, default=54.8, help="DMD superpixel pitch in µm")
    args = parser.parse_args()

    phase, loader = load_phase(args.phase, args.phase_var)
    print(f"Loaded phase: {phase.shape}, [{phase[np.isfinite(phase)].min():.4f}, "
          f"{phase[np.isfinite(phase)].max():.4f}] ({loader})")

    # Handle NaN outside aperture
    nan_count = np.sum(~np.isfinite(phase))
    phase = np.nan_to_num(phase, nan=0.0)
    if nan_count:
        print(f"Replaced {nan_count} NaN values with 0")

    phase_dmd = convert_to_dmd(
        phase,
        dx_doe_um=args.dx_doe,
        N=phase.shape[0],
        lambda_m=args.wavelength,
        f_m=args.f,
        focal_dx_um=args.focal_dx,
        dmd_res=args.dmd_res,
        dmd_pitch_um=args.dmd_pitch,
    )

    print(f"DMD phase: {phase_dmd.shape}, [{phase_dmd.min():.4f}, {phase_dmd.max():.4f}]")
    save_outputs(phase_dmd, args.out, args.label)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
