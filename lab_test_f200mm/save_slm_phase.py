"""Convert computational DOE phase to SLM-loadable 8-bit grayscale PNG.

The computational phase is on an N×N grid with pixel pitch dx_doe (typically
~20.78 um for f=200mm). The SLM has 1920×1080 pixels at 6.4 um pitch. This
script crops the central physical region that fits on the SLM and resamples.

Usage:
  python save_slm_phase.py <phase.npy> --out <output.png>
  python save_slm_phase.py <phase.npy> --out <output.png> --slm-w 1920 --slm-h 1080 --slm-pitch 6.4
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import zoom


def compute_dx_doe(lambda_m: float, f_m: float, N: int, focal_dx_um: float) -> float:
    """Compute computational DOE pixel pitch from optical parameters."""
    doe_extent = lambda_m * f_m / (focal_dx_um * 1e-6)
    return doe_extent / N


def save_slm_phase(
    phase_path: str,
    out_path: str,
    slm_w: int = 1920,
    slm_h: int = 1080,
    slm_pitch_um: float = 6.4,
    dx_doe_um: float | None = None,
    lambda_m: float = 532e-9,
    f_m: float = 200e-3,
    focal_dx_um: float = 2.5,
) -> None:
    """Convert a computational phase array to an SLM-loadable 8-bit PNG.

    Parameters
    ----------
    phase_path : str
        Path to the phase_refined.npy file (N×N float32, [0, 2π)).
    out_path : str
        Output PNG path.
    slm_w, slm_h : int
        SLM resolution in pixels (default 1920×1080).
    slm_pitch_um : float
        SLM pixel pitch in um (default 6.4).
    dx_doe_um : float or None
        Computational pixel pitch in um. If None, computed from lambda, f, N, focal_dx.
    lambda_m, f_m, focal_dx_um : float
        Used to compute dx_doe_um if not explicitly provided.
    """
    phase = np.load(phase_path).astype(np.float64)
    N = phase.shape[0]

    if dx_doe_um is None:
        dx_doe_um = compute_dx_doe(lambda_m, f_m, N, focal_dx_um) * 1e6

    # Physical size of SLM active area
    slm_phys_x_mm = slm_w * slm_pitch_um * 1e-3  # mm
    slm_phys_y_mm = slm_h * slm_pitch_um * 1e-3

    # Number of computational pixels that fit in the SLM physical area
    crop_w = int(slm_phys_x_mm / (dx_doe_um * 1e-3))
    crop_h = int(slm_phys_y_mm / (dx_doe_um * 1e-3))

    print(f"Computational grid: {N}×{N}, dx_doe = {dx_doe_um:.4f} um")
    print(f"SLM active area: {slm_phys_x_mm:.3f} × {slm_phys_y_mm:.3f} mm")
    print(f"Cropping central {crop_w}×{crop_h} computational pixels")

    # Crop central region
    cx, cy = N // 2, N // 2
    x0 = cx - crop_w // 2
    x1 = x0 + crop_w
    y0 = cy - crop_h // 2
    y1 = y0 + crop_h
    phase_cropped = phase[y0:y1, x0:x1]

    # Resample to SLM resolution
    zoom_y = slm_h / crop_h
    zoom_x = slm_w / crop_w
    print(f"Resampling {crop_w}×{crop_h} → {slm_w}×{slm_h} (zoom {zoom_x:.4f}x, {zoom_y:.4f}x)")

    phase_slm = zoom(phase_cropped, (zoom_y, zoom_x), order=3)

    # Wrap to [0, 2π) and convert to 8-bit
    phase_slm = np.mod(phase_slm, 2.0 * np.pi)
    gray = (phase_slm / (2.0 * np.pi) * 255.0).clip(0, 255).astype(np.uint8)

    # Save
    from PIL import Image
    img = Image.fromarray(gray, mode="L")
    img.save(out_path)
    print(f"Saved SLM phase: {out_path} ({slm_w}×{slm_h}, 8-bit grayscale)")
    print(f"  Gray range: [{gray.min()}, {gray.max()}], mean = {gray.mean():.1f}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert computational phase to SLM PNG.")
    parser.add_argument("phase", help="Path to phase_refined.npy")
    parser.add_argument("--out", required=True, help="Output PNG path")
    parser.add_argument("--slm-w", type=int, default=1920, help="SLM width in pixels")
    parser.add_argument("--slm-h", type=int, default=1080, help="SLM height in pixels")
    parser.add_argument("--slm-pitch", type=float, default=6.4, help="SLM pixel pitch in um")
    parser.add_argument("--dx-doe", type=float, default=None, help="Computational pixel pitch in um (auto if omitted)")
    parser.add_argument("--wavelength", type=float, default=532e-9, help="Wavelength in m")
    parser.add_argument("--f", type=float, default=200e-3, help="Focal length in m")
    parser.add_argument("--focal-dx", type=float, default=2.5, help="Focal-plane sampling in um")
    args = parser.parse_args()

    save_slm_phase(
        phase_path=args.phase,
        out_path=args.out,
        slm_w=args.slm_w,
        slm_h=args.slm_h,
        slm_pitch_um=args.slm_pitch,
        dx_doe_um=args.dx_doe,
        lambda_m=args.wavelength,
        f_m=args.f,
        focal_dx_um=args.focal_dx,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
