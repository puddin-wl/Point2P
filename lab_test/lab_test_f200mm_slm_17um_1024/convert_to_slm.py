"""Convert computational DOE phase (2048x2048) to SLM format (1024x1024 @ 17um).

要点:
  1. 根据 dx_doe 和 SLM 物理尺寸, 从 2048×2048 计算网格中心裁切对应区域
  2. 复振幅插值 (cubic spline on exp(iφ)), 避免直接对包裹相位插值产生振铃
  3. 输出 .npy, .mat, .bmp (8-bit 灰度) 三种格式

Usage:
  python convert_to_slm.py <phase.npy> --out <output_dir>
  python convert_to_slm.py <phase.npy> --out <output_dir> --slm-res 1024 --slm-pitch 17.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat
from scipy.ndimage import zoom


def compute_dx_doe(lambda_m: float, f_m: float, N: int, focal_dx_um: float) -> float:
    """根据焦面采样间距反推 DOE 面的像素尺寸 dx = λf / (N × dx_focal)."""
    doe_extent = lambda_m * f_m / (focal_dx_um * 1e-6)
    return doe_extent / N


def load_phase(path: str, var: str | None = None) -> tuple[np.ndarray, str]:
    """从 .npy 或 .mat 加载相位 (自动识别格式). 返回 (phase, loader_name)."""
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


def convert_to_slm(
    phase: np.ndarray,
    dx_doe_um: float | None = None,
    N: int = 2048,
    lambda_m: float = 532e-9,
    f_m: float = 200e-3,
    focal_dx_um: float = 2.5,
    slm_res: int = 1024,
    slm_pitch_um: float = 17.0,
) -> np.ndarray:
    """核心转换: 2048×2048 计算相位 → SLM 1024×1024 (中心裁切 + 复振幅 cubic 插值)."""
    # 计算计算网格的像素尺寸 (dx_doe)
    if dx_doe_um is None:
        dx_doe_um = compute_dx_doe(lambda_m, f_m, N, focal_dx_um) * 1e6

    # SLM 物理尺寸 → 裁切对应的计算网格区域
    phys_mm = slm_res * slm_pitch_um * 1e-3
    crop_px = int(phys_mm / (dx_doe_um * 1e-3))
    cx, cy = N // 2, N // 2

    x0 = cx - crop_px // 2
    y0 = cy - crop_px // 2
    phase_crop = phase[y0:y0 + crop_px, x0:x0 + crop_px]

    zoom_ratio = slm_res / crop_px
    print(f"dx_doe = {dx_doe_um:.4f} um, DOE extent = {N * dx_doe_um * 1e-3:.2f} mm")
    print(f"SLM physical: {phys_mm:.2f} x {phys_mm:.2f} mm")
    print(f"Crop: {crop_px}x{crop_px} px -> {slm_res}x{slm_res} (zoom {zoom_ratio:.4f}x)")

    # 复振幅 cubic 插值: 对 exp(iφ) 的实部/虚部分别做 cubic spline zoom,
    # 再用 arctan2 恢复相位。这样可以完全避免 2π→0 包裹跳变带来的振铃。
    cfield = np.exp(1j * phase_crop)
    real_z = zoom(cfield.real, zoom_ratio, order=3)
    imag_z = zoom(cfield.imag, zoom_ratio, order=3)
    phase_slm = np.arctan2(imag_z, real_z)
    phase_slm = np.mod(phase_slm, 2.0 * np.pi)

    return phase_slm


def save_outputs(phase_slm: np.ndarray, out_dir: str, label: str = "phase_slm_1024x1024") -> None:
    """保存 SLM 相位为 .npy, .mat, .bmp (8-bit 灰度) 三种格式."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # .npy
    npy_path = out / f"{label}.npy"
    np.save(npy_path, phase_slm.astype(np.float32))
    print(f"Saved: {npy_path}")

    # .mat
    mat_path = out / f"{label}.mat"
    savemat(mat_path, {"phase_slm_rad": phase_slm.astype(np.float64)}, do_compression=True)
    print(f"Saved: {mat_path}")

    # .bmp (8-bit grayscale)
    from PIL import Image
    gray = (phase_slm / (2.0 * np.pi) * 255).clip(0, 255).astype(np.uint8)
    bmp_path = out / f"{label}.bmp"
    Image.fromarray(gray, mode="L").save(bmp_path, format="BMP")
    print(f"Saved: {bmp_path}  ({gray.min()}-{gray.max()}, mean={gray.mean():.1f})")


def main():
    """CLI 入口: 加载相位 → 转换 → 保存."""
    parser = argparse.ArgumentParser(
        description="Convert computational DOE phase to SLM format (1024x1024 @ 17um)."
    )
    parser.add_argument("phase", help="Path to phase .npy or .mat (2048x2048)")
    parser.add_argument("--phase-var", default=None, help="Variable name if .mat file")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--label", default="phase_slm_1024x1024", help="Output filename prefix")
    parser.add_argument("--f", type=float, default=200e-3, help="Focal length in m")
    parser.add_argument("--wavelength", type=float, default=532e-9, help="Wavelength in m")
    parser.add_argument("--focal-dx", type=float, default=2.5, help="Focal-plane sampling in um")
    parser.add_argument("--slm-res", type=int, default=1024, help="SLM resolution (square)")
    parser.add_argument("--slm-pitch", type=float, default=17.0, help="SLM pixel pitch in um")
    parser.add_argument("--dx-doe", type=float, default=None, help="Computational dx_doe in um (auto)")
    args = parser.parse_args()

    phase, loader = load_phase(args.phase, args.phase_var)
    print(f"Loaded phase: {phase.shape}, [{phase[np.isfinite(phase)].min():.4f}, "
          f"{phase[np.isfinite(phase)].max():.4f}] ({loader})")

    # Handle NaN outside aperture
    nan_count = np.sum(~np.isfinite(phase))
    phase = np.nan_to_num(phase, nan=0.0)
    if nan_count:
        print(f"Replaced {nan_count} NaN values with 0")

    phase_slm = convert_to_slm(
        phase,
        dx_doe_um=args.dx_doe,
        N=phase.shape[0],
        lambda_m=args.wavelength,
        f_m=args.f,
        focal_dx_um=args.focal_dx,
        slm_res=args.slm_res,
        slm_pitch_um=args.slm_pitch,
    )

    print(f"SLM phase: {phase_slm.shape}, [{phase_slm.min():.4f}, {phase_slm.max():.4f}]")
    save_outputs(phase_slm, args.out, args.label)
    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
