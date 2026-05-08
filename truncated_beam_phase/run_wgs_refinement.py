"""WGS refinement with truncated-beam input field and numerical initial phase.

Loads the numerical initial phase (computed by compute_numerical_phase.py) and
the truncated-beam amplitude at the DOE plane, then runs WGS optimisation.
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
from src.mraf_gs import run_refinement
from src.rtad_target import make_rtad_rect_target
from src.utils import ensure_unique_dir, timestamp, save_json


# ---------------------------------------------------------------------------
# Reuse the same aperture-simulation functions
# ---------------------------------------------------------------------------
def angular_spectrum_propagate(field_in, dx_m, dy_m, wavelength_m, distance_m, xp):
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
    amp = xp.where(transmitted > 0, amp / xp.sqrt(transmitted), amp)
    return amp.astype(xp.complex64, copy=False), float(transmitted) * 100.0 if float(transmitted) > 0 else 0.0


# ---------------------------------------------------------------------------
# Parameters (same as compute_numerical_phase.py)
# ---------------------------------------------------------------------------
PHYSICAL = {"wavelength_m": 532e-9, "focal_length_m": 300e-3}
GRID = {"N": 2048, "focal_dx_um": 2.5, "focal_dy_um": 2.5, "dx_doe_m": 3.1171875e-05}
TARGET_CONFIG = {
    "W50_um": 330.0, "H50_um": 120.0,
    "delta_x_um": 15.0, "delta_y_um": 8.0,
    "guard_x_um": 20.0, "guard_y_um": 12.0,
    "constraint_mode": "truncated_rtad",
    "release_level": float(np.exp(-2.0)),
    "center_x_um": 0.0, "center_y_um": 0.0,
}
D_EXP_MM = 25.0
Z_MM = 100.0
APERTURE_DIAMETER_MM = 5.0
CLEAR_APERTURE_MM = 15.0

# WGS parameters (matching baseline)
WGS_PARAMS = {
    "method": "wgs",
    "num_iters": 200,
    "mraf_factor": 0.4,
    "wgs_strategy": "flat_local",
    "wgs_iters": 200,
    "wgs_feedback_exponent": 0.8,
    "wgs_weight_min": 0.5,
    "wgs_weight_max": 2.0,
    "bg_mode": "attenuate",
    "bg_factor": 0.9,
    "metrics_interval": 10,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="WGS refinement with truncated beam")
    parser.add_argument("--cpu", action="store_true")
    parser.add_argument("--phase-npy", type=str,
                        default=str(THIS_DIR / "artifacts" / "20260508-184654_numerical_phase" / "phase_numerical.npy"))
    parser.add_argument("--amp-npy", type=str,
                        default=str(THIS_DIR / "artifacts" / "20260508-184654_numerical_phase" / "amplitude_at_doe.npy"))
    parser.add_argument("--field-npy", type=str,
                        default=str(THIS_DIR / "artifacts" / "20260508-184654_numerical_phase" / "field_at_doe.npy"))
    parser.add_argument("--iters", type=int, default=200)
    args = parser.parse_args()

    backend = get_backend(use_cupy=not args.cpu, verbose=True)
    xp = backend.xp
    dtype = backend.float_dtype

    N = int(GRID["N"])
    dx_doe_m = float(GRID["dx_doe_m"])
    shape = (N, N)
    wavelength_m = float(PHYSICAL["wavelength_m"])

    # ------------------------------------------------------------------
    # 1. Load numerical phase
    # ------------------------------------------------------------------
    phase_path = Path(args.phase_npy)
    if phase_path.exists():
        print(f"Loading numerical phase from {phase_path}")
        phase_numerical_np = np.load(phase_path)
    else:
        print(f"Phase file not found: {phase_path}")
        print("Run compute_numerical_phase.py first.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 2. Build truncated-beam amplitude at DOE plane
    # ------------------------------------------------------------------
    field_path = Path(args.field_npy)
    if field_path.exists():
        print(f"Loading DOE-plane field from {field_path}")
        field_doe_np = np.load(field_path)
        amp_doe_np = np.abs(field_doe_np).astype(np.float32)
        # Extract diffraction phase
        phase_diffraction_np = np.angle(field_doe_np).astype(np.float32)
    else:
        print("Computing truncated beam field...")
        field_ap, throughput = make_truncated_gaussian_at_aperture(
            shape, dx_doe_m, D_EXP_MM, APERTURE_DIAMETER_MM, xp, dtype
        )
        field_doe = angular_spectrum_propagate(
            field_ap, dx_doe_m, dx_doe_m, wavelength_m, Z_MM * 1e-3, xp
        )
        # 15mm clear aperture
        x_doe = (xp.arange(N, dtype=dtype) - (N // 2)) * dtype(dx_doe_m)
        y_doe = (xp.arange(N, dtype=dtype) - (N // 2)) * dtype(dx_doe_m)
        Xd, Yd = xp.meshgrid(x_doe, y_doe)
        r2 = Xd * Xd + Yd * Yd
        r_cl = dtype(CLEAR_APERTURE_MM * 1e-3 / 2.0)
        field_doe = xp.where(r2 <= r_cl * r_cl, field_doe, dtype(0) + 0j * dtype(0))
        field_doe = field_doe.astype(xp.complex64)
        p = dtype(xp.sum(xp.abs(field_doe) ** 2))
        if float(p) > 0:
            field_doe = field_doe / xp.sqrt(p)

        if hasattr(field_doe, "get"):
            field_doe_np = field_doe.get()
        else:
            field_doe_np = np.asarray(field_doe)
        amp_doe_np = np.abs(field_doe_np).astype(np.float32)
        phase_diffraction_np = np.angle(field_doe_np).astype(np.float32)

    # ------------------------------------------------------------------
    # 3. Prepare input for WGS
    # ------------------------------------------------------------------
    # Amplitude: real truncated-beam amplitude (normalised to unit power)
    input_amp = amp_doe_np

    # Phase: numerical phase + diffraction phase, wrapped to [0, 2π)
    phase0_total = np.mod(phase_numerical_np.astype(np.float64) + phase_diffraction_np.astype(np.float64),
                          np.float64(2.0 * np.pi)).astype(np.float32)

    print(f"  input_amp: shape={input_amp.shape}, max={input_amp.max():.6f}")
    print(f"  phase0: min={phase0_total.min():.4f}, max={phase0_total.max():.4f}")

    # ------------------------------------------------------------------
    # 4. Build RTAD target
    # ------------------------------------------------------------------
    x_um = (np.arange(N, dtype=np.float64) - N // 2) * float(GRID["focal_dx_um"])
    y_um = (np.arange(N, dtype=np.float64) - N // 2) * float(GRID["focal_dy_um"])

    target = make_rtad_rect_target(
        shape=shape, x_um=x_um, y_um=y_um,
        W50_um=TARGET_CONFIG["W50_um"],
        H50_um=TARGET_CONFIG["H50_um"],
        delta_x_um=TARGET_CONFIG["delta_x_um"],
        delta_y_um=TARGET_CONFIG["delta_y_um"],
        guard_x_um=TARGET_CONFIG["guard_x_um"],
        guard_y_um=TARGET_CONFIG["guard_y_um"],
        constraint_mode=TARGET_CONFIG["constraint_mode"],
        release_level=TARGET_CONFIG["release_level"],
        mode="separable",
    )

    # ------------------------------------------------------------------
    # 5. Run WGS refinement
    # ------------------------------------------------------------------
    print(f"\n=== WGS refinement: {args.iters} iterations ===")

    result = run_refinement(
        phase0=phase0_total,
        input_amp=input_amp,
        target_amp=target.A_signal,
        masks=target.masks(),
        x_um=x_um,
        y_um=y_um,
        backend=backend,
        method="wgs",
        num_iters=args.iters,
        wgs_strategy="flat_local",
        wgs_iters=args.iters,
        mraf_factor=0.4,
        wgs_feedback_exponent=0.8,
        wgs_weight_min=0.5,
        wgs_weight_max=2.0,
        bg_mode="attenuate",
        bg_factor=0.9,
        metrics_interval=10,
        show_progress=True,
    )

    # ------------------------------------------------------------------
    # 6. Extract physical DOE phase
    # ------------------------------------------------------------------
    # WGS optimised the total phase φ_total = φ_doe + φ_diffraction
    # Physical DOE phase: φ_doe = φ_total - φ_diffraction
    phase_refined_np = np.asarray(result.phase_refined, dtype=np.float64)
    phase_diff_np = np.asarray(phase_diffraction_np, dtype=np.float64)
    phase_doe_np = np.mod(phase_refined_np - phase_diff_np + np.pi, 2.0 * np.pi) - np.pi
    phase_doe_wrapped = np.mod(phase_doe_np, np.float64(2.0 * np.pi)).astype(np.float32)

    # ------------------------------------------------------------------
    # 7. Save outputs
    # ------------------------------------------------------------------
    stamp = timestamp()
    outdir = ensure_unique_dir(THIS_DIR / "artifacts" / f"{stamp}_wgs_refined")
    print(f"\nOutput: {outdir}")

    np.save(outdir / "phase_refined.npy", phase_doe_wrapped)
    np.save(outdir / "phase_refined_unwrapped.npy", phase_doe_np.astype(np.float32))
    np.save(outdir / "reconstruction_refined.npy",
            np.asarray(result.reconstruction_intensity, dtype=np.float32))
    np.save(outdir / "input_amplitude.npy", input_amp)

    try:
        from scipy.io import savemat
        savemat(outdir / "phase_refined.mat",
                {"phase_refined_rad": phase_doe_wrapped},
                do_compression=True)
    except Exception:
        pass

    # Metrics
    fm = result.final_metrics
    rms_pct = 100.0 * float(fm.get("flat_rms", np.nan)) if np.isfinite(float(fm.get("flat_rms", np.nan))) else np.nan
    print(f"\nFinal metrics:")
    print(f"  flat_uniformity: {fm.get('flat_uniformity', 'N/A')}")
    print(f"  flat_rms: {fm.get('flat_rms', 'N/A')}")
    print(f"  RMS nonuniformity: {rms_pct:.4f}%" if np.isfinite(rms_pct) else "  RMS nonuniformity: N/A")
    print(f"  size50_x: {fm.get('size50_x_um', 'N/A')} um")
    print(f"  size50_y: {fm.get('size50_y_um', 'N/A')} um")
    print(f"  efficiency (e^-2 support): {fm.get('efficiency_support', 'N/A')}")

    config = {
        "physical": PHYSICAL, "grid": GRID, "target": TARGET_CONFIG,
        "D_exp_mm": D_EXP_MM, "z_mm": Z_MM,
        "aperture_diameter_mm": APERTURE_DIAMETER_MM,
        "clear_aperture_mm": CLEAR_APERTURE_MM,
        "wgs_params": WGS_PARAMS,
        "rms_nonuniformity_pct": rms_pct,
    }
    save_json(outdir / "config.json", config)

    print("Done.")


if __name__ == "__main__":
    main()
