"""Verify phase-shift effect on flat-top uniformity via FFT simulation.

Models the real experiment: Gaussian beam is FIXED at center,
phase pattern is SHIFTED (rolled), then propagated to focal plane.
Reports Y-direction top/bottom asymmetry vs shift amount.

Usage:
  python verify_shift_sweep.py <phase.npy> --shifts-y -20,-15,-10,-5,0,5,10,15,20
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from src.propagation import forward_fft, intensity, make_input_gaussian
from src.rtad_target import make_axis_um, make_rtad_rect_target


def centroid_y(I: np.ndarray, y_um: np.ndarray) -> float:
    """Y centroid (center of mass) in μm."""
    total = I.sum()
    if total == 0:
        return 0.0
    row_sum = I.sum(axis=1)
    return float(np.sum(y_um * row_sum) / total)


def top_bottom_ratio(I: np.ndarray, y_um: np.ndarray, mask_flat: np.ndarray,
                     center_y: float) -> dict:
    """Compute top/bottom energy asymmetry in the flat-core region."""
    cy_idx = int(np.argmin(np.abs(y_um - center_y)))
    I_flat = I * mask_flat

    top_mask = np.zeros_like(mask_flat)
    top_mask[cy_idx:, :] = mask_flat[cy_idx:, :]
    bot_mask = np.zeros_like(mask_flat)
    bot_mask[:cy_idx, :] = mask_flat[:cy_idx, :]

    top_sum = float(I_flat[top_mask].sum())
    bot_sum = float(I_flat[bot_mask].sum())

    if bot_sum > 0:
        ratio = top_sum / bot_sum
    else:
        ratio = float("nan")

    return {
        "top_energy": top_sum,
        "bot_energy": bot_sum,
        "top_bot_ratio": ratio,
        "asymmetry_pct": abs(ratio - 1.0) * 100.0 if np.isfinite(ratio) else float("nan"),
    }


def y_profile_asymmetry(I_norm: np.ndarray, x_um: np.ndarray, y_um: np.ndarray,
                        center_x: float, center_y: float) -> dict:
    """Measure Y-profile top/bottom asymmetry around center."""
    ix0 = int(np.argmin(np.abs(x_um - center_x)))
    prof = I_norm[:, ix0].astype(np.float64)
    cy_idx = int(np.argmin(np.abs(y_um - center_y)))

    # top half and bottom half (relative to center)
    top_half = prof[cy_idx:]  # positive y
    bot_half = prof[:cy_idx+1][::-1]  # negative y, flipped

    min_len = min(len(top_half), len(bot_half))
    top_half = top_half[:min_len]
    bot_half = bot_half[:min_len]

    diff = top_half - bot_half
    rms_diff = float(np.sqrt(np.mean(diff**2))) if min_len > 0 else np.nan

    # mean values
    top_mean = float(np.mean(top_half)) if min_len > 0 else np.nan
    bot_mean = float(np.mean(bot_half)) if min_len > 0 else np.nan

    return {
        "top_mean": top_mean,
        "bot_mean": bot_mean,
        "rms_difference": rms_diff,
        "mean_ratio": top_mean / bot_mean if bot_mean and bot_mean > 0 else np.nan,
    }


def verify_one_shift(
    phase: np.ndarray,
    shift_y: int,
    shift_x: int,
    dx_doe_m: float,
    lam: float,
    f_m: float,
    beam_diam_m: float,
    target_cfg: dict,
    focal_dx_um: float,
    focal_dy_um: float,
    amp: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
) -> dict:
    """Propagate shifted phase and compute asymmetry metrics."""
    # Shift phase
    ph = np.roll(np.roll(phase, shift_y, axis=0), shift_x, axis=1)

    # Propagate (Gaussian stays centered — beam is fixed, SLM is shifted)
    field = amp * np.exp(1j * ph).astype(np.complex64)
    I_raw = intensity(forward_fft(field, np), np).astype(np.float64)

    # Build target for mask (always centered — focal plane reference)
    target = make_rtad_rect_target(
        shape=phase.shape, x_um=x_um, y_um=y_um, **target_cfg,
    )
    mask_flat = target.mask_flat

    # Flat-core normalization
    flat_mean = float(np.mean(I_raw[mask_flat]))
    I_norm = (I_raw / flat_mean).astype(np.float64) if flat_mean > 0 else I_raw

    # Centroid
    cy = centroid_y(I_raw, y_um)

    # Top/bottom ratio
    tb = top_bottom_ratio(I_raw, y_um, mask_flat, 0.0)

    # Y profile asymmetry
    ya = y_profile_asymmetry(I_norm, x_um, y_um, 0.0, 0.0)

    # Flat uniformity
    flat_vals = I_norm[mask_flat]
    rms = float(np.std(flat_vals) / np.mean(flat_vals)) * 100.0 if flat_vals.size else np.nan

    return {
        "shift_y_px": shift_y,
        "shift_x_px": shift_x,
        "shift_y_um": round(shift_y * dx_doe_m * 1e6, 2),
        "shift_x_um": round(shift_x * dx_doe_m * 1e6, 2),
        "centroid_y_um": round(cy, 3),
        "rms_nonuniformity_pct": round(rms, 4),
        "top_bot_ratio": round(tb["top_bot_ratio"], 6),
        "asymmetry_pct": round(tb["asymmetry_pct"], 4),
        "y_profile_rms_diff": round(ya["rms_difference"], 6),
        "y_profile_mean_ratio": round(ya["mean_ratio"], 6),
    }


def parse_args():
    p = argparse.ArgumentParser(
        description="Verify phase-shift effect on uniformity"
    )
    p.add_argument("phase", help="Refined phase .npy path")
    p.add_argument("--shifts-y", type=str,
                   default="-20,-15,-10,-5,0,5,10,15,20")
    p.add_argument("--shifts-x", type=str, default="0")
    p.add_argument("--beam", type=float, default=6.5,
                   help="Beam 1/e^2 diameter (mm)")
    p.add_argument("--f", type=float, default=0.429)
    p.add_argument("--wavelength", type=float, default=532e-9)
    p.add_argument("--focal-dx", type=float, default=2.5)
    p.add_argument("--out", default=None)
    return p.parse_args()


def plot_results(results: list[dict], outdir: Path) -> None:
    """Plot shift vs asymmetry."""
    shifts_um = [r["shift_y_um"] for r in results]
    asym = [r["asymmetry_pct"] for r in results]
    rms = [r["rms_nonuniformity_pct"] for r in results]
    ratio = [r["top_bot_ratio"] for r in results]
    profile_diff = [r["y_profile_rms_diff"] for r in results]
    centroid = [r["centroid_y_um"] for r in results]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    # Top-left: asymmetry %
    ax = axes[0, 0]
    ax.plot(shifts_um, asym, "o-", linewidth=1.5, markersize=6)
    ax.axhline(0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Phase shift Y (μm on DOE)")
    ax.set_ylabel("Top/Bottom asymmetry (%)")
    ax.set_title("Flat-core energy asymmetry")
    ax.grid(True, alpha=0.3)

    # Top-right: RMS uniformity
    ax = axes[0, 1]
    ax.plot(shifts_um, rms, "s-", linewidth=1.5, markersize=6, color="tab:red")
    ax.set_xlabel("Phase shift Y (μm on DOE)")
    ax.set_ylabel("RMS nonuniformity (%)")
    ax.set_title("Flat-core RMS")
    ax.grid(True, alpha=0.3)

    # Bottom-left: Y centroid
    ax = axes[1, 0]
    ax.plot(shifts_um, centroid, "^-", linewidth=1.5, markersize=6, color="tab:green")
    ax.axhline(0, color="0.5", linestyle="--", linewidth=0.8)
    ax.set_xlabel("Phase shift Y (μm on DOE)")
    ax.set_ylabel("Centroid Y (μm)")
    ax.set_title("Focal-plane Y centroid")
    ax.grid(True, alpha=0.3)

    # Bottom-right: Y profile RMS difference
    ax = axes[1, 1]
    ax.plot(shifts_um, profile_diff, "D-", linewidth=1.5, markersize=6,
            color="tab:purple")
    ax.set_xlabel("Phase shift Y (μm on DOE)")
    ax.set_ylabel("Top/Bot profile RMS diff")
    ax.set_title("Y profile top/bottom RMS difference")
    ax.grid(True, alpha=0.3)

    fig.savefig(outdir / "shift_sweep_verification.png", dpi=150)
    plt.close(fig)


def plot_profiles(results: list[dict], phase: np.ndarray,
                  amp: np.ndarray, x_um: np.ndarray, y_um: np.ndarray,
                  target_cfg: dict, outdir: Path) -> None:
    """Plot Y center profiles for a subset of shifts."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

    ix0 = int(np.argmin(np.abs(x_um)))

    # Select a few representative shifts for the profile plot
    key_shifts = [-20, -10, 0, 10, 20]
    colors = plt.cm.RdYlBu(np.linspace(0.1, 0.9, len(key_shifts)))

    for shift_y, color in zip(key_shifts, colors):
        ph = np.roll(np.roll(phase, shift_y, axis=0), 0, axis=1)
        field = amp * np.exp(1j * ph).astype(np.complex64)
        I_raw = intensity(forward_fft(field, np), np).astype(np.float64)

        target = make_rtad_rect_target(
            shape=phase.shape, x_um=x_um, y_um=y_um, **target_cfg,
        )
        mask_flat = target.mask_flat
        flat_mean = float(np.mean(I_raw[mask_flat]))
        I_norm = I_raw / flat_mean if flat_mean > 0 else I_raw

        prof_y = I_norm[:, ix0]
        axes[0].plot(y_um, prof_y, color=color, linewidth=1.0,
                     label=f"shift Y={shift_y:+d} px")
        axes[1].plot(y_um, prof_y, color=color, linewidth=1.0,
                     label=f"shift Y={shift_y:+d} px")

    for ax in axes:
        ax.axhline(1.0, color="0.3", linestyle=":", linewidth=0.8)
        ax.axhline(np.exp(-2), color="0.3", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Y (μm)")
        ax.set_ylabel("I / mean(flat core)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7)

    axes[0].set_title("Y center profiles (full view)")
    axes[1].set_title("Y center profiles (zoomed)")
    axes[1].set_xlim(-100, 100)
    axes[1].set_ylim(0.6, 1.4)

    fig.savefig(outdir / "shift_sweep_profiles.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_args()

    phase = np.load(args.phase).astype(np.float64)
    print(f"Loaded phase: {phase.shape}")

    lam = args.wavelength
    f_m = args.f
    focal_dx_um = args.focal_dx
    beam_m = args.beam * 1e-3
    N = phase.shape[0]

    doe_extent = lam * f_m / (focal_dx_um * 1e-6)
    dx_doe = doe_extent / N

    x_um = make_axis_um(N, focal_dx_um)
    y_um = make_axis_um(N, focal_dx_um)

    # Fixed Gaussian beam (centered — this is the "real" beam)
    amp = make_input_gaussian(
        shape=phase.shape,
        dx_doe_m=dx_doe,
        gaussian_1e2_diameter_m=beam_m,
        clear_aperture_m=15e-3,
        xp=np,
        dtype=np.float32,
    )

    target_cfg = {
        "W50_um": 330.0,
        "H50_um": 120.0,
        "delta_x_um": 15.0,
        "delta_y_um": 8.0,
        "guard_x_um": 20.0,
        "guard_y_um": 12.0,
        "constraint_mode": "truncated_rtad",
        "release_level": 0.1353352832366127,
        "target_mode": "separable",
    }

    shifts_y = [int(s.strip()) for s in args.shifts_y.split(",")]
    shifts_x = [int(s.strip()) for s in args.shifts_x.split(",")]

    results = []
    for sy in shifts_y:
        for sx in shifts_x:
            print(f"  shift Y={sy:+d} X={sx:+d} ...", end=" ")
            r = verify_one_shift(
                phase, sy, sx, dx_doe, lam, f_m, beam_m,
                target_cfg, focal_dx_um, focal_dx_um,
                amp, x_um, y_um,
            )
            results.append(r)
            print(f"asym={r['asymmetry_pct']:.2f}%, RMS={r['rms_nonuniformity_pct']:.2f}%")

    # Output
    outdir = Path(args.out) if args.out else Path("shift_verify_output")
    outdir.mkdir(parents=True, exist_ok=True)

    # Save JSON
    with open(outdir / "shift_sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Print table
    print(f"\n{'Shift Y':>8s} {'Shift X':>8s}  "
          f"{'Asym%':>8s} {'RMS%':>8s} {'T/B ratio':>10s} "
          f"{'Prof diff':>10s} {'CentroidY':>10s}")
    print("-" * 72)
    for r in results:
        print(f"{r['shift_y_um']:>+8.1f} {r['shift_x_um']:>+8.1f}  "
              f"{r['asymmetry_pct']:>8.2f} {r['rms_nonuniformity_pct']:>8.2f} "
              f"{r['top_bot_ratio']:>10.4f} "
              f"{r['y_profile_rms_diff']:>10.6f} "
              f"{r['centroid_y_um']:>10.3f}")

    # Plots
    plot_results(results, outdir)
    plot_profiles(results[:9], phase, amp, x_um, y_um, target_cfg, outdir)
    print(f"\nPlots saved to {outdir}")
    print("Done.")


if __name__ == "__main__":
    main()
