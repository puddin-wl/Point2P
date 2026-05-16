"""Simulate beam expansion + aperture clipping effects on DOE focal-plane performance.

Physical model:
  Expanded Gaussian (D_exp) → 5mm hard aperture → angular spectrum propagation
  (z=100mm) → DOE plane × exp(i*phase) → FFT → focal plane.

Compares against the native 5mm Gaussian baseline for which the DOE was designed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = THIS_DIR.parent

# Both lab_test_f300mm and real_world_simulation have a `src/` package.
# Insert PROJECT_DIR first so the full-qualified import below resolves correctly,
# then THIS_DIR for the bare `src.xxx` imports (lab modules).
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(THIS_DIR))

from src.io_mat import load_phase_mat
from src.propagation import forward_fft, intensity, make_input_gaussian
from src.backend import get_backend
from src.utils import ensure_unique_dir, save_json, timestamp

from real_world_simulation.src.metrics import compute_metrics, generate_target, LEVEL_E2
from real_world_simulation.src.plotting import (
    plot_intensity_image,
    plot_profiles_overlay,
    TREND_METRICS,
)


# ---------------------------------------------------------------------------
# Physical parameters (from lab_test_f300mm config)
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
    "constraint_mode": "truncated_rtad",
    "release_level": float(np.exp(-2.0)),
    "center_x_um": 0.0,
    "center_y_um": 0.0,
}

PHASE_MAT = (
    THIS_DIR
    / "artifacts"
    / "20260508-100844_rtad_mraf_gs_truncI0135"
    / "phase_refined.mat"
)

# ---------------------------------------------------------------------------
# Default sweep
# ---------------------------------------------------------------------------
D_EXP_VALUES_MM = [20.0, 22.0, 24.0, 25.0, 26.0, 28.0, 30.0]
Z_MM = 100.0
APERTURE_DIAMETER_MM = 5.0
CLEAR_APERTURE_MM = 15.0  # secondary clear aperture at DOE plane


# ---------------------------------------------------------------------------
# Angular spectrum propagation
# ---------------------------------------------------------------------------
def angular_spectrum_propagate(field_in, dx_m, dy_m, wavelength_m, distance_m, xp):
    """Propagate a complex field between two parallel planes.

    Uses the angular spectrum transfer function.  Evanescent components are
    filtered out (``arg < 0`` → root = 0).
    """
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


# ---------------------------------------------------------------------------
# Truncated Gaussian at the aperture plane
# ---------------------------------------------------------------------------
def make_truncated_gaussian_at_aperture(shape, dx_m, D_exp_mm, aperture_diameter_mm, xp, dtype):
    """Build a Gaussian field truncated by a hard circular aperture.

    The Gaussian is normalized to unit total power *before* truncation, then
    renormalized to unit power after the aperture.
    """
    Ny, Nx = int(shape[0]), int(shape[1])
    x = (xp.arange(Nx, dtype=dtype) - (Nx // 2)) * dtype(dx_m)
    y = (xp.arange(Ny, dtype=dtype) - (Ny // 2)) * dtype(dx_m)
    X, Y = xp.meshgrid(x, y)
    r2 = X * X + Y * Y

    w = dtype(D_exp_mm * 1e-3 / 2.0)
    amp = xp.exp(-r2 / (w * w)).astype(dtype, copy=False)

    # Normalise to unit power before aperture.
    pre_power = dtype(xp.sum(amp * amp))
    amp = xp.where(pre_power > 0, amp / xp.sqrt(pre_power), amp)

    # Hard circular aperture.
    r_ap = dtype(aperture_diameter_mm * 1e-3 / 2.0)
    mask = r2 <= (r_ap * r_ap)
    amp = xp.where(mask, amp, dtype(0.0)).astype(dtype, copy=False)

    # Throughput and final normalisation.
    transmitted = dtype(xp.sum(amp * amp))
    throughput = 100.0 * float(transmitted) if float(transmitted) > 0 else 0.0
    amp = xp.where(transmitted > 0, amp / xp.sqrt(transmitted), amp)

    return amp.astype(xp.complex64, copy=False), throughput


# ---------------------------------------------------------------------------
# Single case simulation
# ---------------------------------------------------------------------------
def simulate_case(phase_np, D_exp_mm, z_mm, target, xp, dtype):
    """Run a single aperture-clipping simulation and return results."""
    shape = phase_np.shape
    dx_doe_m = float(GRID["dx_doe_m"])
    wavelength_m = float(PHYSICAL["wavelength_m"])

    # 1. Aperture-plane field
    field_ap, throughput = make_truncated_gaussian_at_aperture(
        shape, dx_doe_m, D_exp_mm, APERTURE_DIAMETER_MM, xp, dtype
    )

    # 2. Propagate to DOE plane
    z_m = z_mm * 1e-3
    if z_m > 1e-9:
        field_doe = angular_spectrum_propagate(
            field_ap, dx_doe_m, dx_doe_m, wavelength_m, z_m, xp
        )
    else:
        field_doe = field_ap

    # 3. Apply secondary clear aperture (15mm) at DOE plane — same as baseline
    Ny, Nx = shape
    x_doe = (xp.arange(Nx, dtype=dtype) - (Nx // 2)) * dtype(dx_doe_m)
    y_doe = (xp.arange(Ny, dtype=dtype) - (Ny // 2)) * dtype(dx_doe_m)
    Xd, Yd = xp.meshgrid(x_doe, y_doe)
    r2_doe = Xd * Xd + Yd * Yd
    r_clear = dtype(CLEAR_APERTURE_MM * 1e-3 / 2.0)
    clear_mask = r2_doe <= (r_clear * r_clear)
    field_doe = xp.where(clear_mask, field_doe, dtype(0.0) + 0j * dtype(0.0))
    field_doe = field_doe.astype(xp.complex64)

    # Re-normalize after clear aperture.
    p_doe = dtype(xp.sum(xp.abs(field_doe) ** 2))
    if float(p_doe) > 0:
        field_doe = field_doe / xp.sqrt(p_doe)

    # 4. Apply DOE phase and propagate to focal plane
    phase_dev = xp.asarray(phase_np, dtype=xp.float32)
    doe_field = field_doe * xp.exp(1j * phase_dev)
    focal_field = forward_fft(doe_field, xp)
    I_dev = intensity(focal_field, xp)

    # 5. Move intensity to NumPy for metrics
    if hasattr(I_dev, "get"):
        I_np = I_dev.get()
    else:
        I_np = np.asarray(I_dev)

    # 6. Compute metrics
    scalar, details, warnings = compute_metrics(I_np, target, float(throughput))

    label = f"Dexp={D_exp_mm:.0f}mm_z={z_mm:.0f}mm"
    row = {
        "D_exp_mm": D_exp_mm,
        "z_mm": z_mm,
        "label": label,
        **scalar,
    }

    return {
        "intensity": I_np,
        "metrics": scalar,
        "details": details,
        "warnings": warnings,
        "row": row,
        "label": label,
        "D_exp_mm": D_exp_mm,
        "z_mm": z_mm,
    }


# ---------------------------------------------------------------------------
# Baseline: native 5mm Gaussian at DOE plane
# ---------------------------------------------------------------------------
def simulate_baseline(phase_np, target, xp, dtype):
    """Simulate the idealised native 5mm Gaussian input (no aperture clipping)."""
    shape = phase_np.shape
    dx_doe_m = float(GRID["dx_doe_m"])

    amp = make_input_gaussian(
        shape, dx_doe_m,
        gaussian_1e2_diameter_m=5e-3,
        clear_aperture_m=CLEAR_APERTURE_MM * 1e-3,
        xp=xp, dtype=dtype,
    )

    phase_dev = xp.asarray(phase_np, dtype=xp.float32)
    doe_field = amp.astype(xp.complex64) * xp.exp(1j * phase_dev)
    focal_field = forward_fft(doe_field, xp)
    I_dev = intensity(focal_field, xp)

    if hasattr(I_dev, "get"):
        I_np = I_dev.get()
    else:
        I_np = np.asarray(I_dev)

    scalar, details, warnings = compute_metrics(I_np, target, 100.0)

    return {
        "intensity": I_np,
        "metrics": scalar,
        "details": details,
        "warnings": warnings,
        "row": {
            "D_exp_mm": "baseline",
            "z_mm": "baseline",
            "label": "baseline_native_5mm",
            **scalar,
        },
        "label": "baseline (native 5mm Gaussian)",
        "D_exp_mm": None,
        "z_mm": None,
    }


# ---------------------------------------------------------------------------
# Custom metric-trend figure (D_exp on x-axis)
# ---------------------------------------------------------------------------
METRICS_TO_PLOT = [
    ("rms_nonuniformity_percent", "RMS nonuniformity / %"),
    ("efficiency_e2_percent", "e^-2 efficiency / %"),
    ("aperture_throughput_percent", "aperture throughput / %"),
    ("size50_x_um", "size50 x / um"),
    ("size50_y_um", "size50 y / um"),
    ("transition_13p5_90_x_um", "transition x / um"),
]


def make_Dexp_trends_figure(rows):
    """Metric trends with D_exp on the x-axis."""
    import matplotlib.pyplot as plt

    x = np.asarray([float(r["D_exp_mm"]) for r in rows], dtype=np.float64)
    labels = [f"{v:.0f}" for v in x]

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2), constrained_layout=True)
    for ax, (key, ylabel) in zip(axes.ravel(), METRICS_TO_PLOT):
        y = np.asarray([r.get(key, np.nan) for r in rows], dtype=np.float64)
        ax.plot(x, y, marker="o", linewidth=1.2)
        ax.set_title(ylabel)
        ax.set_xlabel("D_exp / mm")
        ax.grid(True, alpha=0.25)
        if len(labels) <= 9:
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Aperture clipping simulation for f=300mm DOE"
    )
    parser.add_argument(
        "--D-exp", nargs="*", type=float,
        help="Expanded beam 1/e² diameters in mm (default: 20 22 24 25 26 28 30)",
    )
    parser.add_argument(
        "--z-mm", type=float, default=Z_MM,
        help=f"Aperture-to-DOE distance in mm (default: {Z_MM})",
    )
    parser.add_argument("--cpu", action="store_true", help="Force CPU/NumPy backend")
    parser.add_argument("--phase-mat", type=str, default=str(PHASE_MAT))
    args = parser.parse_args()

    D_exp_list = args.D_exp or D_EXP_VALUES_MM
    z_mm = args.z_mm
    phase_path = Path(args.phase_mat)

    # ------------------------------------------------------------------
    # Backend
    # ------------------------------------------------------------------
    backend = get_backend(use_cupy=not args.cpu, verbose=True)
    xp = backend.xp
    dtype = backend.float_dtype

    # ------------------------------------------------------------------
    # Load phase
    # ------------------------------------------------------------------
    print(f"\nLoading phase from {phase_path}")
    phase_np, info = load_phase_mat(phase_path, phase_var="phase_refined", swap_xy=False)
    N = int(phase_np.shape[0])
    print(f"  shape: {phase_np.shape}, dtype: {phase_np.dtype}, "
          f"min={phase_np.min():.4f}, max={phase_np.max():.4f}")
    if not info.get("phase_xy_swapped", False) and N != 2048:
        print(f"  WARNING: unexpected phase shape {phase_np.shape}, expected (2048,2048)")

    # ------------------------------------------------------------------
    # Target
    # ------------------------------------------------------------------
    target = generate_target(
        shape=(N, N),
        focal_dx_um=float(GRID["focal_dx_um"]),
        focal_dy_um=float(GRID["focal_dy_um"]),
        target_config=TARGET_CONFIG,
    )

    # ------------------------------------------------------------------
    # Output directory
    # ------------------------------------------------------------------
    stamp = timestamp()
    outdir = ensure_unique_dir(
        THIS_DIR / "artifacts" / f"{stamp}_aperture_simulation"
    )
    print(f"\nOutput directory: {outdir}")

    config_snapshot = {
        "physical": PHYSICAL,
        "grid": GRID,
        "target": TARGET_CONFIG,
        "phase_mat": str(phase_path),
        "phase_info": info,
        "D_exp_mm_list": D_exp_list,
        "z_mm": z_mm,
        "aperture_diameter_mm": APERTURE_DIAMETER_MM,
        "clear_aperture_mm": CLEAR_APERTURE_MM,
    }
    save_json(outdir / "config.json", config_snapshot)

    # ------------------------------------------------------------------
    # Baseline
    # ------------------------------------------------------------------
    print("\n=== Baseline (native 5mm Gaussian) ===")
    baseline = simulate_baseline(phase_np, target, xp, dtype)
    bm = baseline["metrics"]
    print(f"  RMS={bm['rms_nonuniformity_percent']:.4f}%, "
          f"e2_eff={bm['efficiency_e2_percent']:.2f}%, "
          f"size50_x={bm['size50_x_um']:.1f}um, "
          f"size50_y={bm['size50_y_um']:.1f}um")
    if baseline["warnings"]:
        for w in baseline["warnings"]:
            print(f"  WARNING: {w}")

    # ------------------------------------------------------------------
    # Sweep
    # ------------------------------------------------------------------
    print(f"\n=== Sweep: {len(D_exp_list)} D_exp values, z={z_mm} mm ===\n")
    results = []
    all_rows = []

    for D_exp_mm in D_exp_list:
        result = simulate_case(phase_np, D_exp_mm, z_mm, target, xp, dtype)
        results.append(result)
        all_rows.append(result["row"])
        m = result["metrics"]
        print(
            f"  D_exp={D_exp_mm:5.1f} mm  "
            f"RMS={m['rms_nonuniformity_percent']:.4f}%  "
            f"eff={m['efficiency_e2_percent']:.2f}%  "
            f"T={m['aperture_throughput_percent']:.2f}%  "
            f"s50x={m['size50_x_um']:.1f}um  "
            f"s50y={m['size50_y_um']:.1f}um"
        )
        if result["warnings"]:
            for w in result["warnings"]:
                print(f"         WARNING: {w}")

    # Add baseline row at the end
    all_rows.append(baseline["row"])
    results.append(baseline)

    # ------------------------------------------------------------------
    # Save CSV / JSON
    # ------------------------------------------------------------------
    import csv

    csv_fields = [
        "D_exp_mm", "z_mm", "label",
        "rms_nonuniformity_percent", "efficiency_e2_percent",
        "aperture_throughput_percent",
        "size50_x_um", "size50_y_um",
        "size13p5_x_um", "size13p5_y_um",
        "size90_x_um", "size90_y_um",
        "transition_13p5_90_x_um", "transition_13p5_90_y_um",
        "center_offset_x_um", "center_offset_y_um",
        "I_ref_flat_mean",
    ]
    with open(outdir / "summary.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    save_json(outdir / "summary.json", {
        "rows": all_rows,
        "config": config_snapshot,
    })

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    print("\n=== Generating plots ===")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # --- baseline intensity ---
    print("  baseline intensity ...")
    plot_intensity_image(
        baseline["intensity"], target,
        outdir / "intensity_baseline.png",
        "Baseline: native 5mm Gaussian",
    )

    # --- Per-case intensity ---
    for result in results:
        if result["D_exp_mm"] is None:
            continue
        D = result["D_exp_mm"]
        fname = f"intensity_D{D:.0f}_z{z_mm:.0f}mm.png"
        print(f"  {fname} ...")
        plot_intensity_image(
            result["intensity"], target,
            outdir / fname,
            result["label"],
        )

    # --- Center profiles overlay ---
    print("  center profiles overlay ...")
    sweep_results = [r for r in results if r["D_exp_mm"] is not None]
    profile_data = [baseline] + sweep_results
    plot_profiles_overlay(profile_data, target, outdir / "center_profiles_overlay.png")

    # --- Metric trends ---
    print("  metric trends ...")
    sweep_rows = [r["row"] for r in sweep_results]
    fig = make_Dexp_trends_figure(sweep_rows)
    fig.savefig(outdir / "metric_trends.png", dpi=150)
    plt.close(fig)

    # --- Delta intensity grid: each case minus baseline ---
    print("  delta intensity grid ...")
    I_base = baseline["intensity"]
    I_base_norm, _ = _normalize_by_flat(I_base, target.mask_flat)
    n_cases = len(sweep_results)
    ncols = min(4, n_cases)
    nrows = (n_cases + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3.5 * ncols, 3.2 * nrows),
        constrained_layout=True,
        squeeze=False,
    )
    for idx, result in enumerate(sweep_results):
        ax = axes[idx // ncols][idx % ncols]
        I_norm, _ = _normalize_by_flat(result["intensity"], target.mask_flat)
        delta = I_norm - I_base_norm
        extent = [float(target.x_um[0]), float(target.x_um[-1]),
                  float(target.y_um[0]), float(target.y_um[-1])]
        vmax = max(0.3, float(np.nanpercentile(np.abs(delta), 99.5)))
        im = ax.imshow(delta, extent=extent, origin="lower",
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        fig.colorbar(im, ax=ax)
        a2 = float(target.params.get("a2_um", 200.0))
        b2 = float(target.params.get("b2_um", 100.0))
        ax.set_xlim(-a2 - 65, a2 + 65)
        ax.set_ylim(-b2 - 65, b2 + 65)
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        ax.set_title(result["label"])
    # Hide unused axes.
    for idx in range(n_cases, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)
    fig.savefig(outdir / "delta_intensity_grid.png", dpi=150)
    plt.close(fig)

    print(f"\nDone. Output in: {outdir}")


def _normalize_by_flat(intensity, mask_flat):
    """Local helper to avoid cross-module import issue."""
    values = np.asarray(intensity, dtype=np.float64)[mask_flat]
    mean = float(np.mean(values)) if values.size else 1.0
    if not np.isfinite(mean) or mean <= 0:
        return np.zeros_like(intensity, dtype=np.float32), mean
    return (np.asarray(intensity, dtype=np.float64) / mean).astype(np.float32), mean


if __name__ == "__main__":
    main()
