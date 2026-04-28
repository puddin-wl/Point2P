"""Command-line entry point for one RTAD MRAF/GS refinement case."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import savemat

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from config_default import CONFIG
from src.backend import get_backend
from src.io_mat import load_optional_mat_variable, load_phase_mat, save_phase_mat
from src.mraf_gs import run_refinement
from src.plotting import (
    plot_convergence,
    plot_intensity,
    plot_masks,
    plot_phase,
    plot_profiles_compare,
    plot_target,
)
from src.propagation import make_input_gaussian
from src.rtad_target import make_axis_um, make_rtad_rect_target
from src.utils import ensure_unique_dir, format_metrics, save_json, save_metrics_csv, timestamp


def parse_args() -> argparse.Namespace:
    """Parse CLI options."""
    parser = argparse.ArgumentParser(description="Run RTAD rectangular MRAF/GS refinement.")
    parser.add_argument("--phase-mat", default=None, help="Path to phase0 MAT file.")
    parser.add_argument("--phase-var", default=None, help="Phase variable name, e.g. phase0.")
    parser.add_argument("--iters", type=int, default=None, help="Number of refinement iterations.")
    parser.add_argument("--method", default=None, choices=["gs", "mraf", "wgs", "mraf_then_wgs"], help="Refinement method.")
    parser.add_argument("--mraf-factor", type=float, default=None, help="MRAF free-region attenuation factor.")
    parser.add_argument("--delta-x", type=float, default=None, help="RTAD x falling-edge half width in um.")
    parser.add_argument("--delta-y", type=float, default=None, help="RTAD y falling-edge half width in um.")
    parser.add_argument("--outdir", default=None, help="Output directory. If omitted, timestamped artifacts dir is used.")
    parser.add_argument("--use-cupy", dest="use_cupy", action="store_true", default=None, help="Use CuPy if available.")
    parser.add_argument("--no-cupy", dest="use_cupy", action="store_false", help="Force NumPy CPU.")
    parser.add_argument("--transpose-h5", action="store_true", help="Transpose HDF5 MAT arrays after reading.")
    parser.add_argument("--metrics-interval", type=int, default=None, help="Metrics logging interval.")
    parser.add_argument("--wgs-after-iters", type=int, default=None, help="Iteration at which WGS updates begin.")
    parser.add_argument("--feedback-exponent", type=float, default=None, help="WGS feedback exponent.")
    parser.add_argument("--bg-mode", choices=["keep", "attenuate", "zero"], default=None, help="MRAF background handling.")
    parser.add_argument("--bg-factor", type=float, default=None, help="MRAF background attenuation if bg-mode=attenuate.")
    parser.add_argument("--smoke-shape", type=int, default=None, help="Shape for random/zero smoke phase when no MAT is supplied.")
    return parser.parse_args()


def apply_overrides(config: dict, args: argparse.Namespace) -> dict:
    """Apply CLI overrides to a copied config."""
    cfg = copy.deepcopy(config)
    if args.phase_mat is not None:
        cfg["paths"]["phase_mat"] = args.phase_mat
    if args.phase_var is not None:
        cfg["paths"]["phase_var"] = args.phase_var
    if args.transpose_h5:
        cfg["paths"]["transpose_h5"] = True
    if args.iters is not None:
        cfg["refinement"]["num_iters"] = args.iters
    if args.method is not None:
        cfg["refinement"]["method"] = args.method
    if args.mraf_factor is not None:
        cfg["refinement"]["mraf_factor"] = args.mraf_factor
    if args.delta_x is not None:
        cfg["target"]["delta_x_um"] = args.delta_x
    if args.delta_y is not None:
        cfg["target"]["delta_y_um"] = args.delta_y
    if args.metrics_interval is not None:
        cfg["refinement"]["metrics_interval"] = args.metrics_interval
    if args.wgs_after_iters is not None:
        cfg["refinement"]["wgs_after_iters"] = args.wgs_after_iters
    if args.feedback_exponent is not None:
        cfg["refinement"]["feedback_exponent"] = args.feedback_exponent
    if args.bg_mode is not None:
        cfg["refinement"]["bg_mode"] = args.bg_mode
    if args.bg_factor is not None:
        cfg["refinement"]["bg_factor"] = args.bg_factor
    if args.use_cupy is not None:
        cfg["runtime"]["use_cupy"] = bool(args.use_cupy)
    if args.smoke_shape is not None:
        cfg["runtime"]["smoke_shape"] = args.smoke_shape
    return cfg


def load_or_make_phase(cfg: dict) -> tuple[np.ndarray, dict, np.ndarray | None, np.ndarray | None]:
    """Load phase0 from MAT or create a clearly marked smoke-test phase."""
    phase_mat = cfg["paths"]["phase_mat"]
    phase_var = cfg["paths"]["phase_var"]
    transpose_h5 = bool(cfg["paths"].get("transpose_h5", False))

    if phase_mat:
        phase, info = load_phase_mat(phase_mat, phase_var=phase_var, transpose_h5=transpose_h5)
        focal_x_m, _ = load_optional_mat_variable(phase_mat, "focal_x_m", transpose_h5=transpose_h5, squeeze=True)
        focal_y_m, _ = load_optional_mat_variable(phase_mat, "focal_y_m", transpose_h5=transpose_h5, squeeze=True)
        x_um = np.asarray(focal_x_m, dtype=np.float64).reshape(-1) * 1e6 if focal_x_m is not None else None
        y_um = np.asarray(focal_y_m, dtype=np.float64).reshape(-1) * 1e6 if focal_y_m is not None else None
        info["is_smoke_test"] = False
        return phase, info, x_um, y_um

    shape = int(cfg["runtime"]["smoke_shape"])
    mode = cfg["runtime"].get("smoke_phase", "random")
    rng = np.random.default_rng(int(cfg["runtime"].get("random_seed", 1)))
    if mode == "zeros":
        phase = np.zeros((shape, shape), dtype=np.float32)
    elif mode == "random":
        phase = rng.uniform(0, 2 * np.pi, size=(shape, shape)).astype(np.float32)
    else:
        raise ValueError("runtime.smoke_phase must be 'random' or 'zeros'.")
    info = {
        "is_smoke_test": True,
        "message": "No phase MAT supplied; generated a smoke-test phase. This is not an official DOE result.",
        "shape": tuple(int(v) for v in phase.shape),
        "phase_mode": mode,
    }
    return phase, info, None, None


def make_output_dir(cfg: dict, args: argparse.Namespace) -> Path:
    """Create the run output directory."""
    if args.outdir:
        return ensure_unique_dir(args.outdir)
    root = THIS_DIR / cfg["paths"]["out_root"]
    name = f"{timestamp()}_rtad_mraf_gs"
    return ensure_unique_dir(root / name)


def save_target_outputs(outdir: Path, target) -> None:
    """Save target arrays to NPZ and MAT files."""
    np.savez_compressed(
        outdir / "target.npz",
        I_target=target.I_target,
        A_target=target.A_target,
        mask_flat=target.mask_flat,
        mask_edge=target.mask_edge,
        mask_support=target.mask_support,
        mask_bg=target.mask_bg,
        mask_free=target.mask_free,
        x_um=target.x_um,
        y_um=target.y_um,
        params=json.dumps(target.params),
        measured=json.dumps(target.measured),
    )
    savemat(
        outdir / "target.mat",
        {
            "I_target": target.I_target,
            "A_target": target.A_target,
            "mask_flat": target.mask_flat,
            "mask_edge": target.mask_edge,
            "mask_support": target.mask_support,
            "mask_bg": target.mask_bg,
            "mask_free": target.mask_free,
            "x_um": target.x_um,
            "y_um": target.y_um,
            "params_json": json.dumps(target.params, indent=2),
            "measured_json": json.dumps(target.measured, indent=2),
        },
        do_compression=True,
    )


def write_report(outdir: Path, cfg: dict, phase_info: dict, target, result) -> None:
    """Write a plain-text summary report."""
    lines = [
        "RTAD MRAF/GS Python refinement report",
        "",
        f"Backend: {result.backend_description}",
        f"Phase source: {phase_info}",
        "",
        "Target measured size50:",
        f"  W50_x_um = {target.measured['W50_x_um']:.9g}",
        f"  H50_y_um = {target.measured['H50_y_um']:.9g}",
        "",
        "MRAF projection formula:",
        "  signal/support: E' = W * exp(i * angle(E))",
        "  free guard band: E' = mraf_factor * E",
        "  background: keep, attenuate, or zero according to bg_mode",
        "  mraf_factor=0 fully attenuates the free region; mraf_factor=1 keeps it unchanged.",
        "",
        "Final metrics:",
        format_metrics(result.final_metrics),
        "",
        "Config used:",
        json.dumps(cfg, indent=2),
    ]
    (outdir / "report.txt").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    """Run one refinement case."""
    args = parse_args()
    cfg = apply_overrides(CONFIG, args)
    outdir = make_output_dir(cfg, args)

    print(f"Output directory: {outdir}")
    phase0, phase_info, mat_x_um, mat_y_um = load_or_make_phase(cfg)
    if phase_info.get("is_smoke_test"):
        print(phase_info["message"])
    else:
        print(f"Loaded phase: shape={phase0.shape}, reader={phase_info.get('reader')}, var={phase_info.get('variable')}")

    Ny, Nx = phase0.shape
    if Ny != Nx:
        raise ValueError(f"Only square phase grids are supported in this first version; got {phase0.shape}.")

    cfg["grid"]["N"] = int(Nx)
    physical = cfg["physical"]
    focal_dx_um = float(cfg["grid"]["focal_dx_um"])
    focal_dy_um = float(cfg["grid"]["focal_dy_um"])

    if mat_x_um is not None and mat_y_um is not None and mat_x_um.size == Nx and mat_y_um.size == Ny:
        x_um = mat_x_um
        y_um = mat_y_um
        cfg["grid"]["axis_source"] = "phase_mat:focal_x_m/focal_y_m"
        cfg["grid"]["focal_dx_um"] = float(np.median(np.diff(x_um)))
        cfg["grid"]["focal_dy_um"] = float(np.median(np.diff(y_um)))
    else:
        x_um = make_axis_um(Nx, focal_dx_um)
        y_um = make_axis_um(Ny, focal_dy_um)
        cfg["grid"]["axis_source"] = "config:focal_dx_um/focal_dy_um"

    target = make_rtad_rect_target(
        shape=phase0.shape,
        x_um=x_um,
        y_um=y_um,
        **cfg["target"],
    )
    print(
        "Target measured size50: "
        f"W={target.measured['W50_x_um']:.6f} um, H={target.measured['H50_y_um']:.6f} um"
    )

    dx_doe_m = physical["wavelength_m"] * physical["focal_length_m"] / (Nx * (abs(float(np.median(np.diff(x_um)))) * 1e-6))
    backend = get_backend(use_cupy=bool(cfg["runtime"]["use_cupy"]), device_id=int(cfg["runtime"]["device_id"]), verbose=True)
    input_amp = make_input_gaussian(
        shape=phase0.shape,
        dx_doe_m=dx_doe_m,
        gaussian_1e2_diameter_m=physical["input_gaussian_1e2_diameter_m"],
        clear_aperture_m=physical["clear_aperture_m"],
        xp=backend.xp,
        dtype=backend.float_dtype,
    )
    cfg["grid"]["dx_doe_m"] = float(dx_doe_m)

    save_json(outdir / "config_used.json", cfg)
    save_target_outputs(outdir, target)

    result = run_refinement(
        phase0=phase0,
        input_amp=input_amp,
        target_amp=target.A_target,
        masks=target.masks(),
        x_um=target.x_um,
        y_um=target.y_um,
        backend=backend,
        **cfg["refinement"],
    )

    np.save(outdir / "phase_refined.npy", result.phase_refined)
    np.save(outdir / "reconstruction_refined.npy", result.reconstruction_intensity)
    save_metrics_csv(outdir / "metrics.csv", result.metrics_history)
    save_phase_mat(
        outdir / "phase_refined.mat",
        result.phase_refined,
        params=cfg,
        metrics=result.final_metrics,
        reconstruction_intensity=result.reconstruction_intensity,
    )

    masks = target.masks()
    dpi = int(cfg["runtime"]["figure_dpi"])
    plot_target(target, outdir, dpi=dpi)
    plot_masks(target, outdir, dpi=dpi)
    plot_phase(phase0, outdir / "phase0.png", "Initial phase", dpi=dpi)
    plot_phase(result.phase_refined, outdir / "phase_refined.png", "Refined phase", dpi=dpi)
    plot_intensity(
        result.initial_intensity,
        masks,
        target.x_um,
        target.y_um,
        target.params,
        outdir / "initial_reconstruction_intensity.png",
        "Initial reconstruction intensity",
        dpi=dpi,
    )
    plot_intensity(
        result.reconstruction_intensity,
        masks,
        target.x_um,
        target.y_um,
        target.params,
        outdir / "refined_reconstruction_intensity.png",
        "Refined reconstruction intensity",
        dpi=dpi,
    )
    plot_profiles_compare(target, result.initial_intensity, result.reconstruction_intensity, masks, outdir, dpi=dpi)
    plot_convergence(result.metrics_history, outdir, dpi=dpi)
    write_report(outdir, cfg, phase_info, target, result)

    print("Final metrics:")
    print(format_metrics(result.final_metrics))
    print(f"Done. Output directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
