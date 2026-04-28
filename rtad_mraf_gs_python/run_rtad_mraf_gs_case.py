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
from src.diagnostics import run_case_diagnostics
from src.io_mat import load_optional_mat_variable, load_phase_mat, save_phase_mat
from src.mraf_gs import run_refinement
from src.plotting import (
    plot_convergence,
    plot_intensity,
    plot_masks,
    plot_phase,
    plot_profiles_compare,
    plot_profiles_compare_mraf_xy_xonly,
    plot_profiles_compare_initial_mraf_wgs,
    plot_target,
    plot_wgs_x_weights,
    plot_wgs_weights,
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
    parser.add_argument("--mraf-iters", type=int, default=None, help="MRAF iterations for method=mraf_then_wgs.")
    parser.add_argument("--wgs-iters", type=int, default=None, help="WGS iterations for method=mraf_then_wgs.")
    parser.add_argument("--mraf-factor", type=float, default=None, help="MRAF free-region attenuation factor.")
    parser.add_argument("--constraint-mode", default=None, choices=["truncated_rtad", "full_rtad"], help="Target constraint mode.")
    parser.add_argument("--release-level", type=float, default=None, help="Intensity release threshold for truncated RTAD.")
    parser.add_argument("--target-mode", default=None, choices=["separable"], help="Full RTAD template mode.")
    parser.add_argument("--delta-x", type=float, default=None, help="RTAD x falling-edge half width in um.")
    parser.add_argument("--delta-y", type=float, default=None, help="RTAD y falling-edge half width in um.")
    parser.add_argument("--outdir", default=None, help="Output directory. If omitted, timestamped artifacts dir is used.")
    parser.add_argument("--use-cupy", dest="use_cupy", action="store_true", default=None, help="Use CuPy if available.")
    parser.add_argument("--no-cupy", dest="use_cupy", action="store_false", help="Force NumPy CPU.")
    parser.add_argument("--transpose-h5", action="store_true", help="Transpose HDF5 MAT arrays after reading.")
    parser.add_argument("--swap-phase-xy", dest="swap_phase_xy", action="store_true", default=None, help="Transpose the imported phase0 matrix to swap x/y axes.")
    parser.add_argument("--no-swap-phase-xy", dest="swap_phase_xy", action="store_false", help="Do not swap x/y axes on the imported phase0 matrix.")
    parser.add_argument("--metrics-interval", type=int, default=None, help="Metrics logging interval.")
    parser.add_argument("--wgs-after-iters", type=int, default=None, help="Iteration at which WGS updates begin.")
    parser.add_argument("--feedback-exponent", type=float, default=None, help="WGS feedback exponent.")
    parser.add_argument("--wgs-strategy", default=None, choices=["flat_local", "xy_then_x"], help="WGS strategy.")
    parser.add_argument("--wgs-xy-iters", type=int, default=None, help="Flat-local XY WGS iterations for strategy=xy_then_x.")
    parser.add_argument("--wgs-xonly-iters", type=int, default=None, help="X-only WGS iterations for strategy=xy_then_x.")
    parser.add_argument("--wgs-feedback-exponent", type=float, default=None, help="Flat-core WGS amplitude feedback exponent.")
    parser.add_argument("--wgs-xy-feedback-exponent", type=float, default=None, help="XY WGS feedback exponent.")
    parser.add_argument("--wgs-x-feedback-exponent", type=float, default=None, help="X-only WGS feedback exponent.")
    parser.add_argument("--wgs-update-every", type=int, default=None, help="Update WGS weights every N WGS iterations.")
    parser.add_argument("--wgs-xy-update-every", type=int, default=None, help="Update XY WGS weights every N XY iterations.")
    parser.add_argument("--wgs-x-update-every", type=int, default=None, help="Update x-only WGS weights every N x-only iterations.")
    parser.add_argument("--wgs-weight-min", type=float, default=None, help="Minimum WGS flat-core weight.")
    parser.add_argument("--wgs-weight-max", type=float, default=None, help="Maximum WGS flat-core weight.")
    parser.add_argument("--wgs-x-weight-min", type=float, default=None, help="Minimum x-only WGS weight.")
    parser.add_argument("--wgs-x-weight-max", type=float, default=None, help="Maximum x-only WGS weight.")
    parser.add_argument("--wgs-update-mask", default=None, choices=["flat"], help="WGS update mask; first version supports flat only.")
    parser.add_argument("--bg-mode", choices=["keep", "attenuate", "zero"], default=None, help="MRAF background handling.")
    parser.add_argument("--bg-factor", type=float, default=None, help="MRAF background attenuation if bg-mode=attenuate.")
    parser.add_argument("--smoke-shape", type=int, default=None, help="Shape for random/zero smoke phase when no MAT is supplied.")
    parser.add_argument("--skip-diagnostics", action="store_true", help="Do not run lightweight Python diagnostics after refinement.")
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
    if args.swap_phase_xy is not None:
        cfg["paths"]["swap_phase_xy"] = bool(args.swap_phase_xy)
    if args.iters is not None:
        cfg["refinement"]["num_iters"] = args.iters
    if args.mraf_iters is not None:
        cfg["refinement"]["mraf_iters"] = args.mraf_iters
    if args.wgs_iters is not None:
        cfg["refinement"]["wgs_iters"] = args.wgs_iters
    if args.method is not None:
        cfg["refinement"]["method"] = args.method
    if args.mraf_factor is not None:
        cfg["refinement"]["mraf_factor"] = args.mraf_factor
    if args.constraint_mode is not None:
        cfg["target"]["constraint_mode"] = args.constraint_mode
    if args.release_level is not None:
        cfg["target"]["release_level"] = args.release_level
    if args.target_mode is not None:
        cfg["target"]["target_mode"] = args.target_mode
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
    if args.wgs_strategy is not None:
        cfg["refinement"]["wgs_strategy"] = args.wgs_strategy
    if args.wgs_xy_iters is not None:
        cfg["refinement"]["wgs_xy_iters"] = args.wgs_xy_iters
    if args.wgs_xonly_iters is not None:
        cfg["refinement"]["wgs_xonly_iters"] = args.wgs_xonly_iters
    if args.wgs_feedback_exponent is not None:
        cfg["refinement"]["wgs_feedback_exponent"] = args.wgs_feedback_exponent
    if args.wgs_xy_feedback_exponent is not None:
        cfg["refinement"]["wgs_xy_feedback_exponent"] = args.wgs_xy_feedback_exponent
    if args.wgs_x_feedback_exponent is not None:
        cfg["refinement"]["wgs_x_feedback_exponent"] = args.wgs_x_feedback_exponent
    if args.wgs_update_every is not None:
        cfg["refinement"]["wgs_update_every"] = args.wgs_update_every
    if args.wgs_xy_update_every is not None:
        cfg["refinement"]["wgs_xy_update_every"] = args.wgs_xy_update_every
    if args.wgs_x_update_every is not None:
        cfg["refinement"]["wgs_x_update_every"] = args.wgs_x_update_every
    if args.wgs_weight_min is not None:
        cfg["refinement"]["wgs_weight_min"] = args.wgs_weight_min
    if args.wgs_weight_max is not None:
        cfg["refinement"]["wgs_weight_max"] = args.wgs_weight_max
    if args.wgs_x_weight_min is not None:
        cfg["refinement"]["wgs_x_weight_min"] = args.wgs_x_weight_min
    if args.wgs_x_weight_max is not None:
        cfg["refinement"]["wgs_x_weight_max"] = args.wgs_x_weight_max
    if args.wgs_update_mask is not None:
        cfg["refinement"]["wgs_update_mask"] = args.wgs_update_mask
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
    swap_phase_xy = bool(cfg["paths"].get("swap_phase_xy", False))

    if phase_mat:
        phase, info = load_phase_mat(
            phase_mat,
            phase_var=phase_var,
            transpose_h5=transpose_h5,
            swap_xy=swap_phase_xy,
        )
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
    suffix = ""
    method = str(cfg.get("refinement", {}).get("method", "mraf"))
    if method == "mraf_then_wgs":
        suffix += "_mrafwgs"
        if str(cfg.get("refinement", {}).get("wgs_strategy", "")) == "xy_then_x":
            suffix += "_xyx"
    if cfg.get("target", {}).get("constraint_mode") == "truncated_rtad":
        rel = float(cfg.get("target", {}).get("release_level", 0.0))
        suffix += f"_truncI{int(round(rel * 1000)):04d}"
    name = f"{timestamp()}_rtad_mraf_gs{suffix}"
    return ensure_unique_dir(root / name)


def save_target_outputs(outdir: Path, target) -> None:
    """Save target arrays to NPZ and MAT files."""
    np.savez_compressed(
        outdir / "target.npz",
        I_full=target.I_full,
        A_full=target.A_full,
        I_constraint=target.I_constraint,
        A_constraint=target.A_constraint,
        A_signal=target.A_signal,
        I_target=target.I_target,
        A_target=target.A_target,
        mask_flat=target.mask_flat,
        mask_edge_lock=target.mask_edge_lock,
        mask_signal=target.mask_signal,
        mask_template_support=target.mask_template_support,
        mask_guard_window=target.mask_guard_window,
        mask_bg_far=target.mask_bg_far,
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
            "I_full": target.I_full,
            "A_full": target.A_full,
            "I_constraint": target.I_constraint,
            "A_constraint": target.A_constraint,
            "A_signal": target.A_signal,
            "I_target": target.I_target,
            "A_target": target.A_target,
            "mask_flat": target.mask_flat,
            "mask_edge_lock": target.mask_edge_lock,
            "mask_signal": target.mask_signal,
            "mask_template_support": target.mask_template_support,
            "mask_guard_window": target.mask_guard_window,
            "mask_bg_far": target.mask_bg_far,
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
        f"  full_template_size13p5_x_um = {target.measured['full_template_size13p5_x_um']:.9g}",
        f"  full_template_size13p5_y_um = {target.measured['full_template_size13p5_y_um']:.9g}",
        "",
        "Target constraint semantics:",
        "  old/full RTAD constrained every I_full > 0 pixel.",
        "  truncated RTAD constrains only I_full >= release_level.",
        "  the low-intensity full-template tail goes to the MRAF free/noise region.",
        f"  constraint_mode = {target.params.get('constraint_mode')}",
        f"  release_level = {target.params.get('release_level')}",
        f"  mask_signal pixels = {target.params.get('mask_signal_pixels')}",
        f"  mask_free pixels = {target.params.get('mask_free_pixels')}",
        f"  mask_template_support pixels = {target.params.get('mask_template_support_pixels')}",
        "",
        "MRAF projection formula:",
        "  mask_signal: E' = W * exp(i * angle(E))",
        "  mask_free: E' = mraf_factor * E",
        "  mask_bg_far: keep, attenuate by bg_factor, or zero according to bg_mode",
        "  mraf_factor=0 fully attenuates the free region; mraf_factor=1 keeps it unchanged.",
        f"  bg_mode = {cfg['refinement'].get('bg_mode')}, bg_factor = {cfg['refinement'].get('bg_factor')}",
        "",
        "WGS refinement:",
        f"  method = {cfg['refinement'].get('method')}",
        f"  wgs_strategy = {cfg['refinement'].get('wgs_strategy')}",
        f"  mraf_iters = {cfg['refinement'].get('mraf_iters')}",
        f"  wgs_iters = {cfg['refinement'].get('wgs_iters')}",
        f"  wgs_xy_iters = {cfg['refinement'].get('wgs_xy_iters')}",
        f"  wgs_xonly_iters = {cfg['refinement'].get('wgs_xonly_iters')}",
        f"  wgs_update_mask = {cfg['refinement'].get('wgs_update_mask')} (mask_flat only)",
        f"  wgs_feedback = {cfg['refinement'].get('wgs_feedback')}",
        f"  wgs_feedback_exponent = {cfg['refinement'].get('wgs_feedback_exponent')}",
        f"  wgs_xy_feedback_exponent = {cfg['refinement'].get('wgs_xy_feedback_exponent')}",
        f"  wgs_x_feedback_exponent = {cfg['refinement'].get('wgs_x_feedback_exponent')}",
        f"  wgs_update_every = {cfg['refinement'].get('wgs_update_every')}",
        f"  wgs_xy_update_every = {cfg['refinement'].get('wgs_xy_update_every')}",
        f"  wgs_x_update_every = {cfg['refinement'].get('wgs_x_update_every')}",
        f"  wgs_weight_clip = [{cfg['refinement'].get('wgs_weight_min')}, {cfg['refinement'].get('wgs_weight_max')}]",
        f"  wgs_x_weight_clip = [{cfg['refinement'].get('wgs_x_weight_min')}, {cfg['refinement'].get('wgs_x_weight_max')}]",
        f"  wgs_normalize_weights = {cfg['refinement'].get('wgs_normalize_weights')}",
        f"  wgs_x_normalize = {cfg['refinement'].get('wgs_x_normalize')}",
        "  WGS updates only mask_flat weights. It does not update mask_edge_lock, mask_free, or mask_bg_far.",
        "  MRAF first allocates energy between the RTAD signal/free/background regions; WGS then compensates local flat-core brightness nonuniformity by changing target amplitude weights.",
        "  xy_then_x stage 1 applies 2D flat-local XY feedback. Stage 2 freezes y-local structure and updates a 1D x weight from the flat-core column-mean amplitude.",
        "  Diagnostics metrics remain post-processing only; they are not used as constraints.",
        "",
        "Run warnings/info:",
        "\n".join(f"  {w}" for w in (result.warnings or [])) if result.warnings else "  none",
        "",
        "WGS final weight stats:",
        format_metrics(result.wgs_weight_stats_final or {}),
        "",
        "MRAF-end metrics:",
        format_metrics(result.mraf_metrics or {}),
        "",
        "WGS-XY-end metrics:",
        format_metrics(result.wgs_xy_metrics or {}),
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
        target_amp=target.A_signal,
        masks=target.masks(),
        x_um=target.x_um,
        y_um=target.y_um,
        backend=backend,
        **cfg["refinement"],
    )

    np.save(outdir / "phase_refined.npy", result.phase_refined)
    np.save(outdir / "reconstruction_refined.npy", result.reconstruction_intensity)
    if result.phase_after_mraf is not None:
        np.save(outdir / "phase_after_mraf.npy", result.phase_after_mraf)
    if result.reconstruction_after_mraf is not None:
        np.save(outdir / "reconstruction_after_mraf.npy", result.reconstruction_after_mraf)
    if result.phase_after_wgs_xy is not None:
        np.save(outdir / "phase_after_wgs_xy.npy", result.phase_after_wgs_xy)
    if result.reconstruction_after_wgs_xy is not None:
        np.save(outdir / "reconstruction_after_wgs_xy.npy", result.reconstruction_after_wgs_xy)
    if result.wgs_weights_final is not None:
        np.save(outdir / "wgs_weights_final.npy", result.wgs_weights_final)
    if result.wgs_weights_2d_final is not None:
        np.save(outdir / "wgs_weights_2d_final.npy", result.wgs_weights_2d_final)
    if result.wgs_x_weights_final is not None:
        np.save(outdir / "wgs_x_weights_final.npy", result.wgs_x_weights_final)
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
    if result.reconstruction_after_mraf is not None:
        plot_profiles_compare_initial_mraf_wgs(
            target,
            result.initial_intensity,
            result.reconstruction_after_mraf,
            result.reconstruction_intensity,
            masks,
            outdir,
            dpi=dpi,
        )
    if result.reconstruction_after_mraf is not None and result.reconstruction_after_wgs_xy is not None:
        plot_profiles_compare_mraf_xy_xonly(
            target,
            result.reconstruction_after_mraf,
            result.reconstruction_after_wgs_xy,
            result.reconstruction_intensity,
            masks,
            outdir,
            dpi=dpi,
        )
    if result.wgs_weights_2d_final is not None:
        plot_wgs_weights(result.wgs_weights_2d_final, masks["mask_flat"], target.x_um, target.y_um, outdir, dpi=dpi)
    elif result.wgs_weights_final is not None:
        plot_wgs_weights(result.wgs_weights_final, masks["mask_flat"], target.x_um, target.y_um, outdir, dpi=dpi)
    if result.wgs_x_weights_final is not None:
        plot_wgs_x_weights(target.x_um, result.wgs_x_weights_final, masks["mask_flat"], target.params, outdir, dpi=dpi)
    plot_convergence(result.metrics_history, outdir, dpi=dpi)
    write_report(outdir, cfg, phase_info, target, result)
    if not args.skip_diagnostics:
        diag_metrics, diag_outdir = run_case_diagnostics(outdir)
        print(f"Python diagnostics output directory: {diag_outdir}")
        print("Python diagnostics metrics:")
        print(format_metrics(diag_metrics))

    print("Final metrics:")
    print(format_metrics(result.final_metrics))
    print(f"Done. Output directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
