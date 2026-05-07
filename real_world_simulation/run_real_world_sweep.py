"""Run standalone real-world sweeps for a fixed refined DOE phase."""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path
from typing import Any

import csv
import math
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from config_default import ALL_SWEEPS, CONFIG, SUMMARY_FIELDS, SWEEP_PROFILES
from src.field_models import make_real_input_field
from src.metrics import TargetData, compute_metrics, load_or_generate_target
from src.plotting import (
    plot_intensity_image,
    plot_metric_trends,
    plot_profiles_overlay,
    write_pdf_report,
)
from src.propagation import intensity, propagate_after_doe
from src.utils import center_crop, ensure_unique_dir, load_json, save_json, timestamp, value_label, write_csv


def parse_args() -> argparse.Namespace:
    """Parse CLI options."""
    choices = ["all", *ALL_SWEEPS, "nominal"]
    parser = argparse.ArgumentParser(description="Run fixed-DOE real-world tolerance sweeps.")
    parser.add_argument("--sweep", default="nominal", choices=choices, help="Sweep to run.")
    parser.add_argument(
        "--profile",
        default="mild",
        choices=sorted(SWEEP_PROFILES),
        help="Sweep range profile. 'mild' is the default; 'stress' keeps the broad original ranges.",
    )
    parser.add_argument("--smoke-size", type=int, default=None, help="Use a centered phase crop for a fast smoke run.")
    parser.add_argument("--output-root", default=None, help="Override output root directory.")
    parser.add_argument("--timestamp", default=None, help="Override timestamp prefix for reproducible output names.")
    parser.add_argument("--no-pdf", action="store_true", help="Skip simulation_report.pdf generation.")
    return parser.parse_args()


def merged_config() -> dict[str, Any]:
    """Return config with baseline config_used.json values merged in when present."""
    cfg = copy.deepcopy(CONFIG)
    baseline_config = Path(cfg["paths"]["baseline_config"])
    if baseline_config.exists():
        old = load_json(baseline_config)
        for section in ("physical", "grid", "target"):
            cfg.setdefault(section, {})
            cfg[section].update(old.get(section, {}))
    return cfg


def load_phase(path: str | Path, smoke_size: int | None = None) -> np.ndarray:
    """Load the refined phase as read-only input, optionally center-cropped."""
    phase_mem = np.load(path, mmap_mode="r")
    if phase_mem.ndim != 2:
        raise ValueError(f"Expected a 2D phase array, got shape {phase_mem.shape}.")
    if smoke_size is not None:
        return center_crop(phase_mem, int(smoke_size)).astype(np.float32, copy=False)
    return np.asarray(phase_mem, dtype=np.float32)


def compute_dx_doe_m(cfg: dict[str, Any], shape: tuple[int, int], smoke_size: int | None) -> float:
    """Return DOE-plane sampling from config or Fourier relation."""
    if smoke_size is None and cfg.get("grid", {}).get("dx_doe_m"):
        return float(cfg["grid"]["dx_doe_m"])
    nx = int(shape[1])
    wavelength_m = float(cfg["physical"]["wavelength_m"])
    focal_length_m = float(cfg["physical"]["focal_length_m"])
    focal_dx_m = float(cfg["grid"]["focal_dx_um"]) * 1e-6
    return wavelength_m * focal_length_m / (nx * focal_dx_m)


def load_target(cfg: dict[str, Any], shape: tuple[int, int]) -> TargetData:
    """Load matching target masks or generate them for smoke shape."""
    return load_or_generate_target(
        target_npz=cfg["paths"]["target_npz"],
        shape=shape,
        focal_dx_um=float(cfg["grid"]["focal_dx_um"]),
        focal_dy_um=float(cfg["grid"]["focal_dy_um"]),
        target_config=cfg["target"],
    )


def build_case_params(cfg: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    """Merge nominal case parameters with sweep-specific updates."""
    params = dict(cfg["nominal"])
    params.update(updates)
    return params


def is_nominal_params(cfg: dict[str, Any], params: dict[str, Any]) -> bool:
    """Return whether all real-world perturbations are at nominal values."""
    nominal = cfg["nominal"]
    return all(float(params.get(key, 0.0)) == float(value) for key, value in nominal.items())


def simulate_case(
    phase: np.ndarray,
    target: TargetData,
    cfg: dict[str, Any],
    params: dict[str, Any],
    sweep_name: str,
    sweep_parameter: str,
    sweep_value: Any,
    dx_doe_m: float,
) -> dict[str, Any]:
    """Run one fixed-DOE simulation case and compute diagnostics."""
    physical = cfg["physical"]
    focal_dx_m = float(cfg["grid"]["focal_dx_um"]) * 1e-6
    focal_dy_m = float(cfg["grid"]["focal_dy_um"]) * 1e-6
    input_field, input_meta = make_real_input_field(
        shape=phase.shape,
        dx_doe_m=dx_doe_m,
        wavelength_m=float(physical["wavelength_m"]),
        focal_length_m=float(physical["focal_length_m"]),
        params=params,
    )
    doe_field = input_field * np.exp(1j * phase).astype(np.complex64)
    obs_field = propagate_after_doe(
        doe_field,
        focal_dx_m=focal_dx_m,
        focal_dy_m=focal_dy_m,
        wavelength_m=float(physical["wavelength_m"]),
        defocus_m=float(params.get("defocus_mm", 0.0)) * 1e-3,
    )
    I = intensity(obs_field)
    metrics, details, warnings = compute_metrics(
        I,
        target=target,
        aperture_throughput_percent=input_meta["aperture_throughput_percent"],
    )
    case_label = f"{sweep_parameter}={value_label(sweep_value)}"
    row = {
        "sweep_name": sweep_name,
        "sweep_parameter": sweep_parameter,
        "sweep_value": json.dumps(sweep_value),
        "sweep_value_raw": sweep_value,
        "case_label": case_label,
        **metrics,
    }
    return {
        "row": row,
        "details": details,
        "warnings": warnings,
        "intensity": I,
        "params": params,
        "input_meta": input_meta,
        "is_nominal": is_nominal_params(cfg, params),
        "label": case_label,
    }


def select_representatives(results: list[dict[str, Any]]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Select nominal and worst representative cases."""
    nominal = next((result for result in results if result["is_nominal"]), results[0])

    def score(result: dict[str, Any]) -> float:
        value = result["row"].get("efficiency_e2_percent", np.nan)
        return float(value) if np.isfinite(value) else np.inf

    worst = min(results, key=score)
    return nominal, worst


def run_one_sweep(
    sweep_name: str,
    sweeps: dict[str, list[dict[str, Any]]],
    phase: np.ndarray,
    target: TargetData,
    cfg: dict[str, Any],
    dx_doe_m: float,
    run_dir: Path,
    no_pdf: bool = False,
) -> Path:
    """Run and save one sweep."""
    outdir = ensure_unique_dir(run_dir / sweep_name)
    rows: list[dict[str, Any]] = []
    results: list[dict[str, Any]] = []
    all_warnings: list[str] = []

    for case in sweeps[sweep_name]:
        params = build_case_params(cfg, case["updates"])
        result = simulate_case(
            phase=phase,
            target=target,
            cfg=cfg,
            params=params,
            sweep_name=sweep_name,
            sweep_parameter=str(case["parameter"]),
            sweep_value=case["value"],
            dx_doe_m=dx_doe_m,
        )
        rows.append(result["row"])
        results.append(result)
        for warning in result["warnings"]:
            all_warnings.append(f"{result['label']}: {warning}")
        print(
            f"{sweep_name:14s} {result['label']:28s} "
            f"e2={result['row']['efficiency_e2_percent']:.4g}% "
            f"rms={result['row']['rms_nonuniformity_percent']:.4g}% "
            f"throughput={result['row']['aperture_throughput_percent']:.4g}%"
        )

    nominal, worst = select_representatives(results)
    profile_results = [{"label": result["label"], "details": result["details"]} for result in results]

    write_csv(outdir / "summary.csv", rows, SUMMARY_FIELDS)
    save_json(
        outdir / "summary.json",
        {
            "sweep_name": sweep_name,
            "rows": rows,
            "warnings": all_warnings,
            "config": cfg,
            "dx_doe_m": dx_doe_m,
            "phase_shape": list(phase.shape),
            "representatives": {
                "nominal": nominal["label"],
                "worst": worst["label"],
                "worst_selection": "minimum efficiency_e2_percent",
            },
        },
    )
    plot_profiles_overlay(profile_results, target, outdir / "center_profiles_overlay.png", dpi=int(cfg["runtime"]["figure_dpi"]))
    plot_metric_trends(rows, outdir / "metric_trends.png", dpi=int(cfg["runtime"]["figure_dpi"]))
    for result in results:
        label_safe = result["label"].replace("=", "_").replace(".", "p")
        plot_intensity_image(
            result["intensity"],
            target,
            outdir / f"intensity_{label_safe}.png",
            result["label"],
            dpi=int(cfg["runtime"]["figure_dpi"]),
        )
    if not no_pdf:
        write_pdf_report(
            outdir / "simulation_report.pdf",
            sweep_name=sweep_name,
            rows=rows,
            profile_results=profile_results,
            target=target,
            results=results,
            warnings=all_warnings,
        )
    return outdir


def _read_summary_rows(sweep_dir: Path) -> list[dict[str, Any]]:
    with (sweep_dir / "summary.csv").open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(row: dict[str, Any], key: str) -> float:
    try:
        return float(row[key])
    except Exception:
        return math.nan


def write_run_overview(run_dir: Path, sweep_dirs: list[Path]) -> None:
    """Write run-level CSV and overview figures across sweep directories."""
    summary_rows: list[dict[str, Any]] = []
    for sweep_dir in sweep_dirs:
        rows = _read_summary_rows(sweep_dir)
        if not rows:
            continue
        worst_e = min(
            rows,
            key=lambda row: _float(row, "efficiency_e2_percent")
            if math.isfinite(_float(row, "efficiency_e2_percent"))
            else 1e99,
        )
        worst_r = max(
            rows,
            key=lambda row: _float(row, "rms_nonuniformity_percent")
            if math.isfinite(_float(row, "rms_nonuniformity_percent"))
            else -1e99,
        )
        summary_rows.append(
            {
                "sweep": sweep_dir.name,
                "dir": str(sweep_dir),
                "min_e2_parameter": worst_e["sweep_parameter"],
                "min_e2_value": worst_e["sweep_value"],
                "min_e2_percent": _float(worst_e, "efficiency_e2_percent"),
                "max_rms_parameter": worst_r["sweep_parameter"],
                "max_rms_value": worst_r["sweep_value"],
                "max_rms_percent": _float(worst_r, "rms_nonuniformity_percent"),
                "worst_image": str(sweep_dir / "representative_intensity_worst.png"),
                "metric_trends": str(sweep_dir / "metric_trends.png"),
                "profiles_overlay": str(sweep_dir / "center_profiles_overlay.png"),
                "pdf_report": str(sweep_dir / "simulation_report.pdf"),
            }
        )
    if not summary_rows:
        return

    run_summary = run_dir / "run_effect_summary.csv"
    with run_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    labels = [row["sweep"] for row in summary_rows]
    min_e2 = [float(row["min_e2_percent"]) for row in summary_rows]
    max_rms = [float(row["max_rms_percent"]) for row in summary_rows]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(12, 7.5), constrained_layout=True)
    axes[0].bar(x, min_e2, color="tab:blue")
    axes[0].axhline(92.51767331243764, color="0.35", linestyle="--", linewidth=1, label="nominal e^-2")
    axes[0].set_ylabel("minimum e^-2 efficiency / %")
    axes[0].set_ylim(max(0.0, min(min_e2) - 5.0), 100.0)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=25, ha="right")
    axes[0].grid(True, axis="y", alpha=0.25)
    axes[0].legend()
    axes[1].bar(x, max_rms, color="tab:red")
    axes[1].axhline(1.8610980361700058, color="0.35", linestyle="--", linewidth=1, label="nominal RMS")
    axes[1].set_ylabel("maximum RMS nonuniformity / %")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=25, ha="right")
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].legend()
    fig.savefig(run_dir / "run_effect_metrics_overview.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(3, 3, figsize=(15, 10.5), constrained_layout=True)
    for ax, row in zip(axes.ravel(), summary_rows):
        image = plt.imread(row["worst_image"])
        ax.imshow(image)
        ax.set_title(
            f"{row['sweep']}\nmin e2={float(row['min_e2_percent']):.2f}%, max RMS={float(row['max_rms_percent']):.1f}%",
            fontsize=10,
        )
        ax.axis("off")
    for ax in axes.ravel()[len(summary_rows) :]:
        ax.axis("off")
    fig.savefig(run_dir / "run_worst_intensity_montage.png", dpi=150)
    plt.close(fig)


def main() -> int:
    """Run requested sweeps."""
    args = parse_args()
    cfg = merged_config()
    sweeps = SWEEP_PROFILES[args.profile]
    if args.output_root is not None:
        cfg["paths"]["output_root"] = Path(args.output_root)
    output_root = Path(cfg["paths"]["output_root"])
    output_root.mkdir(parents=True, exist_ok=True)

    phase_path = Path(cfg["paths"]["phase_refined"])
    if not phase_path.exists():
        raise FileNotFoundError(f"Missing fixed phase input: {phase_path}")
    phase = load_phase(phase_path, smoke_size=args.smoke_size)
    target = load_target(cfg, phase.shape)
    dx_doe_m = compute_dx_doe_m(cfg, phase.shape, args.smoke_size)
    stamp = args.timestamp or timestamp()
    sweep_label = "all" if args.sweep == "all" else args.sweep
    run_name = f"{stamp}_{args.profile}_{sweep_label}"
    if args.smoke_size is not None:
        run_name = f"{run_name}_smoke{args.smoke_size}"
    run_dir = ensure_unique_dir(output_root / run_name)

    print(f"Fixed phase input: {phase_path}")
    print(f"Phase shape: {phase.shape}, dx_doe_m={dx_doe_m:.10g}")
    print(f"Output root: {output_root}")
    print(f"Run directory: {run_dir}")
    print(f"Sweep profile: {args.profile}")

    sweep_names = ALL_SWEEPS if args.sweep == "all" else [args.sweep]
    outdirs = [
        run_one_sweep(
            sweep_name=name,
            sweeps=sweeps,
            phase=phase,
            target=target,
            cfg=cfg,
            dx_doe_m=dx_doe_m,
            run_dir=run_dir,
            no_pdf=args.no_pdf,
        )
        for name in sweep_names
    ]
    if len(outdirs) > 1:
        write_run_overview(run_dir, outdirs)
    print("Output directories:")
    for outdir in outdirs:
        print(f"  {outdir}")
    if len(outdirs) > 1:
        print(f"Run overview: {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
