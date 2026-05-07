"""Lightweight post-run diagnostics for RTAD MRAF/GS results.

These metrics are evaluation-only. They are not used as GS/MRAF constraints.
The primary normalization is the fixed RTAD flat core mean:

    I_norm = I / mean(I[mask_flat])
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter1d

from .propagation import forward_fft, intensity, make_input_gaussian
from .rtad_target import make_axis_um, make_rtad_rect_target
from .utils import save_json

LEVEL_E2 = float(np.exp(-2.0))
LEVELS = {"90": 0.9, "50": 0.5, "13p5": LEVEL_E2}


@dataclass
class CrossingResult:
    """Width and crossing details for one centered profile level."""

    level: float
    left: float
    right: float
    width: float
    center: float
    valid: bool
    warning: str | None = None


@dataclass
class CaseData:
    """Arrays and configuration loaded from a refinement case directory."""

    case_dir: Path
    config: dict[str, Any]
    intensity_raw: np.ndarray
    x_um: np.ndarray
    y_um: np.ndarray
    mask_flat: np.ndarray
    mask_edge: np.ndarray
    mask_support: np.ndarray
    mask_signal: np.ndarray
    mask_free: np.ndarray
    mask_template_support: np.ndarray
    params: dict[str, Any]
    intensity_source: str


def _loads_jsonish(value: Any) -> dict[str, Any]:
    """Decode JSON saved in NPZ/MAT fields."""
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"U", "S"}:
            value = "".join(str(v) for v in value.ravel())
            return json.loads(value)
        value = value.squeeze()
        if value.shape == ():
            value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        return json.loads(value)
    return dict(value)


def _load_config(case_dir: Path) -> dict[str, Any]:
    """Load ``config_used.json`` if present."""
    path = case_dir / "config_used.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _load_target_from_npz(path: Path) -> dict[str, Any]:
    """Load target axes and masks from ``target.npz``."""
    data = np.load(path, allow_pickle=True)
    params = _loads_jsonish(data["params"]) if "params" in data else {}
    mask_support = np.asarray(data["mask_support"], dtype=bool)
    mask_signal = np.asarray(data["mask_signal"], dtype=bool) if "mask_signal" in data else mask_support
    mask_free = np.asarray(data["mask_free"], dtype=bool) if "mask_free" in data else np.zeros_like(mask_signal)
    mask_template_support = (
        np.asarray(data["mask_template_support"], dtype=bool)
        if "mask_template_support" in data
        else mask_support
    )
    return {
        "x_um": np.asarray(data["x_um"], dtype=np.float64).reshape(-1),
        "y_um": np.asarray(data["y_um"], dtype=np.float64).reshape(-1),
        "mask_flat": np.asarray(data["mask_flat"], dtype=bool),
        "mask_edge": np.asarray(data["mask_edge"], dtype=bool),
        "mask_support": mask_signal,
        "mask_signal": mask_signal,
        "mask_free": mask_free,
        "mask_template_support": mask_template_support,
        "params": params,
    }


def _load_target_from_mat(path: Path) -> dict[str, Any]:
    """Load target axes and masks from ``target.mat``."""
    data = loadmat(path)
    params = _loads_jsonish(data["params_json"]) if "params_json" in data else {}
    mask_support = np.asarray(data["mask_support"], dtype=bool)
    mask_signal = np.asarray(data["mask_signal"], dtype=bool) if "mask_signal" in data else mask_support
    mask_free = np.asarray(data["mask_free"], dtype=bool) if "mask_free" in data else np.zeros_like(mask_signal)
    mask_template_support = (
        np.asarray(data["mask_template_support"], dtype=bool)
        if "mask_template_support" in data
        else mask_support
    )
    return {
        "x_um": np.asarray(data["x_um"], dtype=np.float64).reshape(-1),
        "y_um": np.asarray(data["y_um"], dtype=np.float64).reshape(-1),
        "mask_flat": np.asarray(data["mask_flat"], dtype=bool),
        "mask_edge": np.asarray(data["mask_edge"], dtype=bool),
        "mask_support": mask_signal,
        "mask_signal": mask_signal,
        "mask_free": mask_free,
        "mask_template_support": mask_template_support,
        "params": params,
    }


def _regenerate_target(config: dict[str, Any], shape: tuple[int, int]) -> dict[str, Any]:
    """Regenerate target data from ``config_used.json``."""
    grid = config.get("grid", {})
    target_cfg = config.get("target", {})
    Nx = int(shape[1])
    Ny = int(shape[0])
    x_um = make_axis_um(Nx, float(grid.get("focal_dx_um", 2.5)))
    y_um = make_axis_um(Ny, float(grid.get("focal_dy_um", 2.5)))
    target = make_rtad_rect_target(shape=shape, x_um=x_um, y_um=y_um, **target_cfg)
    return {
        "x_um": target.x_um,
        "y_um": target.y_um,
        "mask_flat": target.mask_flat,
        "mask_edge": target.mask_edge,
        "mask_support": target.mask_support,
        "mask_signal": target.mask_signal,
        "mask_free": target.mask_free,
        "mask_template_support": target.mask_template_support,
        "params": target.params,
    }


def _load_phase_from_case(case_dir: Path) -> np.ndarray:
    """Load refined phase from NPY or MAT."""
    npy = case_dir / "phase_refined.npy"
    if npy.exists():
        return np.asarray(np.load(npy), dtype=np.float32)
    mat = case_dir / "phase_refined.mat"
    if not mat.exists():
        raise FileNotFoundError("No reconstruction_refined.npy, phase_refined.npy, or phase_refined.mat found.")
    data = loadmat(mat)
    for key in ("phase_refined", "phase_refined_rad"):
        if key in data:
            return np.asarray(data[key], dtype=np.float32)
    raise KeyError(f"No phase_refined variable found in {mat}.")


def _recompute_intensity(case_dir: Path, config: dict[str, Any]) -> tuple[np.ndarray, str]:
    """Recompute reconstruction intensity from refined phase using NumPy FFT."""
    phase = _load_phase_from_case(case_dir)
    Ny, Nx = phase.shape
    grid = config.get("grid", {})
    physical = config.get("physical", {})
    focal_dx_um = abs(float(grid.get("focal_dx_um", 2.5)))
    wavelength_m = float(physical.get("wavelength_m", 532e-9))
    focal_length_m = float(physical.get("focal_length_m", 429e-3))
    dx_doe_m = float(grid.get("dx_doe_m", wavelength_m * focal_length_m / (Nx * focal_dx_um * 1e-6)))
    amp = make_input_gaussian(
        shape=phase.shape,
        dx_doe_m=dx_doe_m,
        gaussian_1e2_diameter_m=float(physical.get("input_gaussian_1e2_diameter_m", 5e-3)),
        clear_aperture_m=float(physical.get("clear_aperture_m", 15e-3)),
        xp=np,
        dtype=np.float32,
    )
    field = amp * np.exp(1j * phase).astype(np.complex64)
    return intensity(forward_fft(field, np), np).astype(np.float32), "recomputed from phase_refined"


def load_case_data(case_dir: str | Path) -> CaseData:
    """Load intensity, target masks, axes, and config from a refinement case."""
    case = Path(case_dir)
    config = _load_config(case)

    reconstruction = case / "reconstruction_refined.npy"
    if reconstruction.exists():
        intensity_raw = np.asarray(np.load(reconstruction), dtype=np.float64)
        intensity_source = str(reconstruction)
    else:
        intensity_raw, intensity_source = _recompute_intensity(case, config)
        intensity_raw = np.asarray(intensity_raw, dtype=np.float64)

    if (case / "target.npz").exists():
        target = _load_target_from_npz(case / "target.npz")
    elif (case / "target.mat").exists():
        target = _load_target_from_mat(case / "target.mat")
    else:
        target = _regenerate_target(config, intensity_raw.shape)

    required = ("mask_flat", "mask_edge", "mask_support", "mask_signal", "mask_free", "mask_template_support")
    for key in required:
        if target[key].shape != intensity_raw.shape:
            raise ValueError(f"{key} shape {target[key].shape} does not match intensity shape {intensity_raw.shape}.")

    return CaseData(
        case_dir=case,
        config=config,
        intensity_raw=intensity_raw,
        x_um=np.asarray(target["x_um"], dtype=np.float64).reshape(-1),
        y_um=np.asarray(target["y_um"], dtype=np.float64).reshape(-1),
        mask_flat=np.asarray(target["mask_flat"], dtype=bool),
        mask_edge=np.asarray(target["mask_edge"], dtype=bool),
        mask_support=np.asarray(target["mask_support"], dtype=bool),
        mask_signal=np.asarray(target["mask_signal"], dtype=bool),
        mask_free=np.asarray(target["mask_free"], dtype=bool),
        mask_template_support=np.asarray(target["mask_template_support"], dtype=bool),
        params=dict(target.get("params", {})),
        intensity_source=intensity_source,
    )


def normalize_by_flat_core(intensity_raw: np.ndarray, mask_flat: np.ndarray) -> tuple[np.ndarray, float]:
    """Normalize intensity by fixed RTAD flat core mean."""
    if mask_flat.size == 0 or not np.any(mask_flat):
        raise ValueError("mask_flat is empty; cannot normalize diagnostics.")
    flat_values = np.asarray(intensity_raw, dtype=np.float64)[mask_flat]
    I_ref = float(np.mean(flat_values))
    if not np.isfinite(I_ref) or I_ref <= 0:
        raise ValueError(f"Flat-core mean intensity must be positive; got {I_ref}.")
    return (np.asarray(intensity_raw, dtype=np.float64) / I_ref).astype(np.float32), I_ref


def _interp_crossing(x1: float, y1: float, x2: float, y2: float, level: float) -> float:
    """Linearly interpolate one threshold crossing."""
    if abs(y2 - y1) < 1e-30:
        return 0.5 * (x1 + x2)
    return x1 + (level - y1) * (x2 - x1) / (y2 - y1)


def profile_crossings(axis_um: np.ndarray, profile: np.ndarray, level: float, center_um: float) -> CrossingResult:
    """Find left/right crossings of a centered profile at a level."""
    axis = np.asarray(axis_um, dtype=np.float64).reshape(-1)
    values = np.asarray(profile, dtype=np.float64).reshape(-1)
    center_index = int(np.argmin(np.abs(axis - center_um)))
    if center_index <= 0 or center_index >= axis.size - 1:
        return CrossingResult(level, np.nan, np.nan, np.nan, np.nan, False, "center index is at profile boundary")
    if not np.isfinite(values[center_index]) or values[center_index] < level:
        return CrossingResult(level, np.nan, np.nan, np.nan, np.nan, False, f"profile center is below level {level:.6g}")

    left = np.nan
    for idx in range(center_index, 0, -1):
        y_inner = values[idx]
        y_outer = values[idx - 1]
        if y_inner >= level and y_outer <= level:
            left = _interp_crossing(axis[idx - 1], y_outer, axis[idx], y_inner, level)
            break

    right = np.nan
    for idx in range(center_index, axis.size - 1):
        y_inner = values[idx]
        y_outer = values[idx + 1]
        if y_inner >= level and y_outer <= level:
            right = _interp_crossing(axis[idx], y_inner, axis[idx + 1], y_outer, level)
            break

    if not np.isfinite(left) or not np.isfinite(right) or right <= left:
        return CrossingResult(level, left, right, np.nan, np.nan, False, f"could not find valid level {level:.6g} crossings")

    width = float(right - left)
    center = float(0.5 * (left + right))
    asym = abs((center_um - left) - (right - center_um))
    warning = None
    if asym > max(5.0, 0.05 * width):
        warning = f"level {level:.6g} crossings are asymmetric by {asym:.3f} um"
    return CrossingResult(level, float(left), float(right), width, center, True, warning)


def _center_profiles(I_norm: np.ndarray, x_um: np.ndarray, y_um: np.ndarray, center_x_um: float, center_y_um: float) -> dict[str, Any]:
    """Extract x/y profiles through the expected target center."""
    ix0 = int(np.argmin(np.abs(x_um - center_x_um)))
    iy0 = int(np.argmin(np.abs(y_um - center_y_um)))
    return {
        "x_axis": x_um,
        "y_axis": y_um,
        "x_profile": np.asarray(I_norm[iy0, :], dtype=np.float64),
        "y_profile": np.asarray(I_norm[:, ix0], dtype=np.float64),
        "ix0": ix0,
        "iy0": iy0,
    }


def _outward_profile(axis_um: np.ndarray, profile: np.ndarray, center_um: float, side: str) -> tuple[np.ndarray, np.ndarray]:
    """Return profile from center outward for one side."""
    axis = np.asarray(axis_um, dtype=np.float64)
    values = np.asarray(profile, dtype=np.float64)
    c = int(np.argmin(np.abs(axis - center_um)))
    if side == "left":
        a = axis[: c + 1][::-1]
        v = values[: c + 1][::-1]
        r = center_um - a
    elif side == "right":
        a = axis[c:]
        v = values[c:]
        r = a - center_um
    else:
        raise ValueError("side must be 'left' or 'right'.")
    valid = r >= 0
    return r[valid], v[valid]


def _outward_crossing_radius(r_um: np.ndarray, values: np.ndarray, level: float) -> float:
    """Find first outward radius where profile falls through ``level``."""
    if r_um.size < 2 or values[0] < level:
        return np.nan
    for idx in range(0, r_um.size - 1):
        if values[idx] >= level and values[idx + 1] <= level:
            return _interp_crossing(r_um[idx], values[idx], r_um[idx + 1], values[idx + 1], level)
    return np.nan


def derivative_sidelobe_axis(
    axis_um: np.ndarray,
    profile_norm: np.ndarray,
    center_um: float,
    smooth_sigma_px: float = 2.0,
) -> dict[str, Any]:
    """Compute derivative-based sidelobe scores for one axis.

    The outward profile should decrease monotonically after the 90% crossing.
    Positive ``dI/dr`` in that outer region indicates a re-brightening that may
    correspond to a sidelobe, shoulder, or edge spike.
    """
    details: dict[str, Any] = {}
    scores: list[float] = []
    peaks: list[float] = []
    for side in ("left", "right"):
        r, profile = _outward_profile(axis_um, profile_norm, center_um, side)
        if r.size < 5:
            details[side] = {"r_um": r, "profile": profile, "derivative": np.full_like(r, np.nan), "positive": np.full_like(r, np.nan), "r90": np.nan}
            continue
        smooth = gaussian_filter1d(profile, sigma=float(smooth_sigma_px), mode="nearest")
        deriv = np.gradient(smooth, r, edge_order=1)
        r90 = _outward_crossing_radius(r, smooth, 0.9)
        search = r >= r90 if np.isfinite(r90) else np.zeros_like(r, dtype=bool)
        positive = np.where(search, np.maximum(deriv, 0.0), 0.0)
        if np.any(search):
            span = max(float(r[search][-1] - r[search][0]), 1e-12)
            scores.append(float(np.trapezoid(positive[search], r[search]) / span))
            peaks.append(float(np.max(smooth[search])))
        details[side] = {
            "r_um": r,
            "profile": smooth,
            "derivative": deriv,
            "positive": positive,
            "r90": float(r90),
        }
    return {
        "score": float(np.mean(scores)) if scores else np.nan,
        "peak": float(np.max(peaks)) if peaks else np.nan,
        "details": details,
    }


def compute_diagnostics(data: CaseData) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    """Compute the requested lightweight RTAD diagnostics."""
    I_norm, I_ref = normalize_by_flat_core(data.intensity_raw, data.mask_flat)
    center_x = float(data.params.get("center_x_um", data.config.get("target", {}).get("center_x_um", 0.0)))
    center_y = float(data.params.get("center_y_um", data.config.get("target", {}).get("center_y_um", 0.0)))
    W50 = float(data.params.get("W50_um", data.config.get("target", {}).get("W50_um", 330.0)))
    H50 = float(data.params.get("H50_um", data.config.get("target", {}).get("H50_um", 120.0)))

    profiles = _center_profiles(I_norm, data.x_um, data.y_um, center_x, center_y)
    crossings = {
        "x": {name: profile_crossings(data.x_um, profiles["x_profile"], level, center_x) for name, level in LEVELS.items()},
        "y": {name: profile_crossings(data.y_um, profiles["y_profile"], level, center_y) for name, level in LEVELS.items()},
    }
    warnings: list[str] = []
    for axis_name, axis_crossings in crossings.items():
        for level_name, result in axis_crossings.items():
            if not result.valid:
                warnings.append(f"{axis_name} {level_name}: {result.warning}")
            elif result.warning:
                warnings.append(f"{axis_name} {level_name}: {result.warning}")

    size13x = crossings["x"]["13p5"].width
    size13y = crossings["y"]["13p5"].width
    size90x = crossings["x"]["90"].width
    size90y = crossings["y"]["90"].width

    if np.isfinite(size13x) and np.isfinite(size13y):
        X, Y = np.meshgrid(data.x_um, data.y_um)
        e2_roi = (np.abs(X - center_x) <= size13x / 2.0) & (np.abs(Y - center_y) <= size13y / 2.0)
        total_power = float(np.sum(data.intensity_raw))
        efficiency_e2 = 100.0 * float(np.sum(data.intensity_raw[e2_roi])) / total_power if total_power > 0 else np.nan
    else:
        efficiency_e2 = np.nan
        warnings.append("e^-2 efficiency ROI unavailable because 13.5% size crossing failed")

    I_flat = I_norm[data.mask_flat]
    flat_mean = float(np.mean(I_flat))
    flat_std = float(np.std(I_flat))
    rms_nonuniformity = 100.0 * flat_std / flat_mean if flat_mean > 0 else np.nan

    side_x = derivative_sidelobe_axis(data.x_um, profiles["x_profile"], center_x)
    side_y = derivative_sidelobe_axis(data.y_um, profiles["y_profile"], center_y)

    metrics = {
        "I_ref_flat_mean": float(I_ref),
        "size50_x_um": crossings["x"]["50"].width,
        "size50_y_um": crossings["y"]["50"].width,
        "size50_error_x_um": crossings["x"]["50"].width - W50 if np.isfinite(crossings["x"]["50"].width) else np.nan,
        "size50_error_y_um": crossings["y"]["50"].width - H50 if np.isfinite(crossings["y"]["50"].width) else np.nan,
        "size13p5_x_um": size13x,
        "size13p5_y_um": size13y,
        "size90_x_um": size90x,
        "size90_y_um": size90y,
        "transition_13p5_90_x_um": 0.5 * (size13x - size90x) if np.isfinite(size13x) and np.isfinite(size90x) else np.nan,
        "transition_13p5_90_y_um": 0.5 * (size13y - size90y) if np.isfinite(size13y) and np.isfinite(size90y) else np.nan,
        "rms_nonuniformity_percent": float(rms_nonuniformity),
        "uniformity_rms_percent": float(100.0 - rms_nonuniformity) if np.isfinite(rms_nonuniformity) else np.nan,
        "efficiency_e2_percent": float(efficiency_e2),
        "derivative_sidelobe_score_x": side_x["score"],
        "derivative_sidelobe_score_y": side_y["score"],
        "derivative_sidelobe_peak_x": side_x["peak"],
        "derivative_sidelobe_peak_y": side_y["peak"],
        "center_offset_x_um": crossings["x"]["50"].center - center_x if np.isfinite(crossings["x"]["50"].center) else np.nan,
        "center_offset_y_um": crossings["y"]["50"].center - center_y if np.isfinite(crossings["y"]["50"].center) else np.nan,
        "expected_W50_um": W50,
        "expected_H50_um": H50,
        "expected_center_x_um": center_x,
        "expected_center_y_um": center_y,
    }

    details = {
        "I_norm": I_norm,
        "profiles": profiles,
        "crossings": crossings,
        "derivative_x": side_x,
        "derivative_y": side_y,
        "center_x_um": center_x,
        "center_y_um": center_y,
        "W50_um": W50,
        "H50_um": H50,
    }
    return metrics, details, warnings


def _rect_from_crossings(ax: Any, x_cross: CrossingResult, y_cross: CrossingResult, style: str, color: str, label: str) -> None:
    """Plot a rectangle from x/y crossing results."""
    if not (x_cross.valid and y_cross.valid):
        return
    xx = [x_cross.left, x_cross.right, x_cross.right, x_cross.left, x_cross.left]
    yy = [y_cross.left, y_cross.left, y_cross.right, y_cross.right, y_cross.left]
    ax.plot(xx, yy, style, color=color, linewidth=1.2, label=label)


def _plot_expected_flat(ax: Any, params: dict[str, Any]) -> None:
    """Overlay the fixed RTAD flat core rectangle."""
    a0 = float(params.get("a0_um", params.get("W50_um", 330.0) / 2.0 - params.get("delta_x_um", 15.0)))
    b0 = float(params.get("b0_um", params.get("H50_um", 120.0) / 2.0 - params.get("delta_y_um", 8.0)))
    cx = float(params.get("center_x_um", 0.0))
    cy = float(params.get("center_y_um", 0.0))
    xx = [cx - a0, cx + a0, cx + a0, cx - a0, cx - a0]
    yy = [cy - b0, cy - b0, cy + b0, cy + b0, cy - b0]
    ax.plot(xx, yy, "--", color="white", linewidth=1.1, label="target flat core")


def plot_center_profiles(data: CaseData, details: dict[str, Any], outpath: Path) -> None:
    """Plot x/y center profiles with 90%, 50%, and 13.5% crossings."""
    profiles = details["profiles"]
    crossings = details["crossings"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, axis_name, coord, profile, half50 in [
        (axes[0], "x", data.x_um, profiles["x_profile"], details["W50_um"] / 2.0),
        (axes[1], "y", data.y_um, profiles["y_profile"], details["H50_um"] / 2.0),
    ]:
        ax.plot(coord, profile, linewidth=1.1, label=f"{axis_name} center profile")
        for label, level in [("90%", 0.9), ("50%", 0.5), ("13.5%", LEVEL_E2)]:
            ax.axhline(level, linestyle=":", linewidth=0.9, label=label)
        for level_name, result in crossings[axis_name].items():
            if result.valid:
                ax.axvline(result.left, color="0.25", linestyle="--", linewidth=0.8)
                ax.axvline(result.right, color="0.25", linestyle="--", linewidth=0.8)
        center = details[f"center_{axis_name}_um"]
        ax.axvline(center - half50, color="red", linestyle="-", linewidth=0.9, label="target 50% half-width")
        ax.axvline(center + half50, color="red", linestyle="-", linewidth=0.9)
        ax.set_xlabel(f"{axis_name} / um")
        ax.set_ylabel("I / mean(flat core)")
        ax.grid(True, alpha=0.25)
        ax.set_xlim(center - half50 - 120, center + half50 + 120)
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_intensity_diagnostics(data: CaseData, details: dict[str, Any], outpath: Path) -> None:
    """Plot 2D normalized intensity with target and measured rectangles."""
    I_norm = details["I_norm"]
    crossings = details["crossings"]
    fig, ax = plt.subplots(figsize=(7, 4.6), constrained_layout=True)
    vmax = max(1.5, float(np.nanpercentile(I_norm, 99.5)))
    im = ax.imshow(
        I_norm,
        extent=[float(data.x_um[0]), float(data.x_um[-1]), float(data.y_um[0]), float(data.y_um[-1])],
        origin="lower",
        cmap="magma",
        vmin=0,
        vmax=vmax,
    )
    fig.colorbar(im, ax=ax, label="I / mean(flat core)")
    _plot_expected_flat(ax, data.params)
    _rect_from_crossings(ax, crossings["x"]["50"], crossings["y"]["50"], "-", "cyan", "measured 50%")
    _rect_from_crossings(ax, crossings["x"]["13p5"], crossings["y"]["13p5"], "-", "lime", "measured 13.5%")
    cx = details["center_x_um"]
    cy = details["center_y_um"]
    span_x = max(details["W50_um"] / 2.0 + 150.0, 260.0)
    span_y = max(details["H50_um"] / 2.0 + 120.0, 160.0)
    ax.set_xlim(cx - span_x, cx + span_x)
    ax.set_ylim(cy - span_y, cy + span_y)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("Python diagnostics intensity")
    ax.legend(fontsize=8, loc="upper right")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def plot_derivative_sidelobe(details: dict[str, Any], outpath: Path) -> None:
    """Plot outward profiles and positive derivative sidelobe regions."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for col, axis_name, derivative in [
        (0, "x", details["derivative_x"]),
        (1, "y", details["derivative_y"]),
    ]:
        ax_p = axes[0, col]
        ax_d = axes[1, col]
        for side, color in [("left", "tab:blue"), ("right", "tab:orange")]:
            d = derivative["details"][side]
            r = d["r_um"]
            prof = d["profile"]
            deriv = d["derivative"]
            pos = d["positive"]
            ax_p.plot(r, prof, color=color, linewidth=1.0, label=side)
            ax_d.plot(r, deriv, color=color, linewidth=1.0, label=side)
            ax_d.fill_between(r, 0, pos, where=pos > 0, color=color, alpha=0.22)
            if np.isfinite(d["r90"]):
                ax_p.axvline(d["r90"], color=color, linestyle=":", linewidth=0.9)
                ax_d.axvline(d["r90"], color=color, linestyle=":", linewidth=0.9)
        ax_p.axhline(0.9, color="0.35", linestyle=":", linewidth=0.8)
        ax_p.axhline(0.5, color="0.35", linestyle=":", linewidth=0.8)
        ax_p.axhline(LEVEL_E2, color="0.35", linestyle=":", linewidth=0.8)
        ax_d.axhline(0, color="0.2", linestyle="-", linewidth=0.8)
        ax_p.set_title(f"{axis_name} outward profile")
        ax_d.set_title(f"{axis_name} dI/dr")
        ax_p.set_xlabel("r from center / um")
        ax_d.set_xlabel("r from center / um")
        ax_p.set_ylabel("I / mean(flat core)")
        ax_d.set_ylabel("dI/dr / um^-1")
        ax_p.grid(True, alpha=0.25)
        ax_d.grid(True, alpha=0.25)
        ax_p.legend(fontsize=8)
        ax_d.legend(fontsize=8)
    fig.savefig(outpath, dpi=150)
    plt.close(fig)


def _write_metrics_csv(path: Path, metrics: dict[str, Any]) -> None:
    """Write one-row metrics CSV."""
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics.keys()))
        writer.writeheader()
        writer.writerow(metrics)


def _format_value(value: Any) -> str:
    """Format a report value."""
    if isinstance(value, float):
        return "NaN" if not np.isfinite(value) else f"{value:.9g}"
    return str(value)


def write_diagnostics_report(path: Path, data: CaseData, metrics: dict[str, Any], warnings: list[str]) -> None:
    """Write diagnostics report with definitions and warnings."""
    lines = [
        "RTAD Python diagnostics report",
        "",
        f"Case directory: {data.case_dir}",
        f"Intensity source: {data.intensity_source}",
        "",
        "Normalization:",
        "  I_ref = mean(I[mask_flat]) using the fixed RTAD flat core.",
        "  I_norm = I / I_ref.",
        "  Full-image maximum is not used as the primary normalization.",
        "",
        "Metric definitions:",
        "  size50/size13p5/size90 use x/y center-profile linear interpolation at I_norm levels 0.5, exp(-2), and 0.9.",
        "  RMS uniformity uses only fixed mask_flat: 100 * (1 - std(I_flat) / mean(I_flat)).",
        "  e^-2 efficiency uses raw intensity inside the measured 13.5% rectangle divided by total raw intensity.",
        "  Derivative sidelobe score integrates positive dI/dr after the outward 90% crossing.",
        "  Physical meaning: from center outward, intensity should decrease monotonically; positive derivative means the profile re-brightens and may indicate a side lobe, shoulder, or edge spike.",
        "",
        "Target mask areas:",
        f"  mask_signal pixels: {int(np.count_nonzero(data.mask_signal))}",
        f"  mask_free pixels: {int(np.count_nonzero(data.mask_free))}",
        f"  mask_template_support pixels: {int(np.count_nonzero(data.mask_template_support))}",
        "  mask_support is treated as the constrained signal support in the truncated RTAD Python workflow.",
        "",
        "Metrics:",
    ]
    for key, value in metrics.items():
        lines.append(f"  {key}: {_format_value(value)}")
    lines.append("")
    if warnings:
        lines.append("WARNING:")
        for warning in warnings:
            lines.append(f"  - {warning}")
    else:
        lines.append("Warnings: none")
    path.write_text("\n".join(lines), encoding="utf-8")


def run_case_diagnostics(case_dir: str | Path, output_dir: str | Path | None = None) -> tuple[dict[str, Any], Path]:
    """Run diagnostics for a refinement case and save outputs."""
    data = load_case_data(case_dir)
    metrics, details, warnings = compute_diagnostics(data)
    outdir = Path(output_dir) if output_dir is not None else data.case_dir / "diagnostics_python"
    outdir.mkdir(parents=True, exist_ok=True)

    save_json(outdir / "diagnostics_metrics.json", {"metrics": metrics, "warnings": warnings})
    _write_metrics_csv(outdir / "diagnostics_metrics.csv", metrics)
    write_diagnostics_report(outdir / "diagnostics_report.txt", data, metrics, warnings)
    plot_center_profiles(data, details, outdir / "center_profiles_diagnostics.png")
    plot_intensity_diagnostics(data, details, outdir / "intensity_diagnostics.png")
    plot_derivative_sidelobe(details, outdir / "derivative_sidelobe_diagnostics.png")
    return metrics, outdir
