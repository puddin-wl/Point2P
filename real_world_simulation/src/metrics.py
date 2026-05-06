"""Target loading/generation and simulation metrics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .utils import jsonish_to_dict


LEVEL_E2 = float(np.exp(-2.0))


@dataclass
class CrossingResult:
    """One centered profile crossing result."""

    left: float
    right: float
    width: float
    center: float
    valid: bool
    warning: str = ""


@dataclass
class TargetData:
    """Target axes, masks, and geometry metadata."""

    x_um: np.ndarray
    y_um: np.ndarray
    mask_flat: np.ndarray
    mask_signal: np.ndarray
    mask_edge: np.ndarray
    mask_free: np.ndarray
    mask_bg: np.ndarray
    I_full: np.ndarray
    params: dict[str, Any]


def make_axis_um(length: int, dx_um: float) -> np.ndarray:
    """Create an FFT-shifted coordinate axis with zero at ``N//2``."""
    return (np.arange(length, dtype=np.float64) - length // 2) * float(dx_um)


def raised_cosine_edge(u: np.ndarray, u0: float, u1: float) -> np.ndarray:
    """Evaluate a one-dimensional raised-cosine falling edge."""
    if u1 <= u0:
        raise ValueError(f"u1 must be larger than u0; got {u0} and {u1}.")
    values = np.zeros_like(u, dtype=np.float32)
    values[u <= u0] = 1.0
    idx = (u > u0) & (u < u1)
    t = (u[idx] - u0) / (u1 - u0)
    values[idx] = 0.5 * (1.0 + np.cos(np.pi * t))
    return values


def generate_target(
    shape: tuple[int, int],
    focal_dx_um: float,
    focal_dy_um: float,
    target_config: dict[str, Any],
) -> TargetData:
    """Generate the RTAD rectangular target and masks without external code."""
    ny, nx = int(shape[0]), int(shape[1])
    x_um = make_axis_um(nx, focal_dx_um)
    y_um = make_axis_um(ny, focal_dy_um)
    X, Y = np.meshgrid(x_um, y_um)

    W50 = float(target_config.get("W50_um", 330.0))
    H50 = float(target_config.get("H50_um", 120.0))
    delta_x = float(target_config.get("delta_x_um", 15.0))
    delta_y = float(target_config.get("delta_y_um", 8.0))
    guard_x = float(target_config.get("guard_x_um", 20.0))
    guard_y = float(target_config.get("guard_y_um", 12.0))
    release = float(target_config.get("release_level", LEVEL_E2))
    center_x = float(target_config.get("center_x_um", 0.0))
    center_y = float(target_config.get("center_y_um", 0.0))

    a50 = W50 / 2.0
    b50 = H50 / 2.0
    a0 = a50 - delta_x
    a1 = a50 + delta_x
    b0 = b50 - delta_y
    b1 = b50 + delta_y
    a2 = a1 + guard_x
    b2 = b1 + guard_y
    abs_x = np.abs(X - center_x)
    abs_y = np.abs(Y - center_y)
    Ix = raised_cosine_edge(abs_x.astype(np.float32), a0, a1)
    Iy = raised_cosine_edge(abs_y.astype(np.float32), b0, b1)
    I_full = np.clip(Ix * Iy, 0.0, 1.0).astype(np.float32)

    mask_flat = (abs_x <= a0) & (abs_y <= b0)
    mask_template_support = I_full > 0
    if str(target_config.get("constraint_mode", "truncated_rtad")).lower() == "full_rtad":
        mask_signal = mask_template_support.copy()
    else:
        mask_signal = I_full >= release
    mask_signal = mask_signal | mask_flat
    mask_edge = mask_signal & ~mask_flat
    mask_guard = (abs_x <= a2) & (abs_y <= b2)
    mask_free = mask_guard & ~mask_signal
    mask_bg = ~(mask_signal | mask_free)

    params = dict(target_config)
    params.update(
        {
            "W50_um": W50,
            "H50_um": H50,
            "delta_x_um": delta_x,
            "delta_y_um": delta_y,
            "guard_x_um": guard_x,
            "guard_y_um": guard_y,
            "release_level": release,
            "center_x_um": center_x,
            "center_y_um": center_y,
            "a0_um": a0,
            "a50_um": a50,
            "a1_um": a1,
            "a2_um": a2,
            "b0_um": b0,
            "b50_um": b50,
            "b1_um": b1,
            "b2_um": b2,
        }
    )
    return TargetData(
        x_um=x_um,
        y_um=y_um,
        mask_flat=mask_flat,
        mask_signal=mask_signal,
        mask_edge=mask_edge,
        mask_free=mask_free,
        mask_bg=mask_bg,
        I_full=I_full,
        params=params,
    )


def load_or_generate_target(
    target_npz: str | Path,
    shape: tuple[int, int],
    focal_dx_um: float,
    focal_dy_um: float,
    target_config: dict[str, Any],
) -> TargetData:
    """Load matching baseline target data or regenerate for smoke/crop shape."""
    path = Path(target_npz)
    if path.exists():
        data = np.load(path, allow_pickle=True)
        if "mask_flat" in data and tuple(data["mask_flat"].shape) == tuple(shape):
            params = jsonish_to_dict(data["params"]) if "params" in data else dict(target_config)
            return TargetData(
                x_um=np.asarray(data["x_um"], dtype=np.float64).reshape(-1),
                y_um=np.asarray(data["y_um"], dtype=np.float64).reshape(-1),
                mask_flat=np.asarray(data["mask_flat"], dtype=bool),
                mask_signal=np.asarray(data["mask_signal"], dtype=bool)
                if "mask_signal" in data
                else np.asarray(data["mask_support"], dtype=bool),
                mask_edge=np.asarray(data["mask_edge"], dtype=bool),
                mask_free=np.asarray(data["mask_free"], dtype=bool)
                if "mask_free" in data
                else np.zeros(shape, dtype=bool),
                mask_bg=np.asarray(data["mask_bg"], dtype=bool)
                if "mask_bg" in data
                else np.zeros(shape, dtype=bool),
                I_full=np.asarray(data["I_full"], dtype=np.float32)
                if "I_full" in data
                else np.zeros(shape, dtype=np.float32),
                params=params,
            )
    return generate_target(shape, focal_dx_um, focal_dy_um, target_config)


def normalize_by_flat(intensity: np.ndarray, mask_flat: np.ndarray) -> tuple[np.ndarray, float]:
    """Normalize intensity by fixed flat-core mean."""
    values = np.asarray(intensity, dtype=np.float64)[mask_flat]
    if values.size == 0:
        raise ValueError("mask_flat is empty.")
    mean = float(np.mean(values))
    if not np.isfinite(mean) or mean <= 0:
        return np.zeros_like(intensity, dtype=np.float32), mean
    return (np.asarray(intensity, dtype=np.float64) / mean).astype(np.float32), mean


def _interp_crossing(x1: float, y1: float, x2: float, y2: float, level: float) -> float:
    if abs(y2 - y1) < 1e-30:
        return 0.5 * (x1 + x2)
    return x1 + (level - y1) * (x2 - x1) / (y2 - y1)


def profile_crossings(axis: np.ndarray, profile: np.ndarray, level: float, center: float = 0.0) -> CrossingResult:
    """Find left/right threshold crossings around the expected center."""
    x = np.asarray(axis, dtype=np.float64).reshape(-1)
    y = np.asarray(profile, dtype=np.float64).reshape(-1)
    if x.size != y.size or x.size < 3:
        return CrossingResult(np.nan, np.nan, np.nan, np.nan, False, "bad profile length")
    idx0 = int(np.argmin(np.abs(x - center)))
    if not np.isfinite(y[idx0]) or y[idx0] < level:
        return CrossingResult(np.nan, np.nan, np.nan, np.nan, False, "center below level")

    left = np.nan
    for idx in range(idx0, 0, -1):
        y0 = y[idx]
        y1 = y[idx - 1]
        if (y1 <= level <= y0) or (y0 <= level <= y1):
            left = x[idx - 1] if y1 == level else _interp_crossing(x[idx - 1], y1, x[idx], y0, level)
            break

    right = np.nan
    for idx in range(idx0, x.size - 1):
        y0 = y[idx]
        y1 = y[idx + 1]
        if (y1 <= level <= y0) or (y0 <= level <= y1):
            right = x[idx + 1] if y1 == level else _interp_crossing(x[idx], y0, x[idx + 1], y1, level)
            break

    if not np.isfinite(left) or not np.isfinite(right):
        return CrossingResult(left, right, np.nan, np.nan, False, "crossing not found")
    width = float(right - left)
    crossing_center = float(0.5 * (left + right))
    return CrossingResult(float(left), float(right), width, crossing_center, True, "")


def center_profiles(intensity_norm: np.ndarray, target: TargetData) -> dict[str, np.ndarray]:
    """Return x/y profiles through the expected target center."""
    center_x = float(target.params.get("center_x_um", 0.0))
    center_y = float(target.params.get("center_y_um", 0.0))
    ix0 = int(np.argmin(np.abs(target.x_um - center_x)))
    iy0 = int(np.argmin(np.abs(target.y_um - center_y)))
    return {
        "x_um": target.x_um,
        "y_um": target.y_um,
        "x_profile": np.asarray(intensity_norm[iy0, :], dtype=np.float32),
        "y_profile": np.asarray(intensity_norm[:, ix0], dtype=np.float32),
    }


def compute_metrics(
    intensity: np.ndarray,
    target: TargetData,
    aperture_throughput_percent: float,
) -> tuple[dict[str, Any], dict[str, Any], list[str]]:
    """Compute scalar diagnostics for one simulated case."""
    I = np.asarray(intensity, dtype=np.float64)
    I_norm, flat_mean = normalize_by_flat(I, target.mask_flat)
    profiles = center_profiles(I_norm, target)
    center_x = float(target.params.get("center_x_um", 0.0))
    center_y = float(target.params.get("center_y_um", 0.0))

    crossings_x = {
        "90": profile_crossings(target.x_um, profiles["x_profile"], 0.9, center_x),
        "50": profile_crossings(target.x_um, profiles["x_profile"], 0.5, center_x),
        "13p5": profile_crossings(target.x_um, profiles["x_profile"], LEVEL_E2, center_x),
    }
    crossings_y = {
        "90": profile_crossings(target.y_um, profiles["y_profile"], 0.9, center_y),
        "50": profile_crossings(target.y_um, profiles["y_profile"], 0.5, center_y),
        "13p5": profile_crossings(target.y_um, profiles["y_profile"], LEVEL_E2, center_y),
    }
    warnings: list[str] = []
    for axis_name, crossings in (("x", crossings_x), ("y", crossings_y)):
        for level_name, result in crossings.items():
            if not result.valid:
                warnings.append(f"{axis_name} {level_name}: {result.warning}")

    size13x = crossings_x["13p5"].width
    size13y = crossings_y["13p5"].width
    size90x = crossings_x["90"].width
    size90y = crossings_y["90"].width
    total_power = float(np.sum(I))
    if np.isfinite(size13x) and np.isfinite(size13y) and total_power > 0:
        x_roi = np.abs(target.x_um - center_x) <= size13x / 2.0
        y_roi = np.abs(target.y_um - center_y) <= size13y / 2.0
        efficiency_e2 = 100.0 * float(np.sum(I[np.ix_(y_roi, x_roi)])) / total_power
    else:
        efficiency_e2 = np.nan
        warnings.append("e^-2 efficiency ROI unavailable")

    flat_values = I_norm[target.mask_flat]
    flat_std = float(np.std(flat_values)) if flat_values.size else np.nan
    flat_avg = float(np.mean(flat_values)) if flat_values.size else np.nan
    rms_nonuniformity = 100.0 * flat_std / flat_avg if flat_avg > 0 else np.nan

    metrics = {
        "I_ref_flat_mean": float(flat_mean),
        "size50_x_um": crossings_x["50"].width,
        "size50_y_um": crossings_y["50"].width,
        "size13p5_x_um": size13x,
        "size13p5_y_um": size13y,
        "size90_x_um": size90x,
        "size90_y_um": size90y,
        "transition_13p5_90_x_um": 0.5 * (size13x - size90x)
        if np.isfinite(size13x) and np.isfinite(size90x)
        else np.nan,
        "transition_13p5_90_y_um": 0.5 * (size13y - size90y)
        if np.isfinite(size13y) and np.isfinite(size90y)
        else np.nan,
        "rms_nonuniformity_percent": float(rms_nonuniformity),
        "efficiency_e2_percent": float(efficiency_e2),
        "aperture_throughput_percent": float(aperture_throughput_percent),
        "center_offset_x_um": crossings_x["50"].center - center_x
        if np.isfinite(crossings_x["50"].center)
        else np.nan,
        "center_offset_y_um": crossings_y["50"].center - center_y
        if np.isfinite(crossings_y["50"].center)
        else np.nan,
    }
    details = {
        "I_norm": I_norm,
        "profiles": profiles,
        "crossings_x": crossings_x,
        "crossings_y": crossings_y,
        "warnings": warnings,
    }
    return metrics, details, warnings

