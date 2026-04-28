"""Internal diagnostics for RTAD MRAF/GS refinement."""

from __future__ import annotations

from typing import Any

import numpy as np


def normalize_intensity(intensity: np.ndarray, mask_flat: np.ndarray | None = None, mode: str = "flat_mean") -> np.ndarray:
    """Normalize intensity for profile and display diagnostics."""
    arr = np.asarray(intensity, dtype=np.float64)
    if mode == "max":
        denom = np.nanmax(arr)
    elif mode == "sum":
        denom = np.nansum(arr)
    elif mode == "flat_mean":
        if mask_flat is None or not np.any(mask_flat):
            denom = np.nanmax(arr)
        else:
            denom = np.nanmean(arr[mask_flat])
    else:
        raise ValueError(f"Unknown normalization mode '{mode}'.")
    if not np.isfinite(denom) or abs(denom) < 1e-30:
        return np.zeros_like(arr, dtype=np.float32)
    return (arr / denom).astype(np.float32)


def center_profiles(intensity: np.ndarray, x_um: np.ndarray, y_um: np.ndarray) -> dict[str, np.ndarray]:
    """Return x/y center profiles from an intensity image."""
    x_axis = np.asarray(x_um).reshape(-1)
    y_axis = np.asarray(y_um).reshape(-1)
    ix0 = int(np.argmin(np.abs(x_axis)))
    iy0 = int(np.argmin(np.abs(y_axis)))
    return {
        "x_um": x_axis,
        "y_um": y_axis,
        "x_profile": np.asarray(intensity[iy0, :]),
        "y_profile": np.asarray(intensity[:, ix0]),
    }


def _interp_crossing(x1: float, y1: float, x2: float, y2: float, level: float) -> float:
    """Linearly interpolate one threshold crossing."""
    if abs(y2 - y1) < 1e-30:
        return 0.5 * (x1 + x2)
    return x1 + (level - y1) * (x2 - x1) / (y2 - y1)


def width_at_level(axis: np.ndarray, profile: np.ndarray, level: float = 0.5) -> float:
    """Estimate centered full width where ``profile`` crosses ``level``.

    This assumes a roughly centered, single-lobed profile. It is intentionally
    simple and meant for trend tracking, not final beam metrology.
    """
    x = np.asarray(axis, dtype=np.float64).reshape(-1)
    y = np.asarray(profile, dtype=np.float64).reshape(-1)
    if x.size != y.size or x.size < 3:
        return float("nan")
    center = int(np.argmin(np.abs(x)))
    if not np.isfinite(y[center]) or y[center] < level:
        return float("nan")

    x_left = float("nan")
    for k in range(center, 0, -1):
        if y[k - 1] <= level <= y[k] or y[k] <= level <= y[k - 1]:
            x_left = x[k - 1] if y[k - 1] == level else _interp_crossing(x[k - 1], y[k - 1], x[k], y[k], level)
            break

    x_right = float("nan")
    for k in range(center, x.size - 1):
        if y[k + 1] <= level <= y[k] or y[k] <= level <= y[k + 1]:
            x_right = x[k + 1] if y[k + 1] == level else _interp_crossing(x[k], y[k], x[k + 1], y[k + 1], level)
            break

    if not np.isfinite(x_left) or not np.isfinite(x_right):
        return float("nan")
    return float(x_right - x_left)


def compute_metrics(
    intensity: np.ndarray,
    masks: dict[str, np.ndarray],
    x_um: np.ndarray,
    y_um: np.ndarray,
    iteration: int,
) -> dict[str, Any]:
    """Compute lightweight scalar diagnostics for one reconstruction."""
    I = np.asarray(intensity, dtype=np.float64)
    mask_flat = np.asarray(masks["mask_flat"], dtype=bool)
    mask_edge = np.asarray(masks["mask_edge"], dtype=bool)
    mask_support = np.asarray(masks["mask_support"], dtype=bool)
    mask_free = np.asarray(masks.get("mask_free", np.zeros_like(mask_support)), dtype=bool)
    mask_bg = np.asarray(masks["mask_bg"], dtype=bool) & ~mask_free

    total = float(np.nansum(I))
    flat_vals = I[mask_flat]
    flat_mean = float(np.nanmean(flat_vals)) if flat_vals.size else float("nan")
    flat_std = float(np.nanstd(flat_vals)) if flat_vals.size else float("nan")
    flat_rms = flat_std / flat_mean if np.isfinite(flat_mean) and flat_mean > 0 else float("nan")
    flat_uniformity = 1.0 - flat_rms if np.isfinite(flat_rms) else float("nan")

    In = normalize_intensity(I, mask_flat, mode="flat_mean")
    profiles = center_profiles(In, x_um, y_um)
    edge_max = float(np.nanmax(In[mask_edge])) if np.any(mask_edge) else float("nan")
    flat_max_norm = float(np.nanmax(In[mask_flat])) if np.any(mask_flat) else float("nan")

    return {
        "iteration": int(iteration),
        "flat_uniformity": float(flat_uniformity),
        "flat_rms": float(flat_rms),
        "efficiency_flat": float(np.nansum(I[mask_flat]) / total) if total > 0 else float("nan"),
        "efficiency_support": float(np.nansum(I[mask_support]) / total) if total > 0 else float("nan"),
        "background_fraction": float(np.nansum(I[mask_bg]) / total) if total > 0 else float("nan"),
        "free_fraction": float(np.nansum(I[mask_free]) / total) if total > 0 and np.any(mask_free) else 0.0,
        "max_intensity": float(np.nanmax(I)),
        "mean_intensity_flat": flat_mean,
        "size50_x_um": width_at_level(profiles["x_um"], profiles["x_profile"], level=0.5),
        "size50_y_um": width_at_level(profiles["y_um"], profiles["y_profile"], level=0.5),
        "edge_max_over_flat_mean": edge_max,
        "peak_overshoot": flat_max_norm - 1.0 if np.isfinite(flat_max_norm) else float("nan"),
    }
