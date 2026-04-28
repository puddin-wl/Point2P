"""Raised-cosine RTAD rectangular target generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .metrics import width_at_level


@dataclass
class RTADTarget:
    """Container for target arrays, masks, axes, and geometry metadata."""

    I_full: np.ndarray
    A_full: np.ndarray
    I_constraint: np.ndarray
    A_constraint: np.ndarray
    A_signal: np.ndarray
    I_target: np.ndarray
    A_target: np.ndarray
    mask_flat: np.ndarray
    mask_edge_lock: np.ndarray
    mask_signal: np.ndarray
    mask_template_support: np.ndarray
    mask_guard_window: np.ndarray
    mask_bg_far: np.ndarray
    mask_edge: np.ndarray
    mask_support: np.ndarray
    mask_bg: np.ndarray
    mask_free: np.ndarray
    x_um: np.ndarray
    y_um: np.ndarray
    params: dict[str, Any]
    profiles: dict[str, np.ndarray]
    measured: dict[str, float]

    def masks(self) -> dict[str, np.ndarray]:
        """Return masks in the form expected by the refinement module."""
        return {
            "mask_flat": self.mask_flat,
            "mask_edge_lock": self.mask_edge_lock,
            "mask_signal": self.mask_signal,
            "mask_template_support": self.mask_template_support,
            "mask_guard_window": self.mask_guard_window,
            "mask_bg_far": self.mask_bg_far,
            "mask_edge": self.mask_edge,  # Compatibility alias for edge_lock.
            "mask_support": self.mask_support,  # Compatibility alias for constrained signal support.
            "mask_bg": self.mask_bg,  # Compatibility alias for far background.
            "mask_free": self.mask_free,
        }


def raised_cosine_edge(u: np.ndarray, u0: float, u1: float) -> np.ndarray:
    """Evaluate a one-dimensional raised-cosine falling edge.

    ``C=1`` for ``u <= u0``, ``C=0`` for ``u >= u1``, and a cosine transition
    in between.
    """
    if u1 <= u0:
        raise ValueError(f"u1 must be larger than u0; got u0={u0}, u1={u1}.")
    C = np.zeros_like(u, dtype=np.float32)
    C[u <= u0] = 1.0
    idx = (u > u0) & (u < u1)
    t = (u[idx] - u0) / (u1 - u0)
    C[idx] = 0.5 * (1.0 + np.cos(np.pi * t))
    return C


def make_axis_um(length: int, dx_um: float) -> np.ndarray:
    """Create an FFT-shifted coordinate axis with zero at index ``N//2``."""
    return (np.arange(length, dtype=np.float64) - length // 2) * float(dx_um)


def make_rtad_rect_target(
    shape: tuple[int, int] | int,
    dx_um: float | None = None,
    dy_um: float | None = None,
    x_um: np.ndarray | None = None,
    y_um: np.ndarray | None = None,
    W50_um: float = 330.0,
    H50_um: float = 120.0,
    delta_x_um: float = 15.0,
    delta_y_um: float = 8.0,
    guard_x_um: float = 20.0,
    guard_y_um: float = 12.0,
    release_level: float = float(np.exp(-2.0)),
    constraint_mode: str = "truncated_rtad",
    center_x_um: float = 0.0,
    center_y_um: float = 0.0,
    target_mode: str = "separable",
    mode: str | None = None,
) -> RTADTarget:
    """Generate a separable raised-cosine RTAD rectangular target.

    The requested ``W50_um`` and ``H50_um`` are intensity 50 percent full sizes.
    For ``constraint_mode='truncated_rtad'``, GS/MRAF only constrains the
    ``I_full >= release_level`` signal region. The lower-intensity tail enters
    the MRAF free/noise region.
    Matrix convention is row=y and col=x.
    """
    if isinstance(shape, int):
        Ny = Nx = int(shape)
    else:
        Ny, Nx = int(shape[0]), int(shape[1])

    if x_um is None or y_um is None:
        if dx_um is None or dy_um is None:
            raise ValueError("Either x_um/y_um or dx_um/dy_um must be provided.")
        x_axis = make_axis_um(Nx, dx_um)
        y_axis = make_axis_um(Ny, dy_um)
    else:
        x_axis = np.asarray(x_um, dtype=np.float64).reshape(-1)
        y_axis = np.asarray(y_um, dtype=np.float64).reshape(-1)
        if x_axis.size != Nx or y_axis.size != Ny:
            raise ValueError(
                f"Axis lengths do not match shape: x={x_axis.size}, y={y_axis.size}, shape={(Ny, Nx)}."
            )

    a50 = float(W50_um) / 2.0
    b50 = float(H50_um) / 2.0
    if delta_x_um <= 0 or delta_y_um <= 0:
        raise ValueError("delta_x_um and delta_y_um must be positive.")
    if delta_x_um >= a50 or delta_y_um >= b50:
        raise ValueError("delta_x_um/delta_y_um would make the flat core non-positive.")

    a0 = a50 - float(delta_x_um)
    a1 = a50 + float(delta_x_um)
    b0 = b50 - float(delta_y_um)
    b1 = b50 + float(delta_y_um)
    a2 = a1 + float(guard_x_um)
    b2 = b1 + float(guard_y_um)

    X, Y = np.meshgrid(x_axis - center_x_um, y_axis - center_y_um)
    absX = np.abs(X)
    absY = np.abs(Y)

    mode_l = (mode if mode is not None else target_mode).lower()
    if mode_l != "separable":
        raise ValueError("Only target_mode='separable' is implemented in this Python version.")
    constraint_mode_l = constraint_mode.lower()
    if constraint_mode_l not in {"truncated_rtad", "full_rtad"}:
        raise ValueError("constraint_mode must be 'truncated_rtad' or 'full_rtad'.")
    if release_level <= 0 or release_level >= 1:
        raise ValueError("release_level must be an intensity value between 0 and 1.")

    Ix = raised_cosine_edge(absX, a0, a1)
    Iy = raised_cosine_edge(absY, b0, b1)
    I_full = np.clip(Ix * Iy, 0.0, 1.0).astype(np.float32)
    A_full = np.sqrt(I_full).astype(np.float32)

    mask_flat = (absX <= a0) & (absY <= b0)
    mask_template_support = I_full > 0
    if constraint_mode_l == "full_rtad":
        mask_signal = mask_template_support.copy()
    else:
        mask_signal = I_full >= float(release_level)
    mask_signal = mask_signal | mask_flat
    mask_edge_lock = mask_signal & ~mask_flat
    mask_guard_window = (absX <= a2) & (absY <= b2)
    mask_free = mask_guard_window & ~mask_signal
    mask_bg_far = ~(mask_signal | mask_free)

    I_constraint = I_full.copy()
    I_constraint[~mask_signal] = np.nan
    A_constraint = A_full.copy()
    A_constraint[~mask_signal] = np.nan
    A_signal = A_full.copy()
    A_signal[~mask_signal] = 0.0

    # Compatibility aliases. In the truncated target, mask_support means the
    # constrained signal support, not the full template support.
    I_target = I_constraint
    A_target = A_signal
    mask_edge = mask_edge_lock
    mask_support = mask_signal
    mask_bg = mask_bg_far

    ix0 = int(np.argmin(np.abs(x_axis - center_x_um)))
    iy0 = int(np.argmin(np.abs(y_axis - center_y_um)))
    profile_x_I = I_full[iy0, :].astype(np.float32)
    profile_y_I = I_full[:, ix0].astype(np.float32)
    profile_x_A = A_full[iy0, :].astype(np.float32)
    profile_y_A = A_full[:, ix0].astype(np.float32)

    x_rel = (x_axis - center_x_um).astype(np.float64)
    y_rel = (y_axis - center_y_um).astype(np.float64)
    measured = {
        "full_template_size50_x_um": width_at_level(x_rel, profile_x_I, level=0.5),
        "full_template_size50_y_um": width_at_level(y_rel, profile_y_I, level=0.5),
        "full_template_size13p5_x_um": width_at_level(x_rel, profile_x_I, level=float(np.exp(-2.0))),
        "full_template_size13p5_y_um": width_at_level(y_rel, profile_y_I, level=float(np.exp(-2.0))),
        "signal_width_x_um": width_at_level(x_rel, profile_x_I, level=float(release_level)),
        "signal_width_y_um": width_at_level(y_rel, profile_y_I, level=float(release_level)),
        "W90_x_um": width_at_level(x_rel, profile_x_I, level=0.9),
        "H90_y_um": width_at_level(y_rel, profile_y_I, level=0.9),
    }
    measured["W50_x_um"] = measured["full_template_size50_x_um"]
    measured["H50_y_um"] = measured["full_template_size50_y_um"]
    params = {
        "target_mode": mode_l,
        "mode": mode_l,
        "constraint_mode": constraint_mode_l,
        "release_level": float(release_level),
        "Ny": Ny,
        "Nx": Nx,
        "dx_um": float(np.median(np.diff(x_axis))) if Nx > 1 else float(dx_um or 0),
        "dy_um": float(np.median(np.diff(y_axis))) if Ny > 1 else float(dy_um or 0),
        "W50_um": float(W50_um),
        "H50_um": float(H50_um),
        "a50_um": a50,
        "b50_um": b50,
        "delta_x_um": float(delta_x_um),
        "delta_y_um": float(delta_y_um),
        "a0_um": a0,
        "a1_um": a1,
        "b0_um": b0,
        "b1_um": b1,
        "guard_x_um": float(guard_x_um),
        "guard_y_um": float(guard_y_um),
        "a2_um": a2,
        "b2_um": b2,
        "center_x_um": float(center_x_um),
        "center_y_um": float(center_y_um),
        "mask_signal_pixels": int(np.count_nonzero(mask_signal)),
        "mask_free_pixels": int(np.count_nonzero(mask_free)),
        "mask_bg_far_pixels": int(np.count_nonzero(mask_bg_far)),
        "mask_template_support_pixels": int(np.count_nonzero(mask_template_support)),
    }
    profiles = {
        "x_um": x_rel,
        "y_um": y_rel,
        "x_I": profile_x_I,
        "y_I": profile_y_I,
        "x_A": profile_x_A,
        "y_A": profile_y_A,
    }

    return RTADTarget(
        I_full=I_full,
        A_full=A_full,
        I_constraint=I_constraint,
        A_constraint=A_constraint,
        A_signal=A_signal,
        I_target=I_target,
        A_target=A_target,
        mask_flat=mask_flat,
        mask_edge_lock=mask_edge_lock,
        mask_signal=mask_signal,
        mask_template_support=mask_template_support,
        mask_guard_window=mask_guard_window,
        mask_bg_far=mask_bg_far,
        mask_edge=mask_edge,
        mask_support=mask_support,
        mask_bg=mask_bg,
        mask_free=mask_free,
        x_um=x_axis,
        y_um=y_axis,
        params=params,
        profiles=profiles,
        measured=measured,
    )
