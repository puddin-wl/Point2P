"""Raised-cosine RTAD rectangular target generation."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .metrics import width_at_level


@dataclass
class RTADTarget:
    """Container for target arrays, masks, axes, and geometry metadata."""

    I_target: np.ndarray
    A_target: np.ndarray
    mask_flat: np.ndarray
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
            "mask_edge": self.mask_edge,
            "mask_support": self.mask_support,
            "mask_bg": self.mask_bg,
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
    center_x_um: float = 0.0,
    center_y_um: float = 0.0,
    mode: str = "separable",
) -> RTADTarget:
    """Generate a separable raised-cosine RTAD rectangular target.

    The requested ``W50_um`` and ``H50_um`` are intensity 50 percent full sizes.
    The returned amplitude target is ``sqrt(I_target)`` for MRAF/GS projection.
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

    mode_l = mode.lower()
    if mode_l != "separable":
        raise ValueError("Only mode='separable' is implemented in this Python version.")

    Ix = raised_cosine_edge(absX, a0, a1)
    Iy = raised_cosine_edge(absY, b0, b1)
    I_target = np.clip(Ix * Iy, 0.0, 1.0).astype(np.float32)
    A_target = np.sqrt(I_target).astype(np.float32)

    mask_flat = (absX <= a0) & (absY <= b0)
    mask_support = I_target > 0
    mask_edge = mask_support & ~mask_flat
    mask_bg = I_target == 0
    mask_guard_support = (absX <= a2) & (absY <= b2)
    mask_free = mask_guard_support & ~mask_support

    ix0 = int(np.argmin(np.abs(x_axis - center_x_um)))
    iy0 = int(np.argmin(np.abs(y_axis - center_y_um)))
    profile_x_I = I_target[iy0, :].astype(np.float32)
    profile_y_I = I_target[:, ix0].astype(np.float32)
    profile_x_A = A_target[iy0, :].astype(np.float32)
    profile_y_A = A_target[:, ix0].astype(np.float32)

    x_rel = (x_axis - center_x_um).astype(np.float64)
    y_rel = (y_axis - center_y_um).astype(np.float64)
    measured = {
        "W50_x_um": width_at_level(x_rel, profile_x_I, level=0.5),
        "H50_y_um": width_at_level(y_rel, profile_y_I, level=0.5),
        "W90_x_um": width_at_level(x_rel, profile_x_I, level=0.9),
        "H90_y_um": width_at_level(y_rel, profile_y_I, level=0.9),
        "W10_x_um": width_at_level(x_rel, profile_x_I, level=0.1),
        "H10_y_um": width_at_level(y_rel, profile_y_I, level=0.1),
    }
    params = {
        "mode": mode_l,
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
        I_target=I_target,
        A_target=A_target,
        mask_flat=mask_flat,
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
