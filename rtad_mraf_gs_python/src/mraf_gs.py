"""Simplified GS/MRAF/WGS refinement loop for RTAD flat-top targets."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from tqdm.auto import tqdm

from .backend import ArrayBackend
from .metrics import compute_metrics
from .propagation import backward_fft, forward_fft, intensity, l2_norm


@dataclass
class RefinementResult:
    """Outputs from one refinement run."""

    phase_refined: np.ndarray
    reconstruction_intensity: np.ndarray
    initial_intensity: np.ndarray
    metrics_history: list[dict[str, Any]]
    final_metrics: dict[str, Any]
    backend_description: str


def _normalize_target_amp(target_amp: Any, signal_region: Any, xp: Any, dtype: Any) -> Any:
    """Normalize target amplitude over finite signal/support pixels."""
    target = xp.asarray(target_amp, dtype=dtype)
    target = xp.where(signal_region, target, dtype(0.0))
    return target / l2_norm(target, xp)


def _update_wgs_weights(
    weights: Any,
    base_target: Any,
    farfield_amp: Any,
    signal_region: Any,
    update_region: Any,
    xp: Any,
    feedback_exponent: float,
    clip_min: float,
    clip_max: float,
) -> Any:
    """Apply a Leonardo-style WGS weight update.

    This mirrors the main idea in slmsuite's ``_update_weights_generic``:
    increase target weight where measured/computed amplitude is weak and
    decrease it where amplitude is strong. This local version updates only the
    configured RTAD region, then L2-normalizes the signal weights.
    """
    eps = np.float32(1e-12)
    feedback_scaled = farfield_amp * (l2_norm(base_target[signal_region], xp) / l2_norm(farfield_amp[signal_region], xp))
    correction = xp.power(base_target[update_region] / xp.maximum(feedback_scaled[update_region], eps), feedback_exponent)
    correction = xp.clip(correction, clip_min, clip_max)
    weights[update_region] *= correction
    weights *= l2_norm(base_target[signal_region], xp) / l2_norm(weights[signal_region], xp)
    return weights


def _project_farfield(
    farfield: Any,
    weights: Any,
    masks_b: dict[str, Any],
    xp: Any,
    method: str,
    mraf_factor: float,
    bg_mode: str,
    bg_factor: float,
) -> Any:
    """Apply the farfield amplitude constraint for GS, MRAF, or WGS."""
    phase_ff = xp.angle(farfield)
    phase_factor = xp.exp(1j * phase_ff).astype(farfield.dtype, copy=False)
    signal = masks_b.get("mask_signal", masks_b["mask_support"])

    if method in {"gs", "wgs"}:
        projected = xp.zeros_like(farfield)
        projected[signal] = weights[signal] * phase_factor[signal]
        return projected

    projected = farfield.copy()
    projected[signal] = weights[signal] * phase_factor[signal]

    free = masks_b["mask_free"]
    projected[free] *= mraf_factor

    bg = masks_b.get("mask_bg_far", masks_b["mask_bg"] & ~masks_b["mask_free"])
    if bg_mode == "keep":
        pass
    elif bg_mode == "attenuate":
        projected[bg] *= bg_factor
    elif bg_mode == "zero":
        projected[bg] = 0
    else:
        raise ValueError(f"Unknown bg_mode '{bg_mode}'. Use keep, attenuate, or zero.")
    return projected


def _evaluate_phase(phase: Any, input_amp: Any, backend: ArrayBackend) -> Any:
    """Forward propagate one phase pattern and return intensity."""
    xp = backend.xp
    field = input_amp * xp.exp(1j * phase).astype(backend.complex_dtype, copy=False)
    return intensity(forward_fft(field, xp), xp)


def run_refinement(
    phase0: np.ndarray,
    input_amp: Any,
    target_amp: np.ndarray,
    masks: dict[str, np.ndarray],
    x_um: np.ndarray,
    y_um: np.ndarray,
    backend: ArrayBackend,
    method: str = "mraf",
    num_iters: int = 200,
    mraf_factor: float = 0.4,
    wgs_after_iters: int = 0,
    feedback_exponent: float = 0.7,
    wgs_update_region: str = "flat",
    bg_mode: str = "attenuate",
    bg_factor: float = 0.05,
    wgs_clip_min: float = 0.5,
    wgs_clip_max: float = 2.0,
    metrics_interval: int = 10,
    show_progress: bool = True,
) -> RefinementResult:
    """Run GS/MRAF/WGS refinement from an initial phase.

    MRAF formula used here:

    - signal pixels: ``E' = W * exp(i * angle(E))``
    - free/noise pixels: ``E' = mraf_factor * E``
    - far background: keep, attenuate by ``bg_factor``, or zero by option

    This follows slmsuite's convention that ``mraf_factor=0`` fully attenuates
    a noise/free region while ``mraf_factor=1`` leaves it unchanged. The main
    project difference is the truncated RTAD target: only ``mask_signal``
    follows the target amplitude. The low-intensity full-template tail is
    released into ``mask_free`` and is not forced to match the mathematical
    raised-cosine curve.
    """
    method_l = method.lower()
    if method_l not in {"gs", "mraf", "wgs", "mraf_then_wgs"}:
        raise ValueError("method must be one of: gs, mraf, wgs, mraf_then_wgs.")
    if phase0.shape != target_amp.shape:
        raise ValueError(f"phase0 shape {phase0.shape} does not match target shape {target_amp.shape}.")

    xp = backend.xp
    dtype = backend.float_dtype
    phase = backend.to_backend(np.mod(phase0.astype(np.float32), 2.0 * np.pi), dtype=dtype)
    input_amp_b = backend.to_backend(input_amp, dtype=dtype)
    input_amp_b = input_amp_b / l2_norm(input_amp_b, xp)

    masks_b = {k: backend.to_backend(v, dtype=bool) for k, v in masks.items()}
    signal = masks_b.get("mask_signal", masks_b["mask_support"])
    base_target = _normalize_target_amp(target_amp, signal, xp, dtype)
    weights = base_target.copy()

    if wgs_update_region == "flat":
        update_region = masks_b["mask_flat"]
    elif wgs_update_region == "support":
        update_region = masks_b["mask_support"]
    else:
        raise ValueError("wgs_update_region must be 'flat' or 'support'.")

    initial_I = _evaluate_phase(phase, input_amp_b, backend)
    initial_I_np = backend.to_numpy(initial_I).astype(np.float32)
    metrics_history: list[dict[str, Any]] = [
        compute_metrics(initial_I_np, masks, x_um, y_um, iteration=0)
    ]

    iterator = range(1, int(num_iters) + 1)
    if show_progress:
        iterator = tqdm(iterator, desc=f"{method_l} refinement", unit="iter")

    for it in iterator:
        field = input_amp_b * xp.exp(1j * phase).astype(backend.complex_dtype, copy=False)
        farfield = forward_fft(field, xp)
        farfield_amp = xp.abs(farfield)

        do_wgs = method_l == "wgs" or (method_l == "mraf_then_wgs" and it >= int(wgs_after_iters))
        if do_wgs and it > 1:
            weights = _update_wgs_weights(
                weights,
                base_target,
                farfield_amp,
                signal,
                update_region,
                xp,
                feedback_exponent=float(feedback_exponent),
                clip_min=float(wgs_clip_min),
                clip_max=float(wgs_clip_max),
            )

        project_method = "wgs" if method_l == "wgs" else ("mraf" if method_l in {"mraf", "mraf_then_wgs"} else "gs")
        projected = _project_farfield(
            farfield,
            weights,
            masks_b,
            xp,
            method=project_method,
            mraf_factor=float(mraf_factor),
            bg_mode=bg_mode,
            bg_factor=float(bg_factor),
        )
        nearfield = backward_fft(projected, xp)
        phase = xp.mod(xp.angle(nearfield), 2.0 * np.pi).astype(dtype, copy=False)

        if metrics_interval > 0 and (it % metrics_interval == 0 or it == int(num_iters)):
            I_now = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
            metrics_history.append(compute_metrics(I_now, masks, x_um, y_um, iteration=it))

    refined_I = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
    phase_refined = backend.to_numpy(phase).astype(np.float32)
    final_metrics = compute_metrics(refined_I, masks, x_um, y_um, iteration=int(num_iters))
    if metrics_history[-1]["iteration"] != int(num_iters):
        metrics_history.append(final_metrics)

    return RefinementResult(
        phase_refined=phase_refined,
        reconstruction_intensity=refined_I,
        initial_intensity=initial_I_np,
        metrics_history=metrics_history,
        final_metrics=final_metrics,
        backend_description=backend.describe(),
    )
