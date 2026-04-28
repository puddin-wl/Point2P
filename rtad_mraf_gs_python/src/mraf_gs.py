"""Simplified GS/MRAF/WGS refinement loop for RTAD flat-top targets."""

from __future__ import annotations

import warnings
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
    phase_after_mraf: np.ndarray | None = None
    reconstruction_after_mraf: np.ndarray | None = None
    mraf_metrics: dict[str, Any] | None = None
    wgs_weights_final: np.ndarray | None = None
    wgs_weight_stats_final: dict[str, float] | None = None
    warnings: list[str] | None = None


def _normalize_target_amp(target_amp: Any, signal_region: Any, xp: Any, dtype: Any) -> Any:
    """Normalize target amplitude over finite signal/support pixels."""
    target = xp.asarray(target_amp, dtype=dtype)
    target = xp.where(signal_region, target, dtype(0.0))
    return target / l2_norm(target, xp)


def _empty_weight_stats() -> dict[str, float]:
    """Return NaN placeholders for WGS weight metrics."""
    return {
        "wgs_weight_mean": float("nan"),
        "wgs_weight_std": float("nan"),
        "wgs_weight_min": float("nan"),
        "wgs_weight_max": float("nan"),
        "wgs_weight_p01": float("nan"),
        "wgs_weight_p99": float("nan"),
    }


def _wgs_weight_stats(weights: Any, update_region: Any, backend: ArrayBackend) -> dict[str, float]:
    """Summarize WGS weights inside the update region."""
    vals = backend.to_numpy(weights[update_region]).astype(np.float64).reshape(-1)
    if vals.size == 0:
        return _empty_weight_stats()
    return {
        "wgs_weight_mean": float(np.nanmean(vals)),
        "wgs_weight_std": float(np.nanstd(vals)),
        "wgs_weight_min": float(np.nanmin(vals)),
        "wgs_weight_max": float(np.nanmax(vals)),
        "wgs_weight_p01": float(np.nanpercentile(vals, 1)),
        "wgs_weight_p99": float(np.nanpercentile(vals, 99)),
    }


def _metrics_row(
    intensity_np: np.ndarray,
    masks: dict[str, np.ndarray],
    x_um: np.ndarray,
    y_um: np.ndarray,
    iteration: int,
    stage: str,
    weight_stats: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Compute metrics and attach the refinement stage and WGS stats."""
    row = compute_metrics(intensity_np, masks, x_um, y_um, iteration=iteration)
    row["stage"] = stage
    row.update(weight_stats if weight_stats is not None else _empty_weight_stats())
    return row


def _update_flat_wgs_weights(
    weights: Any,
    farfield_amp: Any,
    update_region: Any,
    xp: Any,
    backend: ArrayBackend,
    feedback_exponent: float,
    clip_min: float,
    clip_max: float,
    normalize_weights: bool,
) -> Any:
    """Apply an amplitude WGS update only inside the RTAD flat core.

    The update follows the requested flat-only feedback formula:

    ``weights <- weights * (mean(abs(E_flat)) / abs(E))**feedback_exponent``

    Only ``update_region`` is modified. For the current project this region is
    ``mask_flat``; edge, free, and far-background pixels keep weight 1.
    """
    eps = np.float32(1e-12)
    amp_flat = farfield_amp[update_region]
    if int(amp_flat.size) == 0:
        return weights, False, "WGS update skipped because mask_flat is empty."

    amp_mean_b = xp.mean(amp_flat)
    amp_mean = float(np.asarray(backend.to_numpy(amp_mean_b)).reshape(()))
    if not np.isfinite(amp_mean) or amp_mean <= 0:
        return weights, False, f"WGS update skipped because flat mean amplitude is {amp_mean}."

    ratio = amp_mean_b / xp.maximum(amp_flat, eps)
    updated = weights[update_region] * xp.power(ratio, float(feedback_exponent))
    updated = xp.clip(updated, float(clip_min), float(clip_max))

    if normalize_weights:
        mean_w_b = xp.mean(updated)
        mean_w = float(np.asarray(backend.to_numpy(mean_w_b)).reshape(()))
        if np.isfinite(mean_w) and mean_w > 0:
            updated = updated / mean_w_b
        else:
            return weights, False, f"WGS weight normalization skipped because mean weight is {mean_w}."

    weights[update_region] = updated
    return weights, True, ""


def _project_farfield(
    farfield: Any,
    target_amp_eff: Any,
    masks_b: dict[str, Any],
    xp: Any,
    method: str,
    mraf_factor: float,
    bg_mode: str,
    bg_factor: float,
) -> Any:
    """Apply the farfield amplitude projection for GS or MRAF semantics."""
    phase_ff = xp.angle(farfield)
    phase_factor = xp.exp(1j * phase_ff).astype(farfield.dtype, copy=False)
    signal = masks_b.get("mask_signal", masks_b["mask_support"])

    if method == "gs":
        projected = xp.zeros_like(farfield)
        projected[signal] = target_amp_eff[signal] * phase_factor[signal]
        return projected

    projected = farfield.copy()
    projected[signal] = target_amp_eff[signal] * phase_factor[signal]

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
    mraf_iters: int | None = None,
    wgs_iters: int = 0,
    mraf_factor: float = 0.4,
    wgs_after_iters: int = 0,
    feedback_exponent: float = 0.7,
    wgs_feedback: str = "amplitude",
    wgs_feedback_exponent: float | None = None,
    wgs_update_every: int = 5,
    wgs_update_mask: str = "flat",
    wgs_weight_min: float | None = None,
    wgs_weight_max: float | None = None,
    wgs_normalize_weights: bool = True,
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
    if wgs_feedback.lower() != "amplitude":
        raise ValueError("Only wgs_feedback='amplitude' is implemented.")
    if int(wgs_update_every) <= 0:
        raise ValueError("wgs_update_every must be a positive integer.")

    if wgs_feedback_exponent is None:
        wgs_feedback_exponent = feedback_exponent
    if wgs_weight_min is None:
        wgs_weight_min = wgs_clip_min
    if wgs_weight_max is None:
        wgs_weight_max = wgs_clip_max
    if wgs_weight_min <= 0 or wgs_weight_max <= 0 or wgs_weight_min > wgs_weight_max:
        raise ValueError("WGS weight clip limits must satisfy 0 < min <= max.")

    xp = backend.xp
    dtype = backend.float_dtype
    phase = backend.to_backend(np.mod(phase0.astype(np.float32), 2.0 * np.pi), dtype=dtype)
    input_amp_b = backend.to_backend(input_amp, dtype=dtype)
    input_amp_b = input_amp_b / l2_norm(input_amp_b, xp)

    masks_b = {k: backend.to_backend(v, dtype=bool) for k, v in masks.items()}
    signal = masks_b.get("mask_signal", masks_b["mask_support"])
    base_target = _normalize_target_amp(target_amp, signal, xp, dtype)
    wgs_weights = xp.ones_like(base_target, dtype=dtype)

    update_mask_name = (wgs_update_mask or wgs_update_region).lower()
    if update_mask_name == "flat":
        update_region = masks_b["mask_flat"]
    else:
        raise ValueError("Only wgs_update_mask='flat' is supported in this first WGS version.")

    if method_l == "mraf_then_wgs":
        mraf_count = int(mraf_iters if mraf_iters is not None else (wgs_after_iters if wgs_after_iters > 0 else num_iters))
        wgs_count = int(wgs_iters)
        total_iters = mraf_count + wgs_count
        if mraf_count < 0 or wgs_count < 0 or total_iters <= 0:
            raise ValueError("mraf_then_wgs requires non-negative mraf_iters/wgs_iters and at least one total iteration.")
    else:
        mraf_count = 0
        wgs_count = int(num_iters) if method_l == "wgs" else 0
        total_iters = int(num_iters)
        if total_iters <= 0:
            raise ValueError("num_iters must be positive.")

    initial_I = _evaluate_phase(phase, input_amp_b, backend)
    initial_I_np = backend.to_numpy(initial_I).astype(np.float32)
    metrics_history: list[dict[str, Any]] = [
        _metrics_row(initial_I_np, masks, x_um, y_um, iteration=0, stage="initial")
    ]

    phase_after_mraf = None
    reconstruction_after_mraf = None
    mraf_metrics = None
    run_warnings: list[str] = []

    iterator = range(1, total_iters + 1)
    if show_progress:
        iterator = tqdm(iterator, desc=f"{method_l} refinement", unit="iter")

    for it in iterator:
        if method_l == "mraf_then_wgs":
            stage = "mraf" if it <= mraf_count else "wgs"
            local_wgs_it = it - mraf_count
            do_wgs = stage == "wgs"
        elif method_l == "wgs":
            stage = "wgs"
            local_wgs_it = it
            do_wgs = True
        else:
            stage = method_l
            local_wgs_it = 0
            do_wgs = False

        field = input_amp_b * xp.exp(1j * phase).astype(backend.complex_dtype, copy=False)
        farfield = forward_fft(field, xp)
        farfield_amp = xp.abs(farfield)

        weights_updated = False
        if do_wgs and local_wgs_it > 0 and (local_wgs_it % int(wgs_update_every) == 0):
            wgs_weights, weights_updated, warning_msg = _update_flat_wgs_weights(
                wgs_weights,
                farfield_amp,
                update_region,
                xp,
                backend,
                feedback_exponent=float(wgs_feedback_exponent),
                clip_min=float(wgs_weight_min),
                clip_max=float(wgs_weight_max),
                normalize_weights=bool(wgs_normalize_weights),
            )
            if warning_msg:
                run_warnings.append(f"iteration {it}: {warning_msg}")
                warnings.warn(warning_msg, RuntimeWarning, stacklevel=2)

        if do_wgs:
            target_eff = base_target.copy()
            target_eff[update_region] = base_target[update_region] * wgs_weights[update_region]
        else:
            target_eff = base_target

        project_method = "gs" if method_l == "gs" else "mraf"
        projected = _project_farfield(
            farfield,
            target_eff,
            masks_b,
            xp,
            method=project_method,
            mraf_factor=float(mraf_factor),
            bg_mode=bg_mode,
            bg_factor=float(bg_factor),
        )
        nearfield = backward_fft(projected, xp)
        phase = xp.mod(xp.angle(nearfield), 2.0 * np.pi).astype(dtype, copy=False)

        need_mraf_snapshot = method_l == "mraf_then_wgs" and it == mraf_count
        need_log = (
            (metrics_interval > 0 and it % metrics_interval == 0)
            or it == total_iters
            or need_mraf_snapshot
            or weights_updated
        )
        I_now = None
        if need_mraf_snapshot or need_log:
            I_now = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
        if need_mraf_snapshot and I_now is not None:
            phase_after_mraf = backend.to_numpy(phase).astype(np.float32)
            reconstruction_after_mraf = I_now
            mraf_metrics = _metrics_row(
                I_now,
                masks,
                x_um,
                y_um,
                iteration=it,
                stage="mraf",
                weight_stats=_wgs_weight_stats(wgs_weights, update_region, backend),
            )
        if need_log and I_now is not None:
            metrics_history.append(
                _metrics_row(
                    I_now,
                    masks,
                    x_um,
                    y_um,
                    iteration=it,
                    stage=stage,
                    weight_stats=_wgs_weight_stats(wgs_weights, update_region, backend),
                )
            )

    refined_I = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
    phase_refined = backend.to_numpy(phase).astype(np.float32)
    final_stage = "wgs" if method_l in {"wgs", "mraf_then_wgs"} and wgs_count > 0 else method_l
    final_weight_stats = _wgs_weight_stats(wgs_weights, update_region, backend)
    final_metrics = _metrics_row(
        refined_I,
        masks,
        x_um,
        y_um,
        iteration=total_iters,
        stage=final_stage,
        weight_stats=final_weight_stats,
    )
    if metrics_history[-1]["iteration"] != total_iters:
        metrics_history.append(final_metrics)

    return RefinementResult(
        phase_refined=phase_refined,
        reconstruction_intensity=refined_I,
        initial_intensity=initial_I_np,
        metrics_history=metrics_history,
        final_metrics=final_metrics,
        backend_description=backend.describe(),
        phase_after_mraf=phase_after_mraf,
        reconstruction_after_mraf=reconstruction_after_mraf,
        mraf_metrics=mraf_metrics,
        wgs_weights_final=backend.to_numpy(wgs_weights).astype(np.float32),
        wgs_weight_stats_final=final_weight_stats,
        warnings=run_warnings,
    )
