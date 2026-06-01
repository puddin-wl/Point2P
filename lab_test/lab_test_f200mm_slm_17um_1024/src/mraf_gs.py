"""WGS flat_local refinement for RTAD flat-top targets — 唯一活跃代码路径.
---
迭代流程 (method=wgs, strategy=flat_local):
  1. DOE → 焦面 (FFT 前向传播)
  2. mask_flat 内 2D 局部权重反馈: w ← w × (mean|E_flat| / |E|)^exponent
  3. MRAF 焦面投影: signal=加权目标, free×mraf_factor, bg×bg_factor
  4. 焦面 → DOE (IFFT 反传), 提取包裹相位
  5. 指标记录

以下代码路径已用 # 注释 (本管线不使用):
  - method='gs' / 'mraf' / 'mraf_then_wgs'
  - wgs_strategy='xy_then_x' (X-only WGS, 已弃用)
  - _update_x_wgs_weights() / _empty_x_weight_stats()
"""

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
    """一次精修运行的完整输出: 相位, 强度, 指标历史, WGS 权重, 各阶段快照."""

    phase_refined: np.ndarray
    reconstruction_intensity: np.ndarray
    initial_intensity: np.ndarray
    metrics_history: list[dict[str, Any]]
    final_metrics: dict[str, Any]
    backend_description: str
    phase_after_mraf: np.ndarray | None = None
    reconstruction_after_mraf: np.ndarray | None = None
    phase_after_wgs_xy: np.ndarray | None = None
    reconstruction_after_wgs_xy: np.ndarray | None = None
    mraf_metrics: dict[str, Any] | None = None
    wgs_xy_metrics: dict[str, Any] | None = None
    wgs_weights_final: np.ndarray | None = None
    wgs_weights_2d_final: np.ndarray | None = None
    wgs_x_weights_final: np.ndarray | None = None
    wgs_weight_stats_final: dict[str, float] | None = None
    warnings: list[str] | None = None


def _normalize_target_amp(target_amp: Any, signal_region: Any, xp: Any, dtype: Any) -> Any:
    """归一化目标振幅: 只在 signal 区非零, L2 范数=1."""
    target = xp.asarray(target_amp, dtype=dtype)
    target = xp.where(signal_region, target, dtype(0.0))
    return target / l2_norm(target, xp)


def _empty_weight_stats() -> dict[str, float]:
    """返回全 NaN 的 WGS 权重统计占位 (未使用 WGS 时)."""
    return {
        "wgs_weight_mean": float("nan"),
        "wgs_weight_std": float("nan"),
        "wgs_weight_min": float("nan"),
        "wgs_weight_max": float("nan"),
        "wgs_weight_p01": float("nan"),
        "wgs_weight_p99": float("nan"),
        "weights_2d_mean": float("nan"),
        "weights_2d_std": float("nan"),
        "weights_2d_min": float("nan"),
        "weights_2d_max": float("nan"),
        "w_x_mean": float("nan"),
        "w_x_std": float("nan"),
        "w_x_min": float("nan"),
        "w_x_max": float("nan"),
    }


def _array_stats(values: Any, backend: ArrayBackend, prefix: str) -> dict[str, float]:
    """计算数组的 mean/std/min/max 统计 (自动排除 NaN/Inf)."""
    vals = backend.to_numpy(values).astype(np.float64).reshape(-1)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {
            f"{prefix}_mean": float("nan"),
            f"{prefix}_std": float("nan"),
            f"{prefix}_min": float("nan"),
            f"{prefix}_max": float("nan"),
        }
    return {
        f"{prefix}_mean": float(np.mean(vals)),
        f"{prefix}_std": float(np.std(vals)),
        f"{prefix}_min": float(np.min(vals)),
        f"{prefix}_max": float(np.max(vals)),
    }


def _wgs_weight_stats(
    weights_2d: Any,
    update_region: Any,
    w_x: Any | None,
    valid_x: Any | None,
    backend: ArrayBackend,
) -> dict[str, float]:
    """汇总 2D flat-local 权重和可选 1D x-only 权重的统计信息."""
    vals_2d = backend.to_numpy(weights_2d[update_region]).astype(np.float64).reshape(-1)
    if vals_2d.size == 0:
        stats = _empty_weight_stats()
    else:
        stats = {
            "wgs_weight_mean": float(np.nanmean(vals_2d)),
            "wgs_weight_std": float(np.nanstd(vals_2d)),
            "wgs_weight_min": float(np.nanmin(vals_2d)),
            "wgs_weight_max": float(np.nanmax(vals_2d)),
            "wgs_weight_p01": float(np.nanpercentile(vals_2d, 1)),
            "wgs_weight_p99": float(np.nanpercentile(vals_2d, 99)),
            "weights_2d_mean": float(np.nanmean(vals_2d)),
            "weights_2d_std": float(np.nanstd(vals_2d)),
            "weights_2d_min": float(np.nanmin(vals_2d)),
            "weights_2d_max": float(np.nanmax(vals_2d)),
            "w_x_mean": float("nan"),
            "w_x_std": float("nan"),
            "w_x_min": float("nan"),
            "w_x_max": float("nan"),
        }
    if w_x is not None and valid_x is not None:
        stats.update(_array_stats(w_x[valid_x], backend, "w_x"))
    return stats


# ======== 以下两个函数仅供 xy_then_x 策略使用 (已弃用) ========

# def _empty_x_weight_stats() -> dict[str, float]:
#     """返回全 NaN 的 x-only WGS 统计占位."""
#     return {
#         "w_x_mean": float("nan"),
#         "w_x_std": float("nan"),
#         "w_x_min": float("nan"),
#         "w_x_max": float("nan"),
#     }


def _metrics_row(
    intensity_np: np.ndarray,
    masks: dict[str, np.ndarray],
    x_um: np.ndarray,
    y_um: np.ndarray,
    iteration: int,
    stage: str,
    wgs_strategy: str = "",
    wgs_substage: str = "",
    weight_stats: dict[str, float] | None = None,
) -> dict[str, Any]:
    """计算单轮指标并附加阶段/策略/权重统计信息."""
    row = compute_metrics(intensity_np, masks, x_um, y_um, iteration=iteration)
    row["stage"] = stage
    row["wgs_strategy"] = wgs_strategy
    row["wgs_substage"] = wgs_substage
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

    WGS 局部权重反馈公式:
      w_new = w × [mean(|E_flat|) / |E_map|]^exponent
    其中 |E_map| 是焦面振幅 (abs(farfield)), mean(|E_flat|) 是 mask_flat 内的平均振幅。
    只修改 update_region (mask_flat), 其他区域权重保持 1.0。
    更新后 clip 到 [clip_min, clip_max], 可选均值归一化。
    """
    eps = np.float32(1e-12)
    amp_flat = farfield_amp[update_region] #提取 mask_flat 区域的振幅
    if int(amp_flat.size) == 0:
        return weights, False, "WGS update skipped because mask_flat is empty."

    amp_mean_b = xp.mean(amp_flat) #计算 mask_flat 区域的平均振幅 (标量)
    amp_mean = float(np.asarray(backend.to_numpy(amp_mean_b)).reshape(())) #转换为 Python float
    if not np.isfinite(amp_mean) or amp_mean <= 0:
        return weights, False, f"WGS update skipped because flat mean amplitude is {amp_mean}."

    ratio = amp_mean_b / xp.maximum(amp_flat, eps)
    updated = weights[update_region] * xp.power(ratio, float(feedback_exponent)) #计算更新的权重
    updated = xp.clip(updated, float(clip_min), float(clip_max)) #将更新的权重裁剪到 [clip_min, clip_max]

    if normalize_weights:
        mean_w_b = xp.mean(updated)
        mean_w = float(np.asarray(backend.to_numpy(mean_w_b)).reshape(()))
        if np.isfinite(mean_w) and mean_w > 0:
            updated = updated / mean_w_b
        else:
            return weights, False, f"WGS weight normalization skipped because mean weight is {mean_w}."

    weights[update_region] = updated #将更新后的权重写回原数组
    return weights, True, ""


# def _update_x_wgs_weights(
#     w_x: Any,
#     farfield_amp: Any,
#     mask_flat: Any,
#     valid_x: Any,
#     flat_counts_x: Any,
#     xp: Any,
#     backend: ArrayBackend,
#     feedback_exponent: float,
#     clip_min: float,
#     clip_max: float,
#     normalize_weights: bool,
# ) -> Any:
#     """Update one-dimensional x weights from y-averaged flat-core amplitude.
#
#     For each x column, the current far-field amplitude is averaged over the
#     fixed RTAD ``mask_flat`` rows. The resulting 1D correction is broadcast
#     along y in the projection step, so this stage mainly changes the x profile
#     while freezing the y-direction local structure inherited from XY WGS.
#
#     已弃用 (xy_then_x). X-only WGS 会冻结 Y 方向结构, 损害均匀性。
#     """
#     eps = np.float32(1e-12)
#     valid_count = int(np.asarray(backend.to_numpy(xp.count_nonzero(valid_x))).reshape(()))
#     if valid_count == 0:
#         return w_x, False, "x-only WGS update skipped because no valid flat-core x columns exist."
#
#     amp_sum_x = xp.sum(xp.where(mask_flat, farfield_amp, 0), axis=0)
#     amp_x = amp_sum_x / xp.maximum(flat_counts_x, 1)
#     amp_valid = amp_x[valid_x]
#     amp_mean_b = xp.mean(amp_valid)
#     amp_mean = float(np.asarray(backend.to_numpy(amp_mean_b)).reshape(()))
#     if not np.isfinite(amp_mean) or amp_mean <= 0:
#         return w_x, False, f"x-only WGS update skipped because mean x amplitude is {amp_mean}."
#
#     ratio = amp_mean_b / xp.maximum(amp_valid, eps)
#     updated = w_x[valid_x] * xp.power(ratio, float(feedback_exponent))
#     updated = xp.clip(updated, float(clip_min), float(clip_max))
#
#     if normalize_weights:
#         mean_w_b = xp.mean(updated)
#         mean_w = float(np.asarray(backend.to_numpy(mean_w_b)).reshape(()))
#         if np.isfinite(mean_w) and mean_w > 0:
#             updated = updated / mean_w_b
#         else:
#             return w_x, False, f"x-only WGS normalization skipped because mean x weight is {mean_w}."
#
#     w_x[valid_x] = updated
#     w_x[~valid_x] = xp.asarray(1.0, dtype=w_x.dtype)
#     return w_x, True, ""


def _make_weighted_target(
    base_target: Any,
    mask_flat: Any,
    weights_2d: Any,
    w_x: Any,
    use_weights_2d: bool,
    use_x_weights: bool,
) -> Any:
    """构建 WGS 投影的有效信号目标振幅: base_target × weights_2d × w_x."""
    target_eff = base_target.copy()
    flat_target = base_target
    if use_weights_2d:
        flat_target = flat_target * weights_2d
    if use_x_weights:
        flat_target = flat_target * w_x[None, :]
    target_eff[mask_flat] = flat_target[mask_flat]
    return target_eff


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
    """Apply the farfield amplitude projection for GS or MRAF semantics.

    GS: 只有 signal 区替换为目标振幅, 其余置零
    MRAF: signal 区 = target_amp × exp(iφ), free 区 = mraf_factor × E, bg 区根据 bg_mode 处理
    """
    phase_ff = xp.angle(farfield)
    phase_factor = xp.exp(1j * phase_ff).astype(farfield.dtype, copy=False)
    signal = masks_b.get("mask_signal", masks_b["mask_support"])

    # --- GS 分支 (未使用, 实际只用 MRAF 投影) ---
    # if method == "gs":
    #     projected = xp.zeros_like(farfield)
    #     projected[signal] = target_amp_eff[signal] * phase_factor[signal]
    #     return projected

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
    """单次前向传播: input_amp × exp(iφ) → FFT → intensity."""
    xp = backend.xp
    field = input_amp * xp.exp(1j * phase).astype(backend.complex_dtype, copy=False)
    return intensity(forward_fft(field, xp), xp)


# ★ 主迭代引擎: GS/MRAF/WGS 统一精修循环
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
    wgs_strategy: str = "flat_local",
    wgs_xy_iters: int = 20,
    wgs_xonly_iters: int = 30,
    wgs_feedback_exponent: float | None = None,
    wgs_xy_feedback_exponent: float | None = None,
    wgs_x_feedback_exponent: float = 0.45,
    wgs_update_every: int = 5,
    wgs_xy_update_every: int | None = None,
    wgs_x_update_every: int = 5,
    wgs_update_mask: str = "flat",
    wgs_weight_min: float | None = None,
    wgs_weight_max: float | None = None,
    wgs_xy_weight_min: float | None = None,
    wgs_xy_weight_max: float | None = None,
    wgs_x_weight_min: float = 0.5,
    wgs_x_weight_max: float = 2.5,
    wgs_normalize_weights: bool = True,
    wgs_x_normalize: bool = True,
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

    迭代循环的核心逻辑:
    1. 前向传播: field → farfield = FFT(input_amp × exp(iφ))
    2. WGS 权重更新 (每隔 update_every 轮): 根据焦面振幅更新 mask_flat 内的权重
    3. 焦面投影: 将 farfield 投影到目标振幅约束 (GS 或 MRAF 语义)
    4. 反传: nearfield = IFFT(projected)
    5. 提取相位: φ = angle(nearfield), wrap 到 [0, 2π)
    6. 指标记录: 每隔 metrics_interval 轮计算 size50, RMS, 效率等
    """
    method_l = method.lower()
    if method_l not in {"gs", "mraf", "wgs", "mraf_then_wgs"}:
        raise ValueError("method must be one of: gs, mraf, wgs, mraf_then_wgs.")
    if phase0.shape != target_amp.shape:
        raise ValueError(f"phase0 shape {phase0.shape} does not match target shape {target_amp.shape}.")
    if wgs_feedback.lower() != "amplitude":
        raise ValueError("Only wgs_feedback='amplitude' is implemented.")
    wgs_strategy_l = wgs_strategy.lower()
    if wgs_strategy_l not in {"flat_local", "xy_then_x"}:
        raise ValueError("wgs_strategy must be 'flat_local' or 'xy_then_x'.")
    if int(wgs_update_every) <= 0 or int(wgs_x_update_every) <= 0:
        raise ValueError("WGS update intervals must be positive integers.")

    if wgs_feedback_exponent is None:
        wgs_feedback_exponent = feedback_exponent
    if wgs_xy_feedback_exponent is None:
        wgs_xy_feedback_exponent = wgs_feedback_exponent
    if wgs_xy_update_every is None:
        wgs_xy_update_every = wgs_update_every
    if wgs_weight_min is None:
        wgs_weight_min = wgs_clip_min
    if wgs_weight_max is None:
        wgs_weight_max = wgs_clip_max
    if wgs_xy_weight_min is None:
        wgs_xy_weight_min = wgs_weight_min
    if wgs_xy_weight_max is None:
        wgs_xy_weight_max = wgs_weight_max
    if wgs_weight_min <= 0 or wgs_weight_max <= 0 or wgs_weight_min > wgs_weight_max:
        raise ValueError("WGS weight clip limits must satisfy 0 < min <= max.")
    if wgs_xy_weight_min <= 0 or wgs_xy_weight_max <= 0 or wgs_xy_weight_min > wgs_xy_weight_max:
        raise ValueError("WGS XY weight clip limits must satisfy 0 < min <= max.")
    if wgs_x_weight_min <= 0 or wgs_x_weight_max <= 0 or wgs_x_weight_min > wgs_x_weight_max:
        raise ValueError("WGS x-only weight clip limits must satisfy 0 < min <= max.")

    xp = backend.xp
    dtype = backend.float_dtype
    phase = backend.to_backend(np.mod(phase0.astype(np.float32), 2.0 * np.pi), dtype=dtype)
    input_amp_b = backend.to_backend(input_amp, dtype=dtype)
    input_amp_b = input_amp_b / l2_norm(input_amp_b, xp)

    masks_b = {k: backend.to_backend(v, dtype=bool) for k, v in masks.items()}
    signal = masks_b.get("mask_signal", masks_b["mask_support"])
    base_target = _normalize_target_amp(target_amp, signal, xp, dtype)
    wgs_weights_2d = xp.ones_like(base_target, dtype=dtype)
    w_x = xp.ones(base_target.shape[1], dtype=dtype)

    update_mask_name = (wgs_update_mask or wgs_update_region).lower()
    if update_mask_name == "flat":
        update_region = masks_b["mask_flat"]
    else:
        raise ValueError("Only wgs_update_mask='flat' is supported in this first WGS version.")
    valid_x = xp.any(update_region, axis=0)
    flat_counts_x = xp.sum(update_region, axis=0).astype(dtype, copy=False)

    # ======== 迭代计数 (只用 WGS flat_local, 其他分支已注释) ========
    mraf_count = 0
    xy_count = 0
    xonly_count = 0
    wgs_count = int(num_iters)
    total_iters = int(num_iters)
    #
    # # --- 已弃用: mraf_then_wgs 两段式 ---
    # if method_l == "mraf_then_wgs":
    #     mraf_count = int(mraf_iters if mraf_iters is not None else (wgs_after_iters if wgs_after_iters > 0 else num_iters))
    #     if wgs_strategy_l == "xy_then_x":
    #         xy_count = int(wgs_xy_iters)
    #         xonly_count = int(wgs_xonly_iters)
    #         wgs_count = xy_count + xonly_count
    #     else:
    #         xy_count = 0; xonly_count = 0
    #         wgs_count = int(wgs_iters)
    #     total_iters = mraf_count + wgs_count
    #     ...
    # else:
    #     mraf_count = 0
    #     # --- 已弃用: xy_then_x 策略 ---
    #     # if method_l == "wgs" and wgs_strategy_l == "xy_then_x":
    #     #     xy_count = int(wgs_xy_iters)
    #     #     xonly_count = int(wgs_xonly_iters)
    #     #     wgs_count = xy_count + xonly_count
    #     #     total_iters = wgs_count
    #     # else:
    #         xy_count = 0; xonly_count = 0
    #         wgs_count = int(num_iters) if method_l == "wgs" else 0
    #         total_iters = int(num_iters)
    #     ...
    if total_iters <= 0:
        raise ValueError("num_iters must be positive.")

    run_warnings: list[str] = []
    # # --- 已弃用: xy_then_x 警告 ---
    # if wgs_strategy_l == "xy_then_x" and int(wgs_iters) != int(wgs_count):
    #     run_warnings.append(...)

    initial_I = _evaluate_phase(phase, input_amp_b, backend)
    initial_I_np = backend.to_numpy(initial_I).astype(np.float32)
    metrics_history: list[dict[str, Any]] = [
        _metrics_row(initial_I_np, masks, x_um, y_um, iteration=0, stage="initial", wgs_strategy=wgs_strategy_l)
    ]

    phase_after_mraf = None
    reconstruction_after_mraf = None
    phase_after_wgs_xy = None
    reconstruction_after_wgs_xy = None
    mraf_metrics = None
    wgs_xy_metrics = None

    iterator = range(1, total_iters + 1)
    if show_progress:
        iterator = tqdm(iterator, desc=f"{method_l} refinement", unit="iter")

    for it in iterator:
        # ======== WGS flat_local 活跃路径 ========
        # 以下被注释掉的代码是 method='mraf_then_wgs' / wgs_strategy='xy_then_x' 的
        # 阶段切换逻辑。当前只用 WGS flat_local:
        #   stage = "wgs", do_wgs_xy = True, do_wgs_xonly = False

        # # --- 已弃用: mraf_then_wgs 阶段切换 ---
        # if method_l == "mraf_then_wgs":
        #     if it <= mraf_count: stage = "mraf" ...
        #     else: stage = "wgs"; ...
        #     if wgs_strategy_l == "xy_then_x": ...
        #     else: do_wgs_xy = True
        # # elif method_l == "wgs":
        #     # --- 已弃用: xy_then_x 子分支 ---
        #     # if wgs_strategy_l == "xy_then_x": ...
        #     # else: ...
        # # --- 已弃用: pure GS/MRAF ---
        # # else: stage = method_l

        stage = "wgs"
        wgs_substage = "wgs_flat_local"
        do_wgs_xy = True
        do_wgs_xonly = False

        # 1. 前向传播
        field = input_amp_b * xp.exp(1j * phase).astype(backend.complex_dtype, copy=False)
        farfield = forward_fft(field, xp)
        farfield_amp = xp.abs(farfield)

        # 2. WGS 权重更新 (每隔 update_every 轮, 活跃路径: 2D flat_local)
        weights_updated = False
        # -- 活跃路径: flat_local WGS --
        # if do_wgs_xy:  (always True)
        update_every = int(wgs_update_every)
        feedback_exp = float(wgs_feedback_exponent)
        clip_min = float(wgs_weight_min)
        clip_max = float(wgs_weight_max)
        if it % update_every == 0:
            wgs_weights_2d, weights_updated, warning_msg = _update_flat_wgs_weights(
                wgs_weights_2d, farfield_amp, update_region, xp, backend,
                feedback_exponent=feedback_exp, clip_min=clip_min, clip_max=clip_max,
                normalize_weights=bool(wgs_normalize_weights),
            )
            if warning_msg:
                run_warnings.append(f"iteration {it}: {warning_msg}")
                warnings.warn(warning_msg, RuntimeWarning, stacklevel=2)
        # # --- 已弃用: x-only WGS 更新 ---
        # elif do_wgs_xonly and ...:
        #     w_x, weights_updated, warning_msg = _update_x_wgs_weights(...)

        # 3. 构建加权目标 (2D weights only, no x-weights)
        # -- 活跃路径: XY weights ---
        # if do_wgs_xy:  (always True)
        target_eff = _make_weighted_target(
            base_target, update_region, wgs_weights_2d, w_x,
            use_weights_2d=True, use_x_weights=False,
        )
        # # --- 已弃用: x-only 加权 ---
        # elif do_wgs_xonly: target_eff = _make_weighted_target(..., use_x_weights=True)
        # # --- 已弃用: 无权重 (纯 GS/MRAF) ---
        # else: target_eff = base_target

        # 4. MRAF 焦面投影 + 反传 + 相位提取
        # 注意: project_method 在 WGS 路径下始终为 "mraf"
        # # project_method = "gs" if method_l == "gs" else "mraf"
        projected = _project_farfield(
            farfield, target_eff, masks_b, xp,
            method="mraf", mraf_factor=float(mraf_factor),
            bg_mode=bg_mode, bg_factor=float(bg_factor),
        )
        nearfield = backward_fft(projected, xp)
        phase = xp.mod(xp.angle(nearfield), 2.0 * np.pi).astype(dtype, copy=False)

        # 5. 指标记录
        # # --- 已弃用: MRAF/XY 阶段快照 ---
        # need_mraf_snapshot = method_l == "mraf_then_wgs" and ...
        # need_xy_snapshot = wgs_strategy_l == "xy_then_x" and ...
        need_log = (
            (metrics_interval > 0 and it % metrics_interval == 0)
            or it == total_iters
            or weights_updated
        )
        if need_log:
            I_now = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
            # # --- 已弃用: MRAF 快照保存 ---
            # if need_mraf_snapshot and I_now is not None:
            #     phase_after_mraf = ...
            # # --- 已弃用: XY 快照保存 ---
            # if need_xy_snapshot and I_now is not None:
            #     phase_after_wgs_xy = ...
            metrics_history.append(
                _metrics_row(
                    I_now, masks, x_um, y_um, iteration=it,
                    stage=stage, wgs_strategy=wgs_strategy_l, wgs_substage=wgs_substage,
                    weight_stats=_wgs_weight_stats(wgs_weights_2d, update_region, w_x, valid_x, backend),
                )
            )

    refined_I = backend.to_numpy(_evaluate_phase(phase, input_amp_b, backend)).astype(np.float32)
    phase_refined = backend.to_numpy(phase).astype(np.float32)
    # WGS flat_local: stage/substage 固定
    final_stage = "wgs"
    final_substage = "wgs_flat_local"
    # # --- 已弃用: 多阶段/xy_then_x final_substage 逻辑 ---
    # final_stage = "wgs" if method_l in {"wgs", "mraf_then_wgs"} and wgs_count > 0 else method_l
    # final_substage = "wgs_xonly" if ... else ("wgs_xy" if ... else final_stage)
    final_weight_stats = _wgs_weight_stats(wgs_weights_2d, update_region, w_x, valid_x, backend)
    final_metrics = _metrics_row(
        refined_I,
        masks,
        x_um,
        y_um,
        iteration=total_iters,
        stage=final_stage,
        wgs_strategy=wgs_strategy_l,
        wgs_substage=final_substage,
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
        phase_after_wgs_xy=phase_after_wgs_xy,
        reconstruction_after_wgs_xy=reconstruction_after_wgs_xy,
        mraf_metrics=mraf_metrics,
        wgs_xy_metrics=wgs_xy_metrics,
        wgs_weights_final=backend.to_numpy(wgs_weights_2d).astype(np.float32),
        wgs_weights_2d_final=backend.to_numpy(wgs_weights_2d).astype(np.float32),
        wgs_x_weights_final=backend.to_numpy(w_x).astype(np.float32),
        wgs_weight_stats_final=final_weight_stats,
        warnings=run_warnings,
    )
