"""Fourier-lens propagation utilities.
---
傅里叶透镜正反传播 + 输入场构建。
矩阵约定: axis 0/row = y, axis 1/col = x (与 MATLAB meshgrid 一致)。
正交 FFT (norm="ortho") 保证离散总功率在变换前后守恒。
"""

from __future__ import annotations

from typing import Any


def forward_fft(field_in: Any, xp: Any) -> Any:
    """DOE 平面 → 焦平面: FFT 前向传播, 正交归一保证功率守恒."""
    # fftshift(fft2(ifftshift(field))), norm="ortho"
    return xp.fft.fftshift(xp.fft.fft2(xp.fft.ifftshift(field_in), norm="ortho"))


def backward_fft(field_out: Any, xp: Any) -> Any:
    """焦平面 → DOE 平面: 逆 FFT 反传, 正交归一保证功率守恒."""
    return xp.fft.fftshift(xp.fft.ifft2(xp.fft.ifftshift(field_out), norm="ortho"))


def intensity(field: Any, xp: Any) -> Any:
    """计算光强 I = |field|²."""
    amp = xp.abs(field)
    return amp * amp


def l2_norm(array: Any, xp: Any, eps: float = 1e-20) -> Any:
    """计算 L2 范数 sqrt(Σ|a|²), 含下限保护防止除零."""
    value = xp.sqrt(xp.sum(xp.abs(array) ** 2))
    return xp.maximum(value, eps)


def normalize_power(array: Any, xp: Any, target_power: float = 1.0, eps: float = 1e-20) -> Any:
    """将振幅/场归一化到指定总功率 (默认 1.0)."""
    norm = l2_norm(array, xp, eps=eps)
    return array * (target_power ** 0.5 / norm)


def normalize_mean_in_mask(values: Any, mask: Any, xp: Any, target_mean: float = 1.0, eps: float = 1e-20) -> Any:
    """将 mask 内均值缩放到 target_mean."""
    mean = xp.mean(values[mask])
    return values * (target_mean / xp.maximum(mean, eps))


def make_input_gaussian(
    shape: tuple[int, int],
    dx_doe_m: float,
    gaussian_1e2_diameter_m: float,
    clear_aperture_m: float,
    xp: Any,
    dtype: Any,
) -> Any:
    """构建 DOE 平面高斯输入振幅: A(r)=exp(-r²/w²), 加通光孔径, 归一化总功率=1."""
    Ny, Nx = int(shape[0]), int(shape[1])
    x = (xp.arange(Nx, dtype=dtype) - (Nx // 2)) * dtype(dx_doe_m)
    y = (xp.arange(Ny, dtype=dtype) - (Ny // 2)) * dtype(dx_doe_m)
    X, Y = xp.meshgrid(x, y)
    r2 = X * X + Y * Y
    w = dtype(gaussian_1e2_diameter_m / 2.0)
    aperture_radius = dtype(clear_aperture_m / 2.0)
    amp = xp.exp(-r2 / (w * w)).astype(dtype, copy=False)
    amp = xp.where(r2 <= aperture_radius * aperture_radius, amp, dtype(0.0))
    return normalize_power(amp, xp, target_power=1.0).astype(dtype, copy=False)
