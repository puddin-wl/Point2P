"""NumPy/CuPy backend selection helpers.
---
GPU/CPU 后端选择模块。通过 ArrayBackend 封装 xp (numpy 或 cupy),
其余代码通过 backend.xp 访问数组模块, 不需要直接 import cupy。
这样保持 CPU/GPU 数组所有权显式可控。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def _decode_gpu_name(raw_name: Any) -> str:
    """从 CuPy 设备属性中解码 GPU 名称字符串."""
    if isinstance(raw_name, bytes):
        return raw_name.decode("utf-8", errors="replace").rstrip("\x00")
    return str(raw_name)


@dataclass
class ArrayBackend:
    """数组后端容器: 封装 xp (numpy/cupy) + 类型/设备信息 + 数据转移方法."""

    xp: Any              # 数组模块: numpy 或 cupy
    name: str            # "numpy" 或 "cupy"
    use_gpu: bool        # 是否使用 GPU
    gpu_name: str | None = None
    device_id: int = 0

    @property
    def float_dtype(self) -> Any:
        """默认浮点类型 (GPU=float32 / CPU=float32)."""
        return self.xp.float32 if self.use_gpu else np.float32

    @property
    def complex_dtype(self) -> Any:
        """默认复数类型 (GPU=complex64 / CPU=complex64)."""
        return self.xp.complex64 if self.use_gpu else np.complex64

    def to_numpy(self, array: Any) -> np.ndarray:
        """将 GPU 数组移回 CPU NumPy (CPU 数组原样返回)."""
        if self.use_gpu and hasattr(array, "get"):
            return array.get()
        return np.asarray(array)

    def to_backend(self, array: Any, dtype: Any | None = None) -> Any:
        """将数组转移到当前后端 (numpy→cupy 或 cupy→numpy), 可选转换 dtype."""
        if dtype is None:
            return self.xp.asarray(array)
        return self.xp.asarray(array, dtype=dtype)

    def describe(self) -> str:
        """后端描述字符串 (GPU: "GPU/CuPy device 0: RTX 5070 Ti" / CPU: "CPU/NumPy")."""
        if self.use_gpu:
            return f"GPU/CuPy device {self.device_id}: {self.gpu_name}"
        return "CPU/NumPy"


def get_backend(use_cupy: bool = True, device_id: int = 0, verbose: bool = True) -> ArrayBackend:
    """获取计算后端: 优先 CuPy/GPU, 失败则回退 NumPy/CPU."""
    if use_cupy:
        try:
            import cupy as cp  # type: ignore

            count = int(cp.cuda.runtime.getDeviceCount())
            if count <= 0:
                raise RuntimeError("CuPy imported, but no CUDA devices were found.")
            if device_id < 0 or device_id >= count:
                raise ValueError(f"Requested CUDA device {device_id}, but only {count} device(s) exist.")
            cp.cuda.Device(device_id).use()
            props = cp.cuda.runtime.getDeviceProperties(device_id)
            gpu_name = _decode_gpu_name(props.get("name", f"device-{device_id}"))
            backend = ArrayBackend(cp, "cupy", True, gpu_name=gpu_name, device_id=device_id)
            if verbose:
                print(f"Using GPU backend: CuPy on {gpu_name}")
            return backend
        except Exception as exc:
            if verbose:
                print(f"CuPy backend unavailable ({exc}); falling back to NumPy CPU.")

    backend = ArrayBackend(np, "numpy", False)
    if verbose:
        print("Using CPU backend: NumPy")
    return backend


def to_numpy(array: Any) -> np.ndarray:
    """将 NumPy/CuPy 数组尽量转为 NumPy (无 GPU 时原样返回)."""
    if hasattr(array, "get"):
        return array.get()
    return np.asarray(array)


def to_backend(array: Any, xp: Any, dtype: Any | None = None) -> Any:
    """将数组转移到指定后端模块 (xp = numpy 或 cupy)."""
    if dtype is None:
        return xp.asarray(array)
    return xp.asarray(array, dtype=dtype)
