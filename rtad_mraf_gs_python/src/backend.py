"""NumPy/CuPy backend selection helpers.

The rest of this package receives a small :class:`ArrayBackend` object instead
of importing CuPy directly. This keeps CPU/GPU array ownership explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def _decode_gpu_name(raw_name: Any) -> str:
    """Return a readable GPU name from CuPy runtime device properties."""
    if isinstance(raw_name, bytes):
        return raw_name.decode("utf-8", errors="replace").rstrip("\x00")
    return str(raw_name)


@dataclass
class ArrayBackend:
    """Container for the active array module and transfer helpers."""

    xp: Any
    name: str
    use_gpu: bool
    gpu_name: str | None = None
    device_id: int = 0

    @property
    def float_dtype(self) -> Any:
        """Return the default real dtype used by the refinement."""
        return self.xp.float32 if self.use_gpu else np.float32

    @property
    def complex_dtype(self) -> Any:
        """Return the default complex dtype used by the refinement."""
        return self.xp.complex64 if self.use_gpu else np.complex64

    def to_numpy(self, array: Any) -> np.ndarray:
        """Move an array to CPU NumPy without changing values."""
        if self.use_gpu and hasattr(array, "get"):
            return array.get()
        return np.asarray(array)

    def to_backend(self, array: Any, dtype: Any | None = None) -> Any:
        """Move an array to the active backend, optionally changing dtype."""
        if dtype is None:
            return self.xp.asarray(array)
        return self.xp.asarray(array, dtype=dtype)

    def describe(self) -> str:
        """Return a short backend description for reports."""
        if self.use_gpu:
            return f"GPU/CuPy device {self.device_id}: {self.gpu_name}"
        return "CPU/NumPy"


def get_backend(use_cupy: bool = True, device_id: int = 0, verbose: bool = True) -> ArrayBackend:
    """Select CuPy when available, otherwise fall back to NumPy.

    Parameters
    ----------
    use_cupy:
        If ``True``, try to use CuPy first.
    device_id:
        CUDA device id to use when CuPy is active.
    verbose:
        Print the selected backend.
    """
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
    """Best-effort conversion of NumPy/CuPy arrays to NumPy."""
    if hasattr(array, "get"):
        return array.get()
    return np.asarray(array)


def to_backend(array: Any, xp: Any, dtype: Any | None = None) -> Any:
    """Convert ``array`` to a given backend module."""
    if dtype is None:
        return xp.asarray(array)
    return xp.asarray(array, dtype=dtype)
