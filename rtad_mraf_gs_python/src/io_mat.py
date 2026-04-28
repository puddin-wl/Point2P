"""MATLAB MAT-file input/output helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy.io import loadmat, savemat


def list_mat_variables(path: str | Path) -> list[str]:
    """List top-level variables in a MATLAB MAT file."""
    path = Path(path)
    try:
        data = loadmat(path, variable_names=None)
        return sorted(k for k in data.keys() if not k.startswith("__"))
    except Exception:
        with h5py.File(path, "r") as h5:
            return sorted(k for k in h5.keys() if not k.startswith("#"))


def _read_h5_variable(path: Path, var_name: str, transpose_h5: bool) -> np.ndarray:
    """Read one variable from a MATLAB v7.3 HDF5 MAT file."""
    with h5py.File(path, "r") as h5:
        if var_name not in h5:
            variables = sorted(k for k in h5.keys() if not k.startswith("#"))
            raise KeyError(f"Variable '{var_name}' not found in {path}. Available: {variables}")
        array = h5[var_name][()]
    if transpose_h5 and array.ndim >= 2:
        array = array.T
    return np.asarray(array)


def load_mat_variable(
    path: str | Path,
    var_name: str,
    transpose_h5: bool = False,
    squeeze: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Read one variable from either classic MAT or v7.3 HDF5 MAT files.

    Parameters
    ----------
    path:
        MAT file path.
    var_name:
        Variable name to read.
    transpose_h5:
        Transpose HDF5 arrays after reading. MATLAB v7.3 files can appear
        transposed depending on how they were written and consumed.
    squeeze:
        If ``True``, remove singleton dimensions after reading.
    """
    path = Path(path)
    info: dict[str, Any] = {
        "path": str(path),
        "variable": var_name,
        "transpose_h5": bool(transpose_h5),
    }
    try:
        data = loadmat(path, variable_names=[var_name])
    except Exception as exc:
        array = _read_h5_variable(path, var_name, transpose_h5=transpose_h5)
        info["reader"] = "h5py"
        info["loadmat_error"] = f"{type(exc).__name__}: {exc}"
    else:
        if var_name not in data:
            variables = list_mat_variables(path)
            raise KeyError(f"Variable '{var_name}' not found in {path}. Available: {variables}")
        array = np.asarray(data[var_name])
        info["reader"] = "scipy.io.loadmat"

    if squeeze:
        array = np.squeeze(array)

    info["shape"] = tuple(int(v) for v in array.shape)
    info["dtype"] = str(array.dtype)
    return array, info


def load_optional_mat_variable(
    path: str | Path,
    var_name: str,
    transpose_h5: bool = False,
    squeeze: bool = True,
) -> tuple[np.ndarray | None, dict[str, Any] | None]:
    """Read a variable if it exists; return ``(None, None)`` otherwise."""
    try:
        return load_mat_variable(path, var_name, transpose_h5=transpose_h5, squeeze=squeeze)
    except Exception:
        return None, None


def load_phase_mat(
    path: str | Path,
    phase_var: str = "phase0",
    transpose_h5: bool = False,
    swap_xy: bool = False,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load a phase matrix in radians and wrap it to ``[0, 2*pi)``.

    The function first tries ``scipy.io.loadmat`` for MATLAB v7 and older
    files. If that fails, it reads MATLAB v7.3 HDF5 files with ``h5py``.
    Set ``swap_xy=True`` when a MATLAB phase matrix is observed to have x/y
    axes swapped after import; this applies a 2D transpose to the phase only.
    """
    phase, info = load_mat_variable(path, phase_var, transpose_h5=transpose_h5, squeeze=False)
    if phase.ndim != 2:
        raise ValueError(f"Phase variable '{phase_var}' must be 2D; got shape {phase.shape}.")
    if not np.issubdtype(phase.dtype, np.number):
        raise TypeError(f"Phase variable '{phase_var}' must be numeric; got dtype {phase.dtype}.")

    info["shape_before_phase_xy_swap"] = tuple(int(v) for v in phase.shape)
    info["phase_xy_swapped"] = bool(swap_xy)
    if swap_xy:
        phase = phase.T

    phase = np.mod(np.asarray(phase, dtype=np.float32), np.float32(2.0 * np.pi))
    nonfinite = ~np.isfinite(phase)
    nonfinite_count = int(np.count_nonzero(nonfinite))
    if nonfinite_count:
        phase[nonfinite] = 0.0
    info["nonfinite_replaced_with_zero"] = nonfinite_count
    info["wrapped_to"] = "[0, 2*pi)"
    info["final_shape"] = tuple(int(v) for v in phase.shape)
    info["final_dtype"] = str(phase.dtype)
    info["min_rad"] = float(np.nanmin(phase))
    info["max_rad"] = float(np.nanmax(phase))
    return phase, info


def _jsonable(value: Any) -> Any:
    """Convert common NumPy objects to JSON-serializable values."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def save_phase_mat(
    path: str | Path,
    phase_refined: np.ndarray,
    params: dict[str, Any] | None = None,
    metrics: dict[str, Any] | None = None,
    reconstruction_intensity: np.ndarray | None = None,
) -> None:
    """Save refined phase, parameters, and summary metrics to a MAT file."""
    path = Path(path)
    payload: dict[str, Any] = {
        "phase_refined": np.asarray(phase_refined, dtype=np.float32),
        "phase_refined_rad": np.asarray(phase_refined, dtype=np.float32),
    }
    if reconstruction_intensity is not None:
        payload["reconstruction_intensity"] = np.asarray(reconstruction_intensity, dtype=np.float32)
    if params is not None:
        payload["params_json"] = json.dumps(_jsonable(params), indent=2)
    if metrics is not None:
        numeric = {
            k: float(v)
            for k, v in metrics.items()
            if isinstance(v, (int, float, np.integer, np.floating)) and np.isfinite(v)
        }
        payload["metrics_json"] = json.dumps(_jsonable(metrics), indent=2)
        payload["metrics_names"] = np.asarray(list(numeric.keys()), dtype=object)
        payload["metrics_values"] = np.asarray(list(numeric.values()), dtype=np.float64)
    savemat(path, payload, do_compression=True)
