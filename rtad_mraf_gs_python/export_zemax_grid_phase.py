"""Export a refined DOE phase map to a Zemax/OpticStudio Grid Phase DAT file.

The exported DAT file follows the same ASCII container format used by Grid Sag
files, as documented by Ansys Optics for OpticStudio Grid Sag/Grid Phase data.
This script:

1. Loads a 2D phase map from .mat or .npy.
2. Infers the computational DOE sampling from ``params_json`` when available.
3. Resamples the centered phase onto a square physical window using complex
   interpolation on ``exp(i*phase)`` to avoid wrap discontinuities.
4. Writes a Grid Phase ``.DAT`` file for use with the Sequential Grid Phase
   surface.

Notes
-----
- The exported phase values are written in radians.
- For a circular DOE aperture, the Grid Phase surface remains a square grid;
  users should still set the OpticStudio clear semi-diameter / aperture
  consistently on the receiving surface.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import loadmat
from scipy.ndimage import map_coordinates

from src.io_mat import load_mat_variable


def compute_dx_doe(lambda_m: float, f_m: float, n: int, focal_dx_um: float) -> float:
    """Infer DOE-plane pixel pitch in meters from focal-plane sampling."""
    return lambda_m * f_m / (n * focal_dx_um * 1e-6)


def _unwrap_mat_string(value: Any) -> str | None:
    """Extract a Python string from common MATLAB string containers."""
    if isinstance(value, str):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"U", "S"}:
            return "".join(np.ravel(value).tolist())
        if value.dtype == object and value.size == 1:
            return _unwrap_mat_string(value.flat[0])
        if value.size == 1:
            return _unwrap_mat_string(value.item())
    return None


def load_phase_and_params(path: Path, phase_var: str) -> tuple[np.ndarray, dict[str, Any] | None]:
    """Load the phase map plus optional params_json metadata from MAT files."""
    phase, _ = load_mat_variable(path, phase_var, squeeze=False)
    phase = np.asarray(phase, dtype=np.float64)
    params: dict[str, Any] | None = None

    if path.suffix.lower() == ".mat":
        try:
            mat = loadmat(path, variable_names=["params_json"])
        except Exception:
            mat = {}
        params_text = _unwrap_mat_string(mat.get("params_json"))
        if params_text:
            try:
                params = json.loads(params_text)
            except json.JSONDecodeError:
                params = None
    return phase, params


def infer_source_dx_mm(
    phase: np.ndarray,
    params: dict[str, Any] | None,
    dx_doe_um: float | None,
    wavelength_m: float,
    focal_length_m: float,
    focal_dx_um: float,
) -> float:
    """Resolve the source sampling pitch in mm."""
    if dx_doe_um is not None:
        return dx_doe_um * 1e-3

    if params:
        grid = params.get("grid", {})
        physical = params.get("physical", {})
        if "dx_doe_m" in grid:
            return float(grid["dx_doe_m"]) * 1e3
        if "wavelength_m" in physical:
            wavelength_m = float(physical["wavelength_m"])
        if "focal_length_m" in physical:
            focal_length_m = float(physical["focal_length_m"])
        if "focal_dx_um" in grid:
            focal_dx_um = float(grid["focal_dx_um"])

    dx_doe_m = compute_dx_doe(
        lambda_m=wavelength_m,
        f_m=focal_length_m,
        n=int(phase.shape[0]),
        focal_dx_um=focal_dx_um,
    )
    return dx_doe_m * 1e3


def resample_phase_window(
    phase: np.ndarray,
    source_dx_mm: float,
    out_width_mm: float,
    out_res: int,
) -> np.ndarray:
    """Resample the centered phase map onto an exact physical output window."""
    if phase.ndim != 2 or phase.shape[0] != phase.shape[1]:
        raise ValueError(f"Expected square 2D phase map, got shape {phase.shape}.")

    n = int(phase.shape[0])
    center_index = n // 2
    coords_mm = np.linspace(-0.5 * out_width_mm, 0.5 * out_width_mm, out_res, dtype=np.float64)
    yy_mm, xx_mm = np.meshgrid(coords_mm, coords_mm, indexing="ij")
    src_y = yy_mm / source_dx_mm + center_index
    src_x = xx_mm / source_dx_mm + center_index
    sample_coords = np.stack([src_y, src_x], axis=0)

    cfield = np.exp(1j * np.mod(phase, 2.0 * np.pi))
    real = map_coordinates(cfield.real, sample_coords, order=3, mode="nearest")
    imag = map_coordinates(cfield.imag, sample_coords, order=3, mode="nearest")
    phase_out = np.mod(np.arctan2(imag, real), 2.0 * np.pi)
    return phase_out


def apply_circular_mask(phase: np.ndarray, out_width_mm: float, radius_mm: float, fill_value: float = 0.0) -> np.ndarray:
    """Zero the phase outside the requested circular clear aperture."""
    coords_mm = np.linspace(-0.5 * out_width_mm, 0.5 * out_width_mm, phase.shape[0], dtype=np.float64)
    yy_mm, xx_mm = np.meshgrid(coords_mm, coords_mm, indexing="ij")
    mask = (xx_mm**2 + yy_mm**2) <= radius_mm**2
    out = np.array(phase, copy=True)
    out[~mask] = fill_value
    return out


def write_grid_phase_dat(
    path: Path,
    phase_rad: np.ndarray,
    width_mm: float,
    unitflag: int = 0,
) -> float:
    """Write a Zemax Grid Phase DAT file and return point spacing in mm."""
    ny, nx = phase_rad.shape
    if nx != ny:
        raise ValueError(f"Expected square phase map, got shape {phase_rad.shape}.")
    delx_mm = width_mm / (nx - 1)
    dely_mm = width_mm / (ny - 1)

    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(f"{nx:d} {ny:d} {delx_mm:.12g} {dely_mm:.12g} {unitflag:d} 0 0\n")
        flat = np.ravel(phase_rad)
        for value in flat:
            f.write(f"{float(value):.16g} 0 0 0 0\n")
    return delx_mm


def write_metadata(
    path: Path,
    *,
    phase_var: str,
    input_path: Path,
    source_shape: tuple[int, int],
    source_dx_mm: float,
    source_extent_mm: float,
    output_res: int,
    output_width_mm: float,
    output_spacing_mm: float,
    circular_aperture_mm: float | None,
) -> None:
    """Save a JSON sidecar for traceability."""
    payload = {
        "input_path": str(input_path),
        "phase_variable": phase_var,
        "source_shape": list(source_shape),
        "source_dx_mm": source_dx_mm,
        "source_extent_mm": source_extent_mm,
        "output_resolution": output_res,
        "output_width_mm": output_width_mm,
        "output_spacing_mm": output_spacing_mm,
        "output_extent_center_to_center_mm": output_spacing_mm * (output_res - 1),
        "circular_aperture_mm": circular_aperture_mm,
        "notes": [
            "Phase values in DAT are written in radians.",
            "Grid Phase uses a square grid; if the physical DOE is circular, set the receiving surface aperture separately in OpticStudio.",
        ],
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export phase_refined to Zemax Grid Phase DAT.")
    parser.add_argument("phase", help="Path to input phase .mat or .npy")
    parser.add_argument("--phase-var", default="phase_refined", help="Variable name when reading a .mat file")
    parser.add_argument("--out", required=True, help="Output DAT file path")
    parser.add_argument("--width-mm", type=float, default=15.0, help="Physical width of the exported square grid")
    parser.add_argument("--res", type=int, default=1024, help="Square output resolution")
    parser.add_argument("--aperture-mm", type=float, default=15.0, help="Optional circular clear-aperture diameter to enforce inside the square grid")
    parser.add_argument("--dx-doe-um", type=float, default=None, help="Override source DOE pitch in microns")
    parser.add_argument("--wavelength-m", type=float, default=532e-9, help="Fallback wavelength when params_json is absent")
    parser.add_argument("--focal-length-m", type=float, default=429e-3, help="Fallback focal length when params_json is absent")
    parser.add_argument("--focal-dx-um", type=float, default=2.5, help="Fallback focal-plane sampling when params_json is absent")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    phase_path = Path(args.phase).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    phase, params = load_phase_and_params(phase_path, args.phase_var)
    source_dx_mm = infer_source_dx_mm(
        phase=phase,
        params=params,
        dx_doe_um=args.dx_doe_um,
        wavelength_m=args.wavelength_m,
        focal_length_m=args.focal_length_m,
        focal_dx_um=args.focal_dx_um,
    )
    source_extent_mm = phase.shape[0] * source_dx_mm

    phase_resampled = resample_phase_window(
        phase=phase,
        source_dx_mm=source_dx_mm,
        out_width_mm=args.width_mm,
        out_res=args.res,
    )

    aperture_mm = float(args.aperture_mm) if args.aperture_mm is not None and args.aperture_mm > 0 else None
    if aperture_mm is not None:
        phase_resampled = apply_circular_mask(
            phase_resampled,
            out_width_mm=args.width_mm,
            radius_mm=0.5 * aperture_mm,
            fill_value=0.0,
        )

    output_spacing_mm = write_grid_phase_dat(
        out_path,
        phase_resampled,
        width_mm=args.width_mm,
        unitflag=0,
    )
    write_metadata(
        out_path.with_suffix(".json"),
        phase_var=args.phase_var,
        input_path=phase_path,
        source_shape=tuple(int(v) for v in phase.shape),
        source_dx_mm=source_dx_mm,
        source_extent_mm=source_extent_mm,
        output_res=args.res,
        output_width_mm=args.width_mm,
        output_spacing_mm=output_spacing_mm,
        circular_aperture_mm=aperture_mm,
    )

    print(f"Loaded: {phase_path}")
    print(f"Source phase shape: {phase.shape}")
    print(f"Source dx_doe: {source_dx_mm * 1e3:.6f} um")
    print(f"Source square extent: {source_extent_mm:.6f} mm")
    print(f"Output DAT: {out_path}")
    print(f"Output grid: {args.res} x {args.res}")
    print(f"Output width: {args.width_mm:.6f} mm")
    print(f"Output spacing (DAT delx=dely): {output_spacing_mm:.12f} mm")
    if aperture_mm is not None:
        print(f"Applied circular aperture mask: {aperture_mm:.6f} mm diameter")
    print(f"Metadata JSON: {out_path.with_suffix('.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
