"""Export a Y-shift sweep of the complete V2 + Z40/Z20 compensation phase.

Unlike the single baseline exporter, this experimental alignment sweep first
combines V2 with the centered Zernike compensation and then shifts the complete
combined phase. The Zernike correction center therefore moves together with
the V2 pattern. The blaze grating is always added after the installation shift.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from convert_to_slm import compute_dx_doe, convert_to_slm, load_phase
from export_zernike_compensated_slm import (
    DEFAULT_EXTENSION_WIDTH_RHO,
    DEFAULT_PUPIL_DIAMETER_MM,
    DEFAULT_Z20_RMS_WAVES,
    DEFAULT_Z40_RMS_WAVES,
    make_zernike_compensation,
    quantize_phase_to_uint8,
    sha256_file,
    verify_zero_compensation_reference,
)
from shift_sweep_slm import add_blaze, format_shift, shift_phase


DEFAULT_Y_SHIFTS = "-10,-5,0,5,10,15,20"


def parse_shift_values(text: str) -> list[float]:
    """Parse a comma-separated shift list and reject duplicates."""
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("At least one Y shift is required")
    if len(set(values)) != len(values):
        raise ValueError(f"Duplicate Y shifts are not allowed: {text}")
    return values


def output_filename(shift_x_px: float, shift_y_px: float) -> str:
    """Return a direct-load BMP filename with explicit absolute shifts."""
    return (
        "phase_V2_posZernikeComp_"
        f"shiftX{format_shift(shift_x_px)}_"
        f"Y{format_shift(shift_y_px)}.bmp"
    )


def render_contact_sheet(
    records: list[dict[str, Any]],
    output_dir: Path,
    output_path: Path,
) -> None:
    """Render all direct-load phase maps in loading order."""
    columns = 4
    rows = int(np.ceil(len(records) / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(4.2 * columns, 4.0 * rows),
        constrained_layout=True,
        squeeze=False,
    )
    for ax in axes.flat:
        ax.axis("off")
    for ax, record in zip(axes.flat, records):
        gray = np.asarray(Image.open(output_dir / record["filename"]).convert("L"))
        ax.imshow(gray, cmap="gray", vmin=0, vmax=255, origin="upper")
        ax.set_title(
            f"X={record['shift_x_px']:+g}, Y={record['shift_y_px']:+g} px\n"
            f"SHA256 {record['sha256'][:12]}...",
            fontsize=10,
        )
        ax.set_xlabel("SLM x / pixel")
        ax.set_ylabel("SLM y / pixel")
        ax.axis("on")
    fig.suptitle(
        "Whole-phase Y alignment sweep: V2 + positive Z40/Z20, then shift, then blaze",
        fontsize=14,
    )
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export a whole-phase Y shift sweep for V2 + Z40/Z20."
    )
    parser.add_argument("phase", help="Source refined V2 phase (.npy or .mat)")
    parser.add_argument("--out", required=True, help="Direct-load output directory")
    parser.add_argument("--phase-var", default=None)
    parser.add_argument("--y-shifts", default=DEFAULT_Y_SHIFTS)
    parser.add_argument("--fixed-x", type=float, default=5.0)
    parser.add_argument("--z40", type=float, default=DEFAULT_Z40_RMS_WAVES)
    parser.add_argument("--z20", type=float, default=DEFAULT_Z20_RMS_WAVES)
    parser.add_argument(
        "--pupil-diameter-mm", type=float, default=DEFAULT_PUPIL_DIAMETER_MM
    )
    parser.add_argument(
        "--extension-width-rho",
        type=float,
        default=DEFAULT_EXTENSION_WIDTH_RHO,
    )
    parser.add_argument("--blaze-x", type=float, default=200.0)
    parser.add_argument("--blaze-y", type=float, default=200.0)
    parser.add_argument("--f", type=float, default=0.429)
    parser.add_argument("--wavelength", type=float, default=532e-9)
    parser.add_argument("--focal-dx", type=float, default=2.5)
    parser.add_argument("--slm-res", type=int, default=1024)
    parser.add_argument("--slm-pitch", type=float, default=17.0)
    parser.add_argument(
        "--verify-zero-reference",
        type=Path,
        default=None,
        help="Historical zero-Zernike X+5/Y+5 BMP for pipeline regression",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    phase_path = Path(args.phase).resolve()
    output_dir = Path(args.out).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    y_shifts = parse_shift_values(args.y_shifts)

    base_phase, loader = load_phase(str(phase_path), args.phase_var)
    base_phase = np.nan_to_num(base_phase.astype(np.float64), nan=0.0)
    if base_phase.ndim != 2 or base_phase.shape[0] != base_phase.shape[1]:
        raise ValueError(f"Expected a square 2D phase, got {base_phase.shape}")
    n = base_phase.shape[0]
    dx_doe_m = compute_dx_doe(args.wavelength, args.f, n, args.focal_dx)
    pupil_diameter_m = args.pupil_diameter_mm * 1e-3
    compensation_rad, _, _ = make_zernike_compensation(
        base_phase.shape,
        dx_doe_m,
        pupil_diameter_m,
        args.z40,
        args.z20,
        extension_width_rho=args.extension_width_rho,
    )
    combined_centered = np.mod(base_phase + compensation_rad, 2.0 * np.pi)

    zero_reference: dict[str, Any] | None = None
    if args.verify_zero_reference is not None:
        zero_reference = verify_zero_compensation_reference(
            args.verify_zero_reference.resolve(),
            base_phase,
            dx_doe_m,
            5.0,
            5.0,
            args.blaze_x,
            args.blaze_y,
            args.wavelength,
            args.f,
            args.focal_dx,
            args.slm_res,
            args.slm_pitch,
        )
        if not zero_reference["pixel_exact_match"]:
            raise RuntimeError("Historical zero-Zernike regression failed")

    records: list[dict[str, Any]] = []
    slm_half_width_m = 0.5 * args.slm_res * args.slm_pitch * 1e-6
    pupil_radius_m = 0.5 * pupil_diameter_m
    for shift_y_px in y_shifts:
        pupil_margin_x_m = slm_half_width_m - (
            pupil_radius_m + abs(args.fixed_x) * dx_doe_m
        )
        pupil_margin_y_m = slm_half_width_m - (
            pupil_radius_m + abs(shift_y_px) * dx_doe_m
        )
        if min(pupil_margin_x_m, pupil_margin_y_m) < 0.0:
            raise ValueError(
                f"X={args.fixed_x:+g}, Y={shift_y_px:+g} moves the 15 mm "
                "normalization pupil outside the physical SLM"
            )
        shifted = shift_phase(combined_centered, shift_y_px, args.fixed_x)
        blazed = add_blaze(
            shifted,
            dx_doe_m,
            args.wavelength,
            args.f,
            args.blaze_x,
            args.blaze_y,
        )
        phase_slm = convert_to_slm(
            blazed,
            dx_doe_um=dx_doe_m * 1e6,
            N=n,
            lambda_m=args.wavelength,
            f_m=args.f,
            focal_dx_um=args.focal_dx,
            slm_res=args.slm_res,
            slm_pitch_um=args.slm_pitch,
        )
        gray = quantize_phase_to_uint8(phase_slm)
        filename = output_filename(args.fixed_x, shift_y_px)
        bmp_path = output_dir / filename
        Image.fromarray(gray, mode="L").save(bmp_path, format="BMP")
        records.append(
            {
                "loading_order": len(records) + 1,
                "filename": filename,
                "shift_x_px": args.fixed_x,
                "shift_y_px": shift_y_px,
                "shift_x_um_on_doe_grid": args.fixed_x * dx_doe_m * 1e6,
                "shift_y_um_on_doe_grid": shift_y_px * dx_doe_m * 1e6,
                "normalization_pupil_margin_x_mm": pupil_margin_x_m * 1e3,
                "normalization_pupil_margin_y_mm": pupil_margin_y_m * 1e3,
                "minimum_normalization_pupil_margin_mm": (
                    min(pupil_margin_x_m, pupil_margin_y_m) * 1e3
                ),
                "gray_min": int(gray.min()),
                "gray_max": int(gray.max()),
                "unique_gray_levels": int(np.unique(gray).size),
                "sha256": sha256_file(bmp_path),
            }
        )
        print(f"Saved Y={shift_y_px:+g}: {bmp_path}")

    hashes = [record["sha256"] for record in records]
    if len(set(hashes)) != len(hashes):
        raise RuntimeError("The Y sweep produced duplicate BMP files")

    contact_sheet = output_dir / "Y_SHIFT_SWEEP_CONTACT_SHEET.png"
    render_contact_sheet(records, output_dir, contact_sheet)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "description": "Y shift sweep of the complete V2 + positive Z40/Z20 phase",
        "source_phase": {
            "path": str(phase_path),
            "loader": loader,
            "sha256": sha256_file(phase_path),
        },
        "shift_semantics": (
            "V2 and the centered Zernike compensation are combined first; "
            "the complete combined phase is then shifted. Positive Y moves "
            "the complete phase downward on the SLM. Blaze is added last."
        ),
        "zernike_compensation": {
            "Z40_primary_spherical_rms_waves": args.z40,
            "Z20_defocus_rms_waves": args.z20,
            "pupil_diameter_mm": args.pupil_diameter_mm,
            "extension_width_rho": args.extension_width_rho,
            "astigmatism_or_coma_added": False,
        },
        "grid": {
            "N": n,
            "dx_doe_um": dx_doe_m * 1e6,
            "wavelength_m": args.wavelength,
            "focal_length_m": args.f,
            "focal_dx_um": args.focal_dx,
        },
        "fixed_x_shift_px": args.fixed_x,
        "y_shifts_px": y_shifts,
        "blaze_grating_focal_shift_um": {
            "x": args.blaze_x,
            "y": args.blaze_y,
        },
        "slm": {
            "resolution": [args.slm_res, args.slm_res],
            "pixel_pitch_um": args.slm_pitch,
            "bit_depth": 8,
        },
        "all_output_hashes_unique": True,
        "zero_compensation_reference_verification": zero_reference,
        "contact_sheet": contact_sheet.name,
        "files": records,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    readme_lines = [
        "V2 + positive Z40/Z20 whole-phase Y shift sweep",
        "",
        f"Z40 = {args.z40:+.5f} RMS waves",
        f"Z20 = {args.z20:+.5f} RMS waves",
        f"Fixed X shift = {args.fixed_x:+g} computational pixels",
        f"Y shifts = {', '.join(format_shift(value) for value in y_shifts)} pixels",
        f"Blaze = X{args.blaze_x:+g}/Y{args.blaze_y:+g} um in the focal plane",
        "",
        "Shift convention: V2 and Zernike compensation move together.",
        "Load in the following order and record the measured result for each BMP:",
        "",
    ]
    for record in records:
        readme_lines.append(
            f"{record['loading_order']}. Y={record['shift_y_px']:+g}: "
            f"{record['filename']}  SHA256={record['sha256']}"
        )
    (output_dir / "README_EXPERIMENT.txt").write_text(
        "\n".join(readme_lines) + "\n",
        encoding="utf-8",
    )
    print(f"Generated {len(records)} unique direct-load BMPs in {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
