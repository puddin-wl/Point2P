"""Export a V2 SLM phase with Z40 spherical and Z20 defocus compensation.

Processing order:

    refined V2 phase
    -> installation shift of the V2 pattern
    -> centered Zernike compensation on the optical pupil
    -> blaze grating
    -> physical SLM crop/resampling
    -> 1024 x 1024, 8-bit BMP

The default positive coefficients are the nominal opposite of the negative
system aberration inferred from the 2026-07-16 hollow-flattop experiment.
Coefficients use Noll-normalized real Zernike definitions and are expressed in
RMS waves on the configured circular pupil.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from convert_to_slm import compute_dx_doe, convert_to_slm, load_phase, save_outputs
from shift_sweep_slm import add_blaze, shift_phase


DEFAULT_Z40_RMS_WAVES = +0.10625
DEFAULT_Z20_RMS_WAVES = +0.25000
DEFAULT_PUPIL_DIAMETER_MM = 15.0
DEFAULT_EXTENSION_WIDTH_RHO = 0.20


def sha256_file(path: Path) -> str:
    """Return the SHA256 digest of one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def format_value(value: float, digits: int = 5) -> str:
    """Format a signed float for a stable, filesystem-safe label."""
    text = f"{value:+.{digits}f}".replace(".", "p")
    return text.rstrip("0").rstrip("p")


def default_label(z40_rms_waves: float, z20_rms_waves: float) -> str:
    """Build the default filename prefix."""
    return (
        "phase_V2_posSphericalComp_"
        f"Z40{format_value(z40_rms_waves)}_"
        f"Z20{format_value(z20_rms_waves)}"
    )


def make_zernike_compensation(
    shape: tuple[int, int],
    dx_doe_m: float,
    pupil_diameter_m: float,
    z40_rms_waves: float,
    z20_rms_waves: float,
    center_x_m: float = 0.0,
    center_y_m: float = 0.0,
    extension_width_rho: float = DEFAULT_EXTENSION_WIDTH_RHO,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create a centered Noll-normalized Z40+Z20 compensation map.

    Returns ``(phase_rad, phase_waves, pupil_mask)``. The compensation is
    unmodified within the normalization pupil:

      Z20 = sqrt(3) * (2 rho^2 - 1)
      Z40 = sqrt(5) * (6 rho^4 - 6 rho^2 + 1)

    Outside ``rho=1``, the radial coordinate is smoothly saturated over
    ``extension_width_rho`` and then held constant. The saturation polynomial
    is C2-continuous at both ends, so the physical SLM phase has no artificial
    circular phase step at the 15 mm normalization boundary.
    """
    ny, nx = shape
    if ny <= 0 or nx <= 0:
        raise ValueError(f"Invalid phase shape: {shape}")
    if dx_doe_m <= 0.0:
        raise ValueError("dx_doe_m must be positive")
    if pupil_diameter_m <= 0.0:
        raise ValueError("pupil_diameter_m must be positive")
    if extension_width_rho <= 0.0:
        raise ValueError("extension_width_rho must be positive")

    x_m = (np.arange(nx, dtype=np.float64) - nx // 2) * dx_doe_m
    y_m = (np.arange(ny, dtype=np.float64) - ny // 2) * dx_doe_m
    x_grid_m, y_grid_m = np.meshgrid(x_m - center_x_m, y_m - center_y_m)
    pupil_radius_m = 0.5 * pupil_diameter_m
    rho = np.sqrt(x_grid_m * x_grid_m + y_grid_m * y_grid_m) / (
        pupil_radius_m
    )
    pupil_mask = rho <= 1.0

    # Preserve the exact Zernike coordinate inside the normalization pupil.
    # Outside, use g(t)=t-t^3+0.5*t^4. It satisfies:
    #   g(0)=0, g'(0)=1, g''(0)=0
    #   g'(1)=0, g''(1)=0
    # and therefore joins the original radial coordinate to a constant with
    # continuous value, slope, and curvature.
    t = np.clip((rho - 1.0) / extension_width_rho, 0.0, 1.0)
    smooth_extension = t - t**3 + 0.5 * t**4
    rho_evaluated = np.where(
        pupil_mask,
        rho,
        1.0 + extension_width_rho * smooth_extension,
    )
    rho2 = rho_evaluated * rho_evaluated

    z20 = math.sqrt(3.0) * (2.0 * rho2 - 1.0)
    z40 = math.sqrt(5.0) * (6.0 * rho2 * rho2 - 6.0 * rho2 + 1.0)
    phase_waves = z40_rms_waves * z40 + z20_rms_waves * z20
    phase_rad = 2.0 * np.pi * phase_waves
    return (
        phase_rad.astype(np.float64),
        phase_waves.astype(np.float64),
        pupil_mask,
    )


def quantize_phase_to_uint8(phase_rad: np.ndarray) -> np.ndarray:
    """Apply the same 8-bit phase quantization used by ``save_outputs``."""
    return (
        np.mod(phase_rad, 2.0 * np.pi) / (2.0 * np.pi) * 255.0
    ).clip(0, 255).astype(np.uint8)


def phase_statistics(phase: np.ndarray) -> dict[str, Any]:
    """Return compact finite/range statistics for a phase array."""
    return {
        "shape": list(phase.shape),
        "dtype": str(phase.dtype),
        "all_finite": bool(np.isfinite(phase).all()),
        "min_rad": float(np.min(phase)),
        "max_rad": float(np.max(phase)),
    }


def pupil_boundary_step_rad(
    phase_rad: np.ndarray,
    dx_doe_m: float,
    pupil_diameter_m: float,
) -> float:
    """Measure the wrapped phase step across +X at the pupil boundary."""
    center_y = phase_rad.shape[0] // 2
    center_x = phase_rad.shape[1] // 2
    pupil_radius_px = 0.5 * pupil_diameter_m / dx_doe_m
    inside_x = center_x + int(np.floor(pupil_radius_px))
    outside_x = inside_x + 1
    difference = np.angle(
        np.exp(
            1j
            * (
                phase_rad[center_y, outside_x]
                - phase_rad[center_y, inside_x]
            )
        )
    )
    return float(abs(difference))


def hard_zero_boundary_step_rad(
    phase_rad: np.ndarray,
    dx_doe_m: float,
    pupil_diameter_m: float,
) -> float:
    """Return the counterfactual boundary step if the exterior were zero."""
    center_y = phase_rad.shape[0] // 2
    center_x = phase_rad.shape[1] // 2
    pupil_radius_px = 0.5 * pupil_diameter_m / dx_doe_m
    inside_x = center_x + int(np.floor(pupil_radius_px))
    difference = np.angle(np.exp(-1j * phase_rad[center_y, inside_x]))
    return float(abs(difference))


def render_preview(
    output_path: Path,
    compensation_waves: np.ndarray,
    phase_slm: np.ndarray,
    dx_doe_m: float,
    pupil_diameter_m: float,
    extension_width_rho: float,
) -> None:
    """Render the compensation map, radial profile, and final SLM phase."""
    ny, nx = compensation_waves.shape
    x_mm = (np.arange(nx) - nx // 2) * dx_doe_m * 1e3
    y_mm = (np.arange(ny) - ny // 2) * dx_doe_m * 1e3
    preview_half_width_mm = 0.6 * pupil_diameter_m * 1e3
    preview_x = np.abs(x_mm) <= preview_half_width_mm
    preview_y = np.abs(y_mm) <= preview_half_width_mm
    preview_waves = compensation_waves[np.ix_(preview_y, preview_x)]
    center_y = ny // 2

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), constrained_layout=True)
    image = axes[0].imshow(
        preview_waves,
        origin="lower",
        extent=[
            x_mm[preview_x][0],
            x_mm[preview_x][-1],
            y_mm[preview_y][0],
            y_mm[preview_y][-1],
        ],
        cmap="coolwarm",
    )
    pupil_radius_mm = 0.5 * pupil_diameter_m * 1e3
    extension_radius_mm = pupil_radius_mm * (1.0 + extension_width_rho)
    axes[0].add_patch(
        plt.Circle(
            (0.0, 0.0),
            pupil_radius_mm,
            fill=False,
            color="black",
            ls="--",
            lw=1.0,
            label="15 mm normalization pupil",
        )
    )
    axes[0].add_patch(
        plt.Circle(
            (0.0, 0.0),
            extension_radius_mm,
            fill=False,
            color="black",
            ls=":",
            lw=0.9,
            label="smooth extension end",
        )
    )
    axes[0].set_title("Z40 + Z20 compensation / waves\n(no hard pupil edge)")
    axes[0].set_xlabel("x / mm")
    axes[0].set_ylabel("y / mm")
    axes[0].legend(loc="upper center", fontsize=8)
    fig.colorbar(image, ax=axes[0], label="wavefront / waves")

    radius_mm = x_mm
    radial = compensation_waves[center_y]
    in_preview = np.abs(radius_mm) <= preview_half_width_mm
    axes[1].plot(radius_mm[in_preview], radial[in_preview], color="black")
    axes[1].axhline(0.0, color="gray", ls=":")
    for signed_radius in (-pupil_radius_mm, pupil_radius_mm):
        axes[1].axvline(signed_radius, color="black", ls="--", lw=0.9)
    for signed_radius in (-extension_radius_mm, extension_radius_mm):
        axes[1].axvline(signed_radius, color="black", ls=":", lw=0.8)
    axes[1].set_title("Horizontal profile with smooth exterior")
    axes[1].set_xlabel("pupil coordinate / mm")
    axes[1].set_ylabel("wavefront / waves")
    axes[1].grid(True, alpha=0.25)

    slm_gray = quantize_phase_to_uint8(phase_slm)
    slm_image = axes[2].imshow(slm_gray, origin="upper", cmap="gray", vmin=0, vmax=255)
    axes[2].set_title("Final 1024 x 1024 SLM phase")
    axes[2].set_xlabel("SLM x / pixel")
    axes[2].set_ylabel("SLM y / pixel")
    fig.colorbar(slm_image, ax=axes[2], label="8-bit gray")
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def verify_zero_compensation_reference(
    reference_bmp: Path,
    base_phase: np.ndarray,
    dx_doe_m: float,
    install_shift_x_px: float,
    install_shift_y_px: float,
    blaze_x_um: float,
    blaze_y_um: float,
    wavelength_m: float,
    focal_length_m: float,
    focal_dx_um: float,
    slm_res: int,
    slm_pitch_um: float,
) -> dict[str, Any]:
    """Rebuild the historical zero-Zernike phase and compare it pixelwise."""
    shifted = shift_phase(base_phase, install_shift_y_px, install_shift_x_px)
    blazed = add_blaze(
        shifted,
        dx_doe_m,
        wavelength_m,
        focal_length_m,
        blaze_x_um,
        blaze_y_um,
    )
    zero_slm = convert_to_slm(
        blazed,
        dx_doe_um=dx_doe_m * 1e6,
        N=base_phase.shape[0],
        lambda_m=wavelength_m,
        f_m=focal_length_m,
        focal_dx_um=focal_dx_um,
        slm_res=slm_res,
        slm_pitch_um=slm_pitch_um,
    )
    rebuilt_gray = quantize_phase_to_uint8(zero_slm)
    reference_gray = np.asarray(Image.open(reference_bmp).convert("L"))
    if reference_gray.shape != rebuilt_gray.shape:
        raise ValueError(
            f"Reference BMP shape {reference_gray.shape} does not match "
            f"rebuilt shape {rebuilt_gray.shape}"
        )
    difference = rebuilt_gray.astype(np.int16) - reference_gray.astype(np.int16)
    return {
        "reference_bmp": str(reference_bmp.resolve()),
        "reference_sha256": sha256_file(reference_bmp),
        "pixel_exact_match": bool(np.array_equal(rebuilt_gray, reference_gray)),
        "different_pixel_count": int(np.count_nonzero(difference)),
        "max_abs_gray_difference": int(np.max(np.abs(difference))),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export V2 with Noll-normalized Z40/Z20 SLM compensation."
    )
    parser.add_argument("phase", help="Source refined V2 phase (.npy or .mat)")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument("--label", default=None, help="Output filename prefix")
    parser.add_argument("--phase-var", default=None, help="Variable name for MAT input")
    parser.add_argument("--z40", type=float, default=DEFAULT_Z40_RMS_WAVES)
    parser.add_argument("--z20", type=float, default=DEFAULT_Z20_RMS_WAVES)
    parser.add_argument(
        "--pupil-diameter-mm", type=float, default=DEFAULT_PUPIL_DIAMETER_MM
    )
    parser.add_argument(
        "--extension-width-rho",
        type=float,
        default=DEFAULT_EXTENSION_WIDTH_RHO,
        help=(
            "C2 smooth-extension width outside the normalization pupil, in "
            "normalized radius (default: 0.20)"
        ),
    )
    parser.add_argument(
        "--zernike-center-x-mm", type=float, default=0.0,
        help="Zernike pupil center relative to the computational grid",
    )
    parser.add_argument(
        "--zernike-center-y-mm", type=float, default=0.0,
        help="Zernike pupil center relative to the computational grid",
    )
    parser.add_argument("--install-shift-x", type=float, default=5.0)
    parser.add_argument("--install-shift-y", type=float, default=5.0)
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
        help="Historical zero-Zernike BMP to rebuild and compare pixelwise",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    phase_path = Path(args.phase).resolve()
    output_dir = Path(args.out).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    base_phase, loader = load_phase(str(phase_path), args.phase_var)
    base_phase = np.nan_to_num(base_phase.astype(np.float64), nan=0.0)
    if base_phase.ndim != 2 or base_phase.shape[0] != base_phase.shape[1]:
        raise ValueError(f"Expected a square 2D phase, got {base_phase.shape}")
    n = base_phase.shape[0]
    dx_doe_m = compute_dx_doe(args.wavelength, args.f, n, args.focal_dx)
    pupil_diameter_m = args.pupil_diameter_mm * 1e-3

    compensation_rad, compensation_waves, pupil_mask = make_zernike_compensation(
        base_phase.shape,
        dx_doe_m,
        pupil_diameter_m,
        args.z40,
        args.z20,
        args.zernike_center_x_mm * 1e-3,
        args.zernike_center_y_mm * 1e-3,
        args.extension_width_rho,
    )

    shifted_v2 = shift_phase(
        base_phase,
        args.install_shift_y,
        args.install_shift_x,
    )
    compensated_pre_blaze = np.mod(
        shifted_v2 + compensation_rad,
        2.0 * np.pi,
    )
    compensated_blazed = add_blaze(
        compensated_pre_blaze,
        dx_doe_m,
        args.wavelength,
        args.f,
        args.blaze_x,
        args.blaze_y,
    )
    phase_slm = convert_to_slm(
        compensated_blazed,
        dx_doe_um=dx_doe_m * 1e6,
        N=n,
        lambda_m=args.wavelength,
        f_m=args.f,
        focal_dx_um=args.focal_dx,
        slm_res=args.slm_res,
        slm_pitch_um=args.slm_pitch,
    )

    label = args.label or default_label(args.z40, args.z20)
    save_outputs(phase_slm, str(output_dir), label=label)
    np.save(
        output_dir / "zernike_compensation_waves_2048.npy",
        compensation_waves.astype(np.float32),
    )
    np.save(
        output_dir / "phase_compensated_pre_blaze_2048.npy",
        compensated_pre_blaze.astype(np.float32),
    )
    np.save(
        output_dir / "phase_compensated_blazed_2048.npy",
        compensated_blazed.astype(np.float32),
    )

    bmp_path = output_dir / f"{label}.bmp"
    preview_path = output_dir / "zernike_compensation_preview.png"
    render_preview(
        preview_path,
        compensation_waves,
        phase_slm,
        dx_doe_m,
        pupil_diameter_m,
        args.extension_width_rho,
    )

    zero_reference: dict[str, Any] | None = None
    compensated_vs_reference: dict[str, Any] | None = None
    if args.verify_zero_reference is not None:
        zero_reference = verify_zero_compensation_reference(
            args.verify_zero_reference.resolve(),
            base_phase,
            dx_doe_m,
            args.install_shift_x,
            args.install_shift_y,
            args.blaze_x,
            args.blaze_y,
            args.wavelength,
            args.f,
            args.focal_dx,
            args.slm_res,
            args.slm_pitch,
        )
        if not zero_reference["pixel_exact_match"]:
            raise RuntimeError(
                "Zero-compensation pipeline does not reproduce the supplied "
                "historical BMP exactly"
            )

    slm_gray = quantize_phase_to_uint8(phase_slm)
    if args.verify_zero_reference is not None:
        reference_gray = np.asarray(
            Image.open(args.verify_zero_reference).convert("L")
        )
        gray_difference = (
            slm_gray.astype(np.int16) - reference_gray.astype(np.int16)
        )
        compensated_vs_reference = {
            "different_pixel_count": int(np.count_nonzero(gray_difference)),
            "different_pixel_fraction": float(
                np.count_nonzero(gray_difference) / gray_difference.size
            ),
            "mean_abs_gray_difference": float(
                np.mean(np.abs(gray_difference))
            ),
            "max_abs_gray_difference": int(
                np.max(np.abs(gray_difference))
            ),
        }
    in_pupil_waves = compensation_waves[pupil_mask]
    nominal_system_rad = -compensation_rad[pupil_mask]
    cancellation_residual = np.angle(
        np.exp(1j * nominal_system_rad)
        * np.exp(1j * compensation_rad[pupil_mask])
    )
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "description": (
            "V2 phase with positive Z40 spherical and positive Z20 defocus "
            "compensation for the inferred negative experimental aberration"
        ),
        "processing_order": [
            "load_refined_v2",
            "shift_v2_for_installation",
            "add_centered_zernike_compensation",
            "add_blaze_grating",
            "crop_and_complex_resample_to_slm",
            "quantize_to_8bit_bmp",
        ],
        "source": {
            "phase": str(phase_path),
            "loader": loader,
            "sha256": sha256_file(phase_path),
            "statistics": phase_statistics(base_phase),
        },
        "grid": {
            "N": n,
            "wavelength_m": args.wavelength,
            "focal_length_m": args.f,
            "focal_dx_um": args.focal_dx,
            "dx_doe_um": dx_doe_m * 1e6,
        },
        "zernike_compensation": {
            "sign_intent": "positive compensation for inferred negative aberration",
            "Z40_primary_spherical_rms_waves": args.z40,
            "Z20_defocus_rms_waves": args.z20,
            "pupil_diameter_mm": args.pupil_diameter_mm,
            "pupil_center_x_mm": args.zernike_center_x_mm,
            "pupil_center_y_mm": args.zernike_center_y_mm,
            "Z20_definition": "sqrt(3)*(2*rho^2-1)",
            "Z40_definition": "sqrt(5)*(6*rho^4-6*rho^2+1)",
            "outside_pupil": (
                "C2-continuous radial-coordinate saturation followed by a "
                "constant piston; no hard phase step"
            ),
            "extension_width_rho": args.extension_width_rho,
            "extension_width_mm": (
                0.5
                * args.pupil_diameter_mm
                * args.extension_width_rho
            ),
            "extension_end_radius_mm": (
                0.5
                * args.pupil_diameter_mm
                * (1.0 + args.extension_width_rho)
            ),
            "sampled_wrapped_step_across_pupil_boundary_rad": (
                pupil_boundary_step_rad(
                    compensation_rad,
                    dx_doe_m,
                    pupil_diameter_m,
                )
            ),
            "counterfactual_hard_zero_boundary_step_rad": (
                hard_zero_boundary_step_rad(
                    compensation_rad,
                    dx_doe_m,
                    pupil_diameter_m,
                )
            ),
            "wavefront_min_waves_inside_pupil": float(np.min(in_pupil_waves)),
            "wavefront_max_waves_inside_pupil": float(np.max(in_pupil_waves)),
            "wavefront_mean_waves_inside_pupil": float(np.mean(in_pupil_waves)),
            "wavefront_rms_waves_inside_pupil": float(
                np.sqrt(np.mean(np.square(in_pupil_waves)))
            ),
            "wavefront_peak_to_valley_waves_inside_pupil": float(
                np.ptp(in_pupil_waves)
            ),
            "nominal_opposite_aberration_cancellation_max_abs_rad": float(
                np.max(np.abs(cancellation_residual))
            ),
        },
        "installation_shift": {
            "x_computational_px": args.install_shift_x,
            "y_computational_px": args.install_shift_y,
            "x_um_on_doe_grid": args.install_shift_x * dx_doe_m * 1e6,
            "y_um_on_doe_grid": args.install_shift_y * dx_doe_m * 1e6,
        },
        "blaze_grating_focal_shift_um": {
            "x": args.blaze_x,
            "y": args.blaze_y,
        },
        "slm_output": {
            "bmp": bmp_path.name,
            "sha256": sha256_file(bmp_path),
            "resolution": [args.slm_res, args.slm_res],
            "pixel_pitch_um": args.slm_pitch,
            "bit_depth": 8,
            "gray_min": int(slm_gray.min()),
            "gray_max": int(slm_gray.max()),
            "unique_gray_levels": int(np.unique(slm_gray).size),
            "phase_statistics": phase_statistics(phase_slm),
        },
        "preview": preview_path.name,
        "zero_compensation_reference_verification": zero_reference,
        "compensated_output_vs_zero_reference": compensated_vs_reference,
        "experimental_note": (
            "These are initial simulation-derived coefficients. Keep Z40/Z20 "
            "coupled during experimental fine tuning. Astigmatism and coma are "
            "not included."
        ),
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    readme_path = output_dir / "README_EXPERIMENT.txt"
    readme_path.write_text(
        "\n".join(
            [
                "V2 positive spherical-aberration compensation phase",
                "",
                f"Direct-load BMP: {bmp_path.name}",
                f"SHA256: {manifest['slm_output']['sha256']}",
                f"Z40 compensation: {args.z40:+.5f} RMS waves",
                f"Z20 compensation: {args.z20:+.5f} RMS waves",
                f"Pupil diameter: {args.pupil_diameter_mm:.3f} mm",
                (
                    "Outside-pupil extension: C2 smooth over "
                    f"{args.extension_width_rho:.3f} rho "
                    f"({0.5 * args.pupil_diameter_mm * args.extension_width_rho:.3f} mm)"
                ),
                (
                    "Installation shift: "
                    f"X={args.install_shift_x:+g}, Y={args.install_shift_y:+g} "
                    "computational pixels"
                ),
                f"Blaze focal shift: X={args.blaze_x:+g} um, Y={args.blaze_y:+g} um",
                "",
                "The Zernike compensation is added after shifting the V2 pattern and",
                "before adding the blaze grating. It remains centered on the optical pupil.",
                "The 15 mm normalization boundary has no hard phase step; the exterior",
                "is smoothly extended to a constant piston phase.",
                "This first coefficient pair comes from simulation and must be fine-tuned",
                "experimentally. Astigmatism and coma are not included.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"Direct-load BMP: {bmp_path}")
    print(f"SHA256: {manifest['slm_output']['sha256']}")
    print(f"Manifest: {manifest_path}")
    print(f"Preview: {preview_path}")
    if zero_reference is not None:
        print("Historical zero-Zernike pipeline reproduction: exact pixel match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
