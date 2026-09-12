"""Prepare the final V2 WGS phase for a 15 mm square, 1024 x 1024 DOE.

The source phase is sampled on the FFT/DOE grid.  Resampling is performed on
exp(1j * phase), rather than directly on wrapped phase, so 0/2*pi boundaries
do not create interpolation artefacts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.io import savemat
from scipy.ndimage import map_coordinates


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SOURCE = (
    PROJECT_ROOT
    / "rtad_mraf_gs_python_test_20260605"
    / "artifacts"
    / "run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm"
    / "phase_refined.npy"
)
DEFAULT_CONFIG = DEFAULT_SOURCE.with_name("config_used.json")
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "v2_wgs_15mm_doe_20260730"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resample_complex_phase(
    phase: np.ndarray,
    source_pitch_m: float,
    target_size_m: float,
    target_n: int,
) -> np.ndarray:
    """Center-sample a physical square and resample its complex phasor."""
    if phase.ndim != 2 or phase.shape[0] != phase.shape[1]:
        raise ValueError(f"Expected a square 2-D phase array, got {phase.shape}")
    source_n = phase.shape[0]
    target_pitch_m = target_size_m / target_n
    target_axis_m = (
        np.arange(target_n, dtype=np.float64) - (target_n - 1) / 2.0
    ) * target_pitch_m
    source_indices = target_axis_m / source_pitch_m + source_n // 2
    yy, xx = np.meshgrid(source_indices, source_indices, indexing="ij")
    phasor = np.exp(1j * np.asarray(phase, dtype=np.float64))
    real = map_coordinates(phasor.real, [yy, xx], order=3, mode="nearest")
    imag = map_coordinates(phasor.imag, [yy, xx], order=3, mode="nearest")
    return np.mod(np.angle(real + 1j * imag), 2.0 * np.pi).astype(np.float32)


def quantize_like_matlab(
    phase: np.ndarray,
    wavelength_nm: float,
    substrate_index: float,
    ambient_index: float,
    etch_step_nm: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Match the supplied MATLAB program's re-zero, wrap and floor rules."""
    phase_rezero = np.asarray(phase, dtype=np.float64) - float(np.min(phase))
    phase_wrapped = np.mod(phase_rezero, 2.0 * np.pi)
    height_nm = (
        wavelength_nm
        / (substrate_index - ambient_index)
        * phase_wrapped
        / (2.0 * np.pi)
    )
    levels = np.floor(height_nm / etch_step_nm).astype(np.uint8)
    levels = np.clip(levels, 0, 15).astype(np.uint8)
    return phase_wrapped.astype(np.float32), levels


def count_run_rectangles(mask: np.ndarray) -> int:
    """Count horizontal runs, the polygons emitted by the MATLAB exporter."""
    padded = np.pad(mask.astype(np.int8), ((0, 0), (1, 1)))
    transitions = np.diff(padded, axis=1)
    return int(np.count_nonzero(transitions == 1))


def write_gds(
    levels: np.ndarray,
    size_mm: float,
    etch_step_nm: float,
    output: Path,
) -> dict[str, int]:
    """Write four binary etch masks as horizontal-run rectangles."""
    try:
        import gdstk
    except ImportError as exc:
        raise RuntimeError("Writing GDS requires the Python package 'gdstk'.") from exc

    nrows, ncols = levels.shape
    pitch_um = size_mm * 1000.0 / ncols
    size_um = size_mm * 1000.0
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    cell = library.new_cell("DOE_V2_WGS_15MM")
    rectangle_counts: dict[str, int] = {}

    for layer, bit in enumerate((8, 4, 2, 1), start=1):
        mask = (levels & bit) != 0
        count = 0
        for row in range(nrows):
            transitions = np.diff(np.r_[False, mask[row], False].astype(np.int8))
            starts = np.flatnonzero(transitions == 1)
            stops = np.flatnonzero(transitions == -1)
            y0 = row * pitch_um
            y1 = (row + 1) * pitch_um
            polygons = [
                gdstk.rectangle(
                    (start * pitch_um, y0),
                    (stop * pitch_um, y1),
                    layer=layer,
                    datatype=0,
                )
                for start, stop in zip(starts, stops)
            ]
            if polygons:
                cell.add(*polygons)
            count += len(polygons)
        rectangle_counts[str(bit)] = count
        print(
            f"GDS layer {layer}: bit {bit}, depth increment "
            f"{bit * etch_step_nm:.1f} nm, {count} rectangles"
        )

    # Layer 100 is a non-etch reference outline for the exact 15 mm square.
    outline_width_um = 1.0
    cell.add(
        gdstk.rectangle((0, 0), (size_um, outline_width_um), layer=100),
        gdstk.rectangle(
            (0, size_um - outline_width_um), (size_um, size_um), layer=100
        ),
        gdstk.rectangle((0, outline_width_um), (outline_width_um, size_um), layer=100),
        gdstk.rectangle(
            (size_um - outline_width_um, outline_width_um),
            (size_um, size_um),
            layer=100,
        ),
    )
    library.write_gds(output)
    return rectangle_counts


def save_preview(
    phase: np.ndarray,
    levels: np.ndarray,
    layer_masks: list[np.ndarray],
    size_mm: float,
    output: Path,
) -> None:
    extent = [-size_mm / 2, size_mm / 2, -size_mm / 2, size_mm / 2]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5), constrained_layout=True)
    im = axes[0, 0].imshow(phase, cmap="twilight", origin="lower", extent=extent)
    axes[0, 0].set_title("V2 WGS phase (wrapped)")
    fig.colorbar(im, ax=axes[0, 0], label="rad")

    im = axes[0, 1].imshow(
        levels, cmap="viridis", origin="lower", extent=extent, vmin=0, vmax=15
    )
    axes[0, 1].set_title("Quantized DOE level")
    fig.colorbar(im, ax=axes[0, 1], label="level")

    axes[0, 2].axis("off")
    axes[0, 2].text(
        0.05,
        0.95,
        "Square DOE\n15.000 mm x 15.000 mm\n"
        f"{phase.shape[1]} x {phase.shape[0]} pixels\n"
        f"{size_mm * 1000 / phase.shape[1]:.7f} um/pixel",
        va="top",
        fontsize=13,
    )

    for axis, mask, bit in zip(axes[1], layer_masks[:3], (8, 4, 2)):
        axis.imshow(mask, cmap="gray", origin="lower", extent=extent, vmin=0, vmax=1)
        axis.set_title(f"Mask bit {bit}")
    for axis in axes.flat:
        if axis.has_data():
            axis.set_xlabel("x (mm)")
            axis.set_ylabel("y (mm)")
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--size-mm", type=float, default=15.0)
    parser.add_argument("--pixels", type=int, default=1024)
    parser.add_argument("--wavelength-nm", type=float, default=532.0)
    parser.add_argument("--substrate-index", type=float, default=1.458)
    parser.add_argument("--ambient-index", type=float, default=1.00029)
    parser.add_argument("--etch-step-nm", type=float, default=77.5)
    parser.add_argument(
        "--skip-gds", action="store_true", help="Prepare MAT/PNG files only"
    )
    args = parser.parse_args()

    with args.config.open("r", encoding="utf-8") as stream:
        config = json.load(stream)
    source_pitch_m = float(config["grid"]["dx_doe_m"])
    phase_source = np.load(args.source)
    phase_resampled = resample_complex_phase(
        phase_source,
        source_pitch_m=source_pitch_m,
        target_size_m=args.size_mm * 1e-3,
        target_n=args.pixels,
    )
    doe_phase, doe_level = quantize_like_matlab(
        phase_resampled,
        wavelength_nm=args.wavelength_nm,
        substrate_index=args.substrate_index,
        ambient_index=args.ambient_index,
        etch_step_nm=args.etch_step_nm,
    )
    bits = (8, 4, 2, 1)
    layer_masks = [((doe_level & bit) != 0) for bit in bits]

    args.output.mkdir(parents=True, exist_ok=True)
    npy_path = args.output / "doe_phase_v2_wgs_15mm_1024.npy"
    mat_path = args.output / "doe_phase_v2_wgs_15mm_1024.mat"
    level_path = args.output / "doe_levels_v2_wgs_15mm_1024.npy"
    phase_png = args.output / "doe_phase_v2_wgs_15mm_16bit.png"
    level_png = args.output / "doe_levels_v2_wgs_15mm.png"
    preview_path = args.output / "doe_v2_wgs_15mm_preview.png"
    gds_path = args.output / "doe_v2_wgs_15mm_square_4mask.gds"

    np.save(npy_path, doe_phase)
    np.save(level_path, doe_level)
    savemat(
        mat_path,
        {
            "doe_phase": doe_phase,
            "doe_phase_rad": doe_phase,
            "doe_level": doe_level,
            "doe_size_mm": np.array([[args.size_mm]], dtype=np.float64),
            "pixel_pitch_um": np.array(
                [[args.size_mm * 1000.0 / args.pixels]], dtype=np.float64
            ),
            "source_dx_um": np.array([[source_pitch_m * 1e6]], dtype=np.float64),
            "wavelength_nm": np.array([[args.wavelength_nm]], dtype=np.float64),
            "substrate_index": np.array([[args.substrate_index]], dtype=np.float64),
            "ambient_index": np.array([[args.ambient_index]], dtype=np.float64),
            "etch_step_nm": np.array([[args.etch_step_nm]], dtype=np.float64),
        },
        do_compression=True,
    )
    Image.fromarray(np.round(doe_phase / (2 * np.pi) * 65535).astype(np.uint16)).save(
        phase_png
    )
    Image.fromarray((doe_level * 17).astype(np.uint8)).save(level_png)
    for bit, mask in zip(bits, layer_masks):
        Image.fromarray((mask * 255).astype(np.uint8)).save(
            args.output / f"mask_bit_{bit}.png"
        )
    save_preview(doe_phase, doe_level, layer_masks, args.size_mm, preview_path)

    gds_rectangle_counts = None
    if not args.skip_gds:
        gds_rectangle_counts = write_gds(
            doe_level, args.size_mm, args.etch_step_nm, gds_path
        )

    unique, counts = np.unique(doe_level, return_counts=True)
    source_resolved = args.source.resolve()
    try:
        source_record = source_resolved.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        source_record = str(source_resolved)

    metadata = {
        "source_phase": source_record,
        "source_sha256": sha256(source_resolved),
        "source_shape": list(phase_source.shape),
        "source_dx_um_per_pixel": source_pitch_m * 1e6,
        "source_full_grid_size_mm": phase_source.shape[1] * source_pitch_m * 1e3,
        "source_span_for_15mm_pixels": args.size_mm * 1e-3 / source_pitch_m,
        "conversion": "centered physical-coordinate sampling of exp(i*phase), cubic interpolation",
        "installation_shift_applied": False,
        "blaze_applied": False,
        "circular_mask_applied": False,
        "doe_shape": "square",
        "doe_size_mm": [args.size_mm, args.size_mm],
        "output_shape": [args.pixels, args.pixels],
        "output_pixel_pitch_um": args.size_mm * 1000.0 / args.pixels,
        "phase_range_rad": [float(doe_phase.min()), float(doe_phase.max())],
        "quantization": {
            "method": "same floor rule as kalyout_doemake.m",
            "wavelength_nm": args.wavelength_nm,
            "substrate_index": args.substrate_index,
            "ambient_index": args.ambient_index,
            "full_2pi_depth_nm": args.wavelength_nm
            / (args.substrate_index - args.ambient_index),
            "etch_step_nm": args.etch_step_nm,
            "binary_mask_bits": list(bits),
            "level_histogram": {str(int(k)): int(v) for k, v in zip(unique, counts)},
            "horizontal_run_rectangles_per_layer": {
                str(bit): count_run_rectangles(mask)
                for bit, mask in zip(bits, layer_masks)
            },
            "gds_rectangles_per_layer": gds_rectangle_counts,
            "gds_layers": {
                "1": {"bit": 8, "etch_depth_increment_nm": 620.0},
                "2": {"bit": 4, "etch_depth_increment_nm": 310.0},
                "3": {"bit": 2, "etch_depth_increment_nm": 155.0},
                "4": {"bit": 1, "etch_depth_increment_nm": 77.5},
                "100": {"purpose": "15 mm square reference outline; not an etch mask"},
            },
        },
        "files": {},
    }
    for path in [npy_path, mat_path, level_path, phase_png, level_png, preview_path]:
        metadata["files"][path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    for bit in bits:
        path = args.output / f"mask_bit_{bit}.png"
        metadata["files"][path.name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    if gds_path.exists():
        metadata["files"][gds_path.name] = {
            "sha256": sha256(gds_path),
            "bytes": gds_path.stat().st_size,
        }
    with (args.output / "metadata.json").open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, ensure_ascii=False, indent=2)

    print(json.dumps(metadata, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
