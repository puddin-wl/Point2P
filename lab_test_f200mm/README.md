# lab_test_f200mm — DOE refinement for f=200mm lab setup

Wavelength 532nm, focal length 200mm, input beam 5/6/7mm, target 330×120μm.

Three beam diameters are supported to cover different experimental conditions.
All other parameters (wavelength, target size, grid, refinement algorithm) are
identical to the baseline lab test folders.

## Quick start

### 1. Generate Romero-Dickey initial phase (MATLAB)

```matlab
cd E:\program\Point2P\lab_test_f200mm\matlab

% Single beam (edit default_initial_phase_config.m first):
run_initial_phase_generation

% All three beams (5/6/7 mm):
run_all_beam_diameters
```

This produces `phase0.mat` under `matlab/artifacts/<timestamp>/` for each beam.

### 2. Run DOE refinement (Python)

```bash
conda activate slmrtad
cd E:\program\Point2P\lab_test_f200mm

# Single beam:
python run_rtad_mraf_gs_case.py \
    --phase-mat matlab/artifacts/<timestamp>/phase0.mat \
    --beam-diameter 5 --iters 200

# All three beams:
python run_all_beam_diameters.py \
    --phase-mat-5mm matlab/artifacts/<ts5>/phase0.mat \
    --phase-mat-6mm matlab/artifacts/<ts6>/phase0.mat \
    --phase-mat-7mm matlab/artifacts/<ts7>/phase0.mat \
    -- --iters 200
```

### 3. Convert to SLM-loadable PNG

```bash
python save_slm_phase.py artifacts/<ts5>_.../phase_refined.npy \
    --out artifacts/<ts5>_.../slm_phase_f200mm_d5mm.png
python save_slm_phase.py artifacts/<ts6>_.../phase_refined.npy \
    --out artifacts/<ts6>_.../slm_phase_f200mm_d6mm.png
python save_slm_phase.py artifacts/<ts7>_.../phase_refined.npy \
    --out artifacts/<ts7>_.../slm_phase_f200mm_d7mm.png
```

Output: 1920×1080 8-bit grayscale PNG, phase [0, 2π) mapped to [0, 255].

### 4. Smoke test (no MATLAB needed)

```bash
python run_rtad_mraf_gs_case.py --smoke-shape 2048 --no-cupy --beam-diameter 5
```

## Key parameters (f=200mm vs baseline)

| Parameter | f=100mm | f=200mm (this) | f=300mm | f=429mm (baseline) |
|-----------|---------|----------------|---------|---------------------|
| focal_length_m | 100e-3 | 200e-3 | 300e-3 | 429e-3 |
| dx_doe | 10.4 μm | 20.8 μm | 31.2 μm | 44.6 μm |
| DOE grid extent | 21.3 mm | 42.6 mm | 63.8 mm | 91.3 mm |

All other parameters (λ=532nm, N=2048, focal_dx=2.5μm, target=330×120μm,
clear aperture=15mm) are identical.

## Beta values at f=200mm

| Beam D (1/e²) | β_x | β_y | Regime |
|---------------|-----|-----|--------|
| 5 mm | ~19.4 | ~7.1 | x: geometric; y: diffraction-influenced |
| 6 mm | ~23.3 | ~8.5 | x: deep geometric; y: intermediate |
| 7 mm | ~27.2 | ~9.9 | x: deep geometric; y: near threshold |

β = 2π · r₁ₑ · Rₒ / (λf). Larger β → stationary-phase approximation is
more accurate. β > 10 is the geometric-optics threshold.

At f=200mm with 5mm beam, the x-direction (β_x ≈ 19) is well into the
geometric regime while the y-direction (β_y ≈ 7) still benefits from
WGS refinement. Larger beam diameters push both betas higher.

## Output

Results go to `artifacts/<timestamp>_rtad_mraf_gs_.../`:
- `phase_refined.npy` / `phase_refined.mat` — 2048×2048 float32, [0, 2π)
- `reconstruction_refined.npy` — focal-plane intensity
- `target.npz` / `target.mat` — target masks and parameters
- `report.txt` — summary metrics
- `metrics.csv` — per-iteration metrics
- `config_used.json` — effective config with CLI overrides
- `diagnostics_python/` — detailed focal-plane diagnostics
- `slm_phase_f200mm_d*.png` — SLM-loadable 1920×1080 8-bit PNG (after Step 3)

## SLM conversion details

| Parameter | Value |
|-----------|-------|
| SLM resolution | 1920 × 1080 |
| SLM pixel pitch | 6.4 μm |
| SLM active area | 12.288 × 6.912 mm |
| Cropped computational region | ~591 × 333 px (at dx_doe=20.78μm) |
| Resampling method | cubic spline (scipy.ndimage.zoom, order=3) |
| Output format | 8-bit grayscale PNG, phase/2π × 255 |

The computational phase (2048×2048, dx_doe≈20.78μm, extent≈42.56mm) is
cropped to the central physical region that fits on the SLM, then
resampled to SLM native resolution.
