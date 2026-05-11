# lab_test_f200mm — f=200mm DOE Beam Shaping Pipeline

Wavelength 532nm, focal length 200mm, 3 beam diameters (5/6/7mm), target flat-top 330×120μm.

This folder is the **reference implementation** for the standard DOE design pipeline.
To adapt for a new focal length, copy this folder and change 2 numbers in 2 config files.

---

## Pipeline Overview

```
Stage 1                  Stage 2                  Stage 3
MATLAB or Python         Python (GPU)             Python
    │                        │                       │
    ▼                        ▼                       ▼
RD initial phase  ──→  WGS refinement  ──→  SLM phase conversion
phase0.mat             phase_refined.npy     1920×1080 8-bit PNG
(2048×2048)            (2048×2048)           (SLM native res)
```

| Stage | Input | Output | Description |
|-------|-------|--------|-------------|
| 1 | Physical params | `phase0.mat` | Romero-Dickey separable initial phase |
| 2 | `phase0.mat` + beam diameter | `phase_refined.npy` | WGS iterative flat-top optimization |
| 3 | `phase_refined.npy` | `slm_phase_f200mm_d*.png` | Crop + complex-field resample to SLM grid |

---

## File Architecture

```
lab_test_f200mm/
├── README.md                          # This file
├── .gitignore
│
├── matlab/                            # Stage 1: MATLAB initial phase (also works in Python)
│   ├── default_initial_phase_config.m # ← EDIT for new focal length (f_m, beam diameter)
│   ├── generate_initial_phase.m       # Core RD phase generator (no changes needed)
│   ├── run_initial_phase_generation.m # Single-beam entry point
│   └── run_all_beam_diameters.m       # Batch: 5/6/7mm
│
├── config_default.py                  # ← EDIT for new focal length (focal_length_m, beam diameter)
├── make_phase0.py                     # Stage 1 alt: Python RD phase generator (no MATLAB)
├── run_rtad_mraf_gs_case.py           # Stage 2: WGS refinement (accepts --beam-diameter)
├── run_all_beam_diameters.py          # Stage 2 batch: 5/6/7mm
├── run_diagnostics_case.py            # Standalone diagnostics re-run
├── save_slm_phase.py                  # Stage 3: SLM pixel conversion
│
├── src/                               # Shared library (identical across all f folders)
│   ├── mraf_gs.py                     # WGS/MRAF iterative refinement engine
│   ├── propagation.py                 # Fourier-lens forward/backward propagation
│   ├── rtad_target.py                 # RTAD raised-cosine flat-top target builder
│   ├── diagnostics.py                 # Focal-plane metrics computation
│   ├── plotting.py                    # Intensity/phase/profile visualization
│   ├── metrics.py                     # Low-level width/uniformity/efficiency metrics
│   ├── io_mat.py                      # MATLAB .mat ↔ numpy I/O
│   ├── backend.py                     # CuPy/NumPy backend abstraction
│   └── utils.py                       # Timestamps, JSON, CSV helpers
│
└── artifacts/                         # Runtime outputs (gitignored)
    └── <timestamp>_rtad_mraf_gs_.../
        ├── phase_refined.npy/.mat     # Final DOE phase (2048×2048 float32, [0, 2π))
        ├── reconstruction_refined.npy # Focal-plane intensity
        ├── slm_phase_f200mm_d*.png    # SLM-loadable 8-bit PNG (1920×1080)
        ├── target.npz/.mat            # RTAD target masks & profiles
        ├── report.txt                 # Human-readable summary with all metrics
        ├── metrics.csv                # Per-iteration convergence log
        ├── config_used.json           # Effective config with CLI overrides applied
        └── diagnostics_python/        # Detailed diagnostic plots & metrics JSON
```

---

## Physical Parameters

### Constants (same for all focal lengths)

| Symbol | Value | Description |
|--------|-------|-------------|
| λ | 532 nm | Laser wavelength |
| Clear aperture | 15 mm diameter | DOE physical aperture |
| N | 2048 | Computational grid (square) |
| focal_dx | 2.5 μm/pixel | Focal-plane sampling |
| Target W50 × H50 | 330 × 120 μm | Flat-top 50%-intensity size |
| Target delta x/y | 15 / 8 μm | Raised-cosine edge half-width |
| Target guard x/y | 20 / 12 μm | Guard band beyond edge |
| SLM resolution | 1920 × 1080 | SLM native pixel count |
| SLM pixel pitch | 6.4 μm | SLM physical pixel size |

### f=200mm Derived Parameters

| Parameter | Formula | Value |
|-----------|---------|-------|
| DOE grid extent | λf / focal_dx | 42.56 mm |
| Computational pixel dx_doe | extent / N | 20.78 μm |
| SLM active area | 1920×6.4μm × 1080×6.4μm | 12.288 × 6.912 mm |
| SLM crop (computational px) | area / dx_doe | 591 × 332 px |
| SLM resample ratio | SLM px / crop px | 3.25× (x), 3.25× (y) |

### β Values (Romero-Dickey quality parameter)

β = 2π · r₁ₑ · Rₒ / (λf)  —  larger β → stationary-phase approximation more accurate.

β > 10: geometric-optics regime (RD phase nearly exact).
β < 10: diffraction-influenced (WGS refinement provides meaningful improvement).

| Beam D (1/e²) | r₁ₑ (mm) | β_x | β_y | x regime | y regime |
|---------------|-----------|-----|-----|----------|----------|
| 5 mm | 1.77 | 19.4 | 7.1 | geometric | diffraction |
| 6 mm | 2.12 | 23.3 | 8.5 | deep geometric | intermediate |
| 7 mm | 2.47 | 27.2 | 9.9 | deep geometric | near threshold |

---

## Stage 1: Initial Phase Generation

The Romero-Dickey phase is an analytical solution to the stationary-phase energy
conservation problem for a Gaussian → rectangular flat-top mapping.

### Option A: MATLAB (original)

```matlab
cd E:\program\Point2P\lab_test_f200mm\matlab

% Single beam:
run_initial_phase_generation

% All three (5/6/7mm):
run_all_beam_diameters
```

Output per beam: `matlab/artifacts/<timestamp>/phase0.mat` containing:
- `phase0_wrapped_rad`: 2048×2048 double, [0, 2π), NaN outside aperture
- `phase0_unwrapped_rad`: continuous phase
- `focal_x_m`, `focal_y_m`, `x_m`, `y_m`: coordinate axes
- `phase_info`: beta values, method, phase range

### Option B: Python (no MATLAB required)

```bash
python make_phase0.py --beam 5 --out make_phase0_output_5mm/
python make_phase0.py --beam 6 --out make_phase0_output_6mm/
python make_phase0.py --beam 7 --out make_phase0_output_7mm/
```

This replicates the MATLAB logic exactly: coordinate grid → Gaussian field →
Romero-Dickey 1D phases → separable 2D phase → wrap to [0, 2π) → save as .mat.

---

## Stage 2: WGS Refinement

### Standard Command

```bash
conda activate slmrtad
cd E:\program\Point2P\lab_test_f200mm

python run_rtad_mraf_gs_case.py \
    --phase-mat make_phase0_output/phase0.mat \
    --phase-var phase0_wrapped_rad \
    --beam-diameter 7 \
    --method wgs \
    --wgs-strategy flat_local \
    --iters 200 \
    --wgs-feedback-exponent 0.8 \
    --wgs-weight-min 0.5 \
    --wgs-weight-max 1.5 \
    --bg-factor 0.9 \
    --no-swap-phase-xy
```

### Key Parameter Rationale

| Parameter | Value | Why |
|-----------|-------|-----|
| `--method wgs` | Direct WGS | No MRAF warmup — direct WGS gives best flat uniformity |
| `--wgs-strategy flat_local` | 2D local weights | Weighted feedback per pixel in the flat-core region |
| `--wgs-feedback-exponent 0.8` | Stronger correction | Faster convergence than default 0.3 |
| `--wgs-weight-max 1.5` | Moderate range | `max≥2` showed no further benefit; 1.5 is sweet spot |
| `--wgs-weight-min 0.5` | Default lower bound | `min≥0.7` slightly worsened uniformity |
| `--bg-factor 0.9` | Mild background attenuation | 0.05 (default) is too aggressive for flat_local WGS |
| `--iters 200` | 200 iterations | Sufficient for f=200mm; more gives diminishing returns |
| `--no-swap-phase-xy` | No transpose | Required when phase0 comes from Python (numpy convention) |

### What WGS Does

Weighted Gerchberg-Saxton: after each forward-backward FFT pair, compare the
focal-plane amplitude to the RTAD target. Pixels that are too dim get their
target weight increased; pixels that are too bright get it decreased. The
weights are updated every 5 iterations on the flat-core region only.

The RTAD target is a separable raised-cosine rectangle:
- Flat core: |x| ≤ 150 μm, |y| ≤ 52 μm
- Cosine edge: 150→180 μm (x), 52→68 μm (y)
- Guard band: 180→200 μm (x), 68→80 μm (y)

Truncated constraint mode only enforces the target on pixels where I ≥ e⁻²
(≈13.5%), leaving the low-intensity tail as a free optimization variable.

### GPU Acceleration

The pipeline auto-detects CuPy. With RTX 5070 Ti, 200 iterations take ~3 seconds.
CPU-only mode via `--no-cupy` is available but much slower.

### Batch Mode

```bash
python run_all_beam_diameters.py \
    --phase-mat-5mm <path5>/phase0.mat \
    --phase-mat-6mm <path6>/phase0.mat \
    --phase-mat-7mm <path7>/phase0.mat \
    -- --iters 200
```

All extra arguments after `--` are forwarded to each `run_rtad_mraf_gs_case.py` call.

---

## Stage 3: SLM Phase Conversion

The computational phase lives on a 2048×2048 grid at 20.78 μm/pixel.
The SLM has 1920×1080 pixels at 6.4 μm/pixel. Conversion requires cropping
the central physical region and resampling.

```bash
python save_slm_phase.py artifacts/<ts>_.../phase_refined.npy \
    --out artifacts/<ts>_.../slm_phase_f200mm_d7mm.png
```

### Conversion Math

```
SLM physical area = 1920 × 6.4μm = 12.288 mm  (x)
                  = 1080 × 6.4μm =  6.912 mm  (y)

Computational pixels in this area:
  crop_w = 12.288 mm / 20.78 μm = 591 px
  crop_h =  6.912 mm / 20.78 μm = 332 px

Resample: 591×332 → 1920×1080  (×3.25 in both directions)
Quantize: phase [0, 2π) → grayscale [0, 255]
```

### Critical: Complex-Field Interpolation

**Never interpolate the wrapped phase directly.**

The phase has 2π→0 discontinuities where the physical phase exceeds 2π and wraps.
Cubic spline interpolation across these discontinuities produces ringing artifacts.

Correct approach: interpolate the complex field `exp(iφ)`, which is continuous
everywhere (`exp(i·0) = exp(i·2π) = 1`), then extract the phase.

```python
# WRONG — ringing at phase wraps
phase_slm = zoom(phase_cropped, (zy, zx), order=3)

# CORRECT — continuous across wraps
complex_field = np.exp(1j * phase_cropped)
real_part = zoom(complex_field.real, (zy, zx), order=3)
imag_part = zoom(complex_field.imag, (zy, zx), order=3)
phase_slm = np.arctan2(imag_part, real_part)
```

This was a hard-won lesson from the initial f=200mm run:
the direct interpolation produced 83% of pixels with gradient > 10.
After switching to complex-field interpolation, it dropped to levels
comparable to the validated f=300mm SLM phases.

---

## Results: f=200mm, All Beam Diameters

Refinement: direct WGS, flat_local strategy, 200 iterations, RTX 5070 Ti (~3s).

| Beam | β_x / β_y | RMS | size50_x | size50_y | e⁻² eff. |
|------|-----------|-----|----------|----------|----------|
| 4.5mm | 17.5 / 6.4 | 0.25% | 328.9 | 118.1 | 96.8% |
| 5.0mm | 19.4 / 7.1 | 0.12% | 328.2 | 117.8 | 96.5% |
| 5.5mm | 21.4 / 7.8 | 0.13% | 326.5 | 118.1 | 95.7% |
| 6.0mm | 23.3 / 8.5 | 0.09% | 326.9 | 118.1 | 95.7% |
| 7.0mm | 27.2 / 9.9 | 0.18% | 327.1 | 118.3 | 97.6% |

All beams meet the <2% RMS target. Best uniformity: 6.0mm (0.09%).
Worst: 4.5mm (0.25%) — β_y=6.4 is deepest in the diffraction regime.

Output artifacts per beam:
```
artifacts/20260511-185912_rtad_mraf_gs_truncI0135/  (5.5mm)
artifacts/20260511-190932_rtad_mraf_gs_truncI0135/  (4.5mm)
artifacts/20260511-165932_rtad_mraf_gs_truncI0135/  (5.0mm)
artifacts/20260511-165943_rtad_mraf_gs_truncI0135/  (6.0mm)
artifacts/20260511-160944_rtad_mraf_gs_truncI0135/  (7.0mm)
```

---

## Diagnostics

After refinement, detailed metrics are auto-computed. To re-run diagnostics:

```bash
python run_diagnostics_case.py artifacts/<timestamp>_rtad_mraf_gs_.../
```

Key diagnostic outputs in `diagnostics_python/`:
- `summary.json` — all metrics (size50/13.5/90, RMS, efficiency, transition widths, etc.)
- `center_profiles.png` — x/y center-line profiles with level crossings marked
- `intensity_2d.png` — 2D focal-plane intensity heatmap
- `profile_overlay.png` — comparison of initial vs refined center profiles

---

## How to Adapt for a New Focal Length

1. Copy this entire folder: `cp -r lab_test_f200mm lab_test_fXXXmm`
2. Edit **2 numbers** in `config_default.py`:
   - `"focal_length_m": XXXe-3`
   - `"input_gaussian_1e2_diameter_m": <your beam diameter>`
3. Edit **2 numbers** in `matlab/default_initial_phase_config.m`:
   - `cfg.f_m = XXXe-3`
   - `cfg.input_1e2_diameter_m = <your beam diameter>`
4. Update `save_slm_phase.py` default `--f` value (or always pass `--f XXXe-3` via CLI)
5. Update this README with the new β values and results

Everything else (src/, run scripts, refinement parameters) stays the same.

---

## Dependencies

- **Environment**: `conda activate slmrtad`
- **GPU**: CuPy 13.6 + NVIDIA GPU (tested on RTX 5070 Ti)
- **CPU fallback**: NumPy (pass `--no-cupy`)
- **Core packages**: numpy, scipy, Pillow (for SLM PNG export)
- **Optional**: MATLAB (for Stage 1 Option A; Python `make_phase0.py` can replace it)
