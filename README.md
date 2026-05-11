# Point2P DOE Baseline Modules

This `main` branch keeps the reusable baseline modules only:

- `initial_phase_generation/` — standalone Romero-Dickey initial phase generator.
- `result_diagnostics/` — standalone focal-plane result diagnostics.

The larger historical MATLAB DOE/MRAF experiment project was moved off `main` to:

```text
codex/doe-romero-dickey-matlab
```

## Baseline Physical Parameters

These values are the current reference parameters for generating the baseline Romero-Dickey initial phase.

| Quantity | Variable | Value | Notes |
|---|---:|---:|---|
| Wavelength | `lambda_m` | `532e-9 m` | 532 nm |
| Physical focal length | `f_m` | `429e-3 m` | 429 mm lens |
| Clear aperture diameter | `aperture_diameter_m` | `15e-3 m` | clear aperture / pupil, not beam diameter |
| Clear aperture radius | `aperture_radius_m` | `7.5e-3 m` | half of clear aperture diameter |
| Input Gaussian 1/e² intensity diameter | `input_1e2_diameter_m` | `5e-3 m` | illuminated beam diameter |
| Input Gaussian 1/e² intensity radius | `input_1e2_radius_m` | `2.5e-3 m` | half of 5 mm |
| Input Gaussian 1/e amplitude radius | `input_1e_radius_m` | `1.76776695297e-3 m` | `input_1e2_radius_m / sqrt(2)` |
| Target x size | `target_size_x_m` | `330e-6 m` | 330 um flat-top reference size |
| Target y size | `target_size_y_m` | `120e-6 m` | 120 um flat-top reference size |
| RD output x scale | `Ro_x_m` | `1.86182562571e-4 m` | `target_size_x_m / sqrt(pi)` |
| RD output y scale | `Ro_y_m` | `6.77027500257e-5 m` | `target_size_y_m / sqrt(pi)` |
| Grid size | `N` | `2048` | square DOE/focal array |
| Requested focal sampling | `requested_focal_dx_m` | `2.5e-6 m/pixel` | 2.5 um/pixel |
| DOE grid extent | `doe_grid_extent_m` | `0.0912912 m` | computed from `lambda*f/dx_focal` |
| DOE sampling | `dx_doe_m` | `4.457578125e-5 m/pixel` | `doe_grid_extent_m / N` |
| Actual focal sampling | `focal_dx_m` | `2.5e-6 m/pixel` | derived FFT sampling |
| Phase method | `phase_method` | `romero_dickey_separable` | x/y separable analytical RD phase |
| Phase sign | `phase_sign` | `1` | Fourier sign convention |
| Phase x/y scale | `phase_scale_x/y` | `1 / 1` | no extra scaling |
| RD beta x | `beta_x` | `9.060975545` | reference baseline beta |
| RD beta y | `beta_y` | `3.294900198` | y direction is the limiting direction |

## Baseline Artifacts

Current verified baseline outputs are recorded in `BASELINES.md`.

Important baseline phase file:

```matlab
load('E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat')
```

Corresponding diagnostics:

```text
E:\program\Point2P\result_diagnostics\artifacts\20260428-141916
```

## Run Initial Phase Generation

```matlab
cd E:\program\Point2P\initial_phase_generation
run_initial_phase_generation
```

This generates:

- `phase0.mat`
- `phase0.png`
- `initial_intensity.png`
- `initial_x_profile.png`
- `initial_y_profile.png`
- `config_snapshot.mat`

## Run Result Diagnostics

```matlab
cd E:\program\Point2P\result_diagnostics
run_diagnostics_example
```

This computes center profiles, 90%/50%/13.5% sizes, transition width, RMS/PV, shoulder, and side-lobe metrics for an initial-phase focal-plane result.

## Standard Pipeline (All Focal Lengths)

The complete DOE design pipeline follows three stages, implemented independently for each lab test folder.

### Lab Test Folders

| Folder | f | dx_doe | DOE extent | Notes |
|--------|---|--------|------------|-------|
| `lab_test_f100mm/` | 100 mm | 10.4 μm | 21.3 mm | β_x=39, β_y=14 — deep geometric |
| `lab_test_f200mm/` | 200 mm | 20.8 μm | 42.6 mm | β_x≈19, β_y≈7 — intermediate |
| `lab_test_f300mm/` | 300 mm | 31.2 μm | 63.8 mm | β_x=13, β_y≈5 — intermediate |
| `rtad_mraf_gs_python/` | 429 mm | 44.6 μm | 91.3 mm | β_x=9, β_y=3 — diffraction |

All folders share: λ=532nm, target=330×120μm, N=2048, focal_dx=2.5μm, clear aperture=15mm.

### Stage 1: Initial Phase Generation

Generate the Romero-Dickey analytical phase for a given beam diameter.

**MATLAB** (native):
```matlab
cd E:\program\Point2P\lab_test_f200mm\matlab
run_initial_phase_generation    % uses default_initial_phase_config.m
```

**Python** (no MATLAB required):
```bash
cd lab_test_f200mm
python make_phase0.py --beam 7 --out make_phase0_output/
```

Output: `phase0.mat` containing `phase0_wrapped_rad`(2048×2048, [0,2π)).

### Stage 2: WGS Refinement

**Standard settings** (validated across f=100/200/300/429mm):
```bash
conda activate slmrtad
cd lab_test_f200mm
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

Key settings rationale (from WGS parameter sweep trials):
- `method=wgs`: direct WGS without MRAF warmup gives the best flat uniformity
- `wgs-strategy=flat_local`: 2D local weight feedback on the flat core
- `bg-factor=0.9`: mild background attenuation (0.05 = aggressive)
- `wgs-feedback-exponent=0.8`: stronger amplitude correction per iteration
- `wgs-weight-max=1.5`: practical sweet spot; max≥2 gave no further benefit
- `--no-swap-phase-xy`: needed when phase0 is generated by Python (numpy convention)

Output per beam diameter:
```
artifacts/<timestamp>_rtad_mraf_gs_truncI0135/
  phase_refined.npy              # 2048×2048 float32 [0, 2π)
  phase_refined.mat              # MATLAB format
  reconstruction_refined.npy     # focal-plane intensity
  report.txt / metrics.csv       # summary & per-iteration metrics
  diagnostics_python/            # detailed focal-plane diagnostics
```

### Stage 3: SLM Phase Conversion

Convert the 2048×2048 computational phase to a 1920×1080 8-bit grayscale PNG matching the SLM pixel grid (6.4 μm pitch).

```bash
python save_slm_phase.py artifacts/<ts>_.../phase_refined.npy \
    --out artifacts/<ts>_.../slm_phase_f200mm_d7mm.png
```

Conversion details: crop central 12.288×6.912mm (SLM active area) → cubic-spline resample to 1920×1080 → quantize [0,2π) to [0,255].

### Batch Mode

Run all three beam diameters (5/6/7mm) at once:

```matlab
% MATLAB
cd matlab && run_all_beam_diameters
```

```bash
# Python
python run_all_beam_diameters.py \
    --phase-mat-5mm <path5> --phase-mat-6mm <path6> --phase-mat-7mm <path7> \
    -- --iters 200
```

### Results Summary

| f | Beam | RMS | size50_x | size50_y | Efficiency |
|----|------|-----|----------|----------|------------|
| 100mm | 5mm | 0.27% | 330.8 | 116.9 | 99.1% |
| 200mm | 7mm | **0.18%** | 327.1 | 118.3 | 97.6% |
| 300mm | 6mm | 1.22% | 329.3 | 118.1 | 94.9% |
| 429mm | 5mm | 1.86% | ~330 | ~124 | 92.5% |

Lower RMS at shorter f is expected — larger β values push the system toward the geometric-optics regime where the Romero-Dickey phase is more accurate.

## RTAD MRAF/GS Python Trial Notes

The independent Python RTAD + MRAF/GS/WGS program lives in:

```text
E:\program\Point2P\rtad_mraf_gs_python
```

Detailed implementation notes and trial logs are in:

```text
rtad_mraf_gs_python\README.md
```

Recent WGS observations from the current `phase0_wrapped_rad` case:

- Direct flat-local WGS without an MRAF warmup gave the best platform
  uniformity so far.
- With direct WGS and `wgs_weight_max > 1`, the x-direction center profile
  changed only mildly as `wgs_weight_max` was increased. `max=1.0` was too
  restrictive, but `max=1.5`, `2.0`, `2.5`, and `3.0` were broadly similar in
  x. The current practical setting is `wgs_weight_max = 1.5`.
- In y, changing `wgs_weight_max` also did not radically change the center
  profile shape, but `size50_y` gradually increased and the shoulder became
  higher as max was increased.
- With `wgs_weight_max = 1.5`, changing `wgs_weight_min` from `0.5` to `0.7`
  and `0.8` had almost no useful effect on y. `min=0.8` began to limit the
  correction and slightly worsened flat uniformity. Keep `wgs_weight_min = 0.5`
  for now.

Summary figures:

```text
rtad_mraf_gs_python\artifacts\wgs_weight_max_summary\center_profiles_compare_wgs_weight_max_1p0_1p5_2p0_2p5_3p0.png
rtad_mraf_gs_python\artifacts\wgs_weight_min_summary\center_profiles_compare_wgs_weight_min_0p5_0p7_0p8_max1p5.png
```

---

## 2026-05-06 Session — Real-World Tolerance Sweeps & Error Fitting

### Divergence Sweeps

Ran stress-level divergence sweep (`--sweep divergence --profile stress`, ±0.5 mrad)
and a fine sweep (±0.01 mrad, ±0.005 mrad).  The convergent (negative) direction
is significantly more sensitive than the divergent direction.  At ±0.005 mrad the
RMS non-uniformity already doubles from the nominal 1.86 %.

All divergence sweep outputs are archived in
`real_world_simulation\artifacts\divergence_sweeps\`.

### Sweep Output Enhancements

- Modified `run_real_world_sweep.py` and `src/plotting.py` so that every sweep
  case saves an individual focal-plane intensity image (not just nominal+worst).
- All images are included in the per-sweep PDF report.
- Center-profile overlay legends are moved below the figure for sweeps with more
  than 10 cases.

### Error Fitting — Five Observed Spot Problems

Identified root causes for five typical flat-top spot errors by running
single-parameter and compound sweeps across the 7 error types (defocus,
beam offset, beam size, divergence, pointing, aperture, ellipticity):

| # | Problem | Root Cause |
|---|---------|-----------|
| 1 | Long-edge concavity (长边中间内凹) | Ellipticity Dx > Dy |
| 2 | Long-edge energy loss + short-edge tilt | beam_offset_x + beam_offset_y |
| 3 | All four edges concave (四条边都内凹) | Negative defocus + beam size ≥ 5.5 mm |
| 4 | Energy piling at top/bottom | beam_size (already known) |
| 5 | Wrong aspect ratio | Pending further analysis |

All fitting results are documented in:
- `real_world_simulation\artifacts\error_fitting_summary.md`
- `presentation\error_fitting_summary.md`

### Presentation

Created a Beamer PDF presentation (24 pages, XeLaTeX) summarizing the entire
project: Romero-Dickey theory, RTAD target, WGS optimization, tolerance
assessment, and the five error-fitting results.  Source and PDF are in:

```text
presentation\
    Point2P_Project_Report.tex
    Point2P_Project_Report.pdf         (24 pages)
    Point2P_Project_Summary.md
    Presentation_Outline.md
    error_fitting_summary.md
    RTAD\          (baseline images)
    PROBLEM1\      (ellipticity comparison)
    PROBLEM2\      (beam offset)
    PROBLEM3\      (defocus + beam size)
    BEAMSIZE\      (beam size worst case)
    OVERLAY\       (profile overlays)
```

### LaTeX Environment

TeX Live 2025 was installed locally at `C:\texlive\2025\`. The presentation
compiles with `xelatex Point2P_Project_Report.tex` from the `presentation\`
directory.  It can also be compiled on Overleaf (upload the `.tex` file and
the image subfolders).
