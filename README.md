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
