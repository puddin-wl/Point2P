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

