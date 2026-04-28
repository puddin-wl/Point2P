# RTAD MRAF/GS Python Refinement

This folder contains a new, independent Python implementation for refining an
existing Point2P initial phase (`phase0`) into a rectangular RTAD flat-top DOE
result. It does not modify the older `initial_phase_generation`,
`result_diagnostics`, or MATLAB flows.

## Purpose

The program:

- loads an existing point-to-point initial phase in radians
- builds a raised-cosine rectangular RTAD target
- runs GS, MRAF, WGS, or MRAF-then-WGS style refinement
- saves refined phase, reconstructed intensity, target/mask files, plots,
  metrics, and a text report

The current physical defaults are:

- wavelength: 532 nm
- focal length: 429 mm
- input Gaussian 1/e^2 intensity diameter: 5 mm
- clear aperture: 15 mm
- target intensity size50: 330 um x 120 um

## Relation To Point2P Phase0

`phase0` is treated as the initial DOE phase in radians. It is wrapped to
`[0, 2*pi)` after loading. No resizing is performed. If the phase shape does not
match the generated target shape, the run stops.

The input amplitude is regenerated from the project physical parameters:

```text
I_input(r) = exp(-2 r^2 / w^2)
A_input(r) = exp(-r^2 / w^2)
w = 2.5 mm
```

A circular 15 mm clear aperture is then applied, and the amplitude is L2
normalized.

## Relation To slmsuite

The implementation was written after reading slmsuite 0.4.1, especially:

- `holography/algorithms/_hologram.py`
- `holography/algorithms/_feedback.py`
- `holography/toolbox`
- `docs/source/examples.py` and `examples.rst`

Key semantics copied conceptually, not as source code:

- target arrays are target amplitude, not intensity
- input amplitude and target amplitude are L2 normalized
- FFT propagation uses centered FFTs with `norm="ortho"`
- MRAF free/noise regions are represented in slmsuite by `NaN` target pixels
- slmsuite `mraf_factor=0` fully attenuates a noise region, while
  `mraf_factor=1` leaves it unchanged
- WGS feedback increases weights where computed amplitude is weak and decreases
  weights where it is strong, then normalizes weights again

This local version keeps explicit masks instead of using `NaN` in the target.

## RTAD Target

The target is defined first as intensity. The full raised-cosine template is:

```text
I_full(x,y) = C(abs(x); a0, a1) * C(abs(y); b0, b1)
A_full = sqrt(I_full)
```

The intensity size50 remains fixed to 330 um x 120 um. With `a50 = W50/2`
and `b50 = H50/2`, the separable raised-cosine edge is:

```text
a0 = a50 - delta_x_um
a1 = a50 + delta_x_um
b0 = b50 - delta_y_um
b1 = b50 + delta_y_um

C(u; u0, u1) = 1                                      for u <= u0
             = 0.5 * (1 + cos(pi*(u-u0)/(u1-u0)))     for u0 < u < u1
             = 0                                      for u >= u1

```

The current default constraint mode is truncated RTAD:

```text
release_level = exp(-2) = 0.1353352832366127
mask_flat = abs(x) <= a0 and abs(y) <= b0
mask_signal = (I_full >= release_level) OR mask_flat
mask_edge_lock = mask_signal AND NOT mask_flat
mask_template_support = I_full > 0
mask_free = guard_window AND NOT mask_signal
mask_bg_far = NOT (mask_signal OR mask_free)
```

The old target constrained every `I_full > 0` pixel. The new target only
constrains `mask_signal`. The low-intensity tail between `release_level` and 0
is released into the MRAF free/noise region. This avoids forcing the DOE result
to follow a full mathematical tail all the way to zero, which can over-smooth
the edge and suppress physically useful sidelobe/halo structure.

`mask_support` is kept as a compatibility alias, but now means constrained
signal support, not full template support. Flat uniformity metrics use only the
fixed `mask_flat`, not the edge or free region.

## Refinement Formula

GS projection:

```text
E_far = FFT(A_input * exp(i phase))
E_far' = A_target * exp(i angle(E_far))
phase = angle(IFFT(E_far'))
```

MRAF projection in this implementation:

```text
mask_signal: E_far' = W * exp(i angle(E_far))
mask_free:   E_far' = mraf_factor * E_far
mask_bg_far: keep, attenuate by bg_factor, or zero according to bg_mode
```

Default `bg_mode="attenuate"` and `bg_factor=0.05` weakly constrain the far
background without hard-zeroing everything outside the support. This is
intentional for RTAD edge behavior.

`mraf_then_wgs` is now the default refinement mode. It first runs MRAF, then
keeps the same truncated-RTAD signal/free/background projection while updating
only flat-core target weights:

```text
MRAF stage:
  mask_signal: E_far' = A_target * exp(i angle(E_far))
  mask_free:   E_far' = mraf_factor * E_far
  mask_bg_far: E_far' = bg_factor * E_far

WGS stage:
  target_amp_eff = target_amp
  target_amp_eff[mask_flat] = target_amp[mask_flat] * weights[mask_flat]

  every wgs_update_every WGS iterations:
    amp_mean = mean(abs(E_far)[mask_flat])
    weights[mask_flat] *= (amp_mean / (abs(E_far)[mask_flat] + eps))^wgs_feedback_exponent
    weights[mask_flat] = clip(weights[mask_flat], wgs_weight_min, wgs_weight_max)
    weights[mask_flat] /= mean(weights[mask_flat])
```

WGS does not update `mask_edge_lock`, `mask_free`, or `mask_bg_far`. It is only
a flat-core uniformity feedback; it does not use size50, e^-2 efficiency, or
diagnostics-side metrics as constraints.

Two WGS strategies are available:

```text
flat_local:
  Run all WGS iterations with the 2D local flat-core update above.

xy_then_x:
  Run WGS-XY first with the same 2D local flat-core update.
  Then freeze that y-local structure and update a separate 1D x weight:

    amp_x[j] = mean(abs(E_far)[mask_flat[:, j]])
    amp_x_mean = mean(amp_x over valid flat-core x columns)
    w_x[j] *= (amp_x_mean / (amp_x[j] + eps))^wgs_x_feedback_exponent
    w_x = clip(w_x, wgs_x_weight_min, wgs_x_weight_max)
    w_x /= mean(w_x)

  During x-only projection:
    target_amp_eff[mask_flat] =
      target_amp[mask_flat] * weights_2d[mask_flat] * w_x_broadcast[mask_flat]
```

`xy_then_x` is meant to avoid continuing to push the y direction after the
flat core has become reasonably uniform. The second stage mainly changes the x
profile because the x weight is broadcast along y and is used only inside
`mask_flat`.

## How To Run

From this folder:

```powershell
cd E:\program\Point2P\rtad_mraf_gs_python
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --phase-mat "E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat" --phase-var phase0_wrapped_rad --iters 200 --method mraf --mraf-factor 0.4
```

Equivalent explicit truncated-target run:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --phase-mat "E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat" --phase-var phase0_wrapped_rad --iters 200 --method mraf --mraf-factor 0.4 --constraint-mode truncated_rtad --release-level 0.1353352832366127 --bg-factor 0.05
```

Recommended MRAF followed by XY then x-only WGS:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --phase-mat "E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat" --phase-var phase0_wrapped_rad --method mraf_then_wgs --mraf-iters 150 --wgs-strategy xy_then_x --wgs-xy-iters 20 --wgs-xonly-iters 30 --mraf-factor 0.4 --wgs-xy-feedback-exponent 0.3 --wgs-x-feedback-exponent 0.45 --release-level 0.1353352832366127 --bg-factor 0.05
```

By default this program applies `swap_phase_xy = True` to imported `phase0`
files, because the current MATLAB v7.3 `phase0` orientation appears x/y swapped
when read into Python. Use `--no-swap-phase-xy` only if a future input file has
already been corrected.

Smoke test without a phase file:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --iters 20 --smoke-shape 256
```

Useful overrides:

- `--phase-mat`
- `--phase-var`
- `--iters`
- `--method`
- `--mraf-factor`
- `--constraint-mode`
- `--release-level`
- `--delta-x`
- `--delta-y`
- `--bg-factor`
- `--wgs-strategy`
- `--wgs-xy-iters`
- `--wgs-xonly-iters`
- `--wgs-xy-feedback-exponent`
- `--wgs-x-feedback-exponent`
- `--wgs-x-weight-min`
- `--wgs-x-weight-max`
- `--outdir`
- `--use-cupy` / `--no-cupy`
- `--swap-phase-xy` / `--no-swap-phase-xy`

## Outputs

Each run creates a timestamped directory under `artifacts/`, unless `--outdir`
is supplied. Saved files include:

- `config_used.json`
- `target.mat`
- `target.npz`
- `phase_refined.mat`
- `phase_refined.npy`
- `phase_after_mraf.npy`
- `phase_after_wgs_xy.npy`
- `reconstruction_refined.npy`
- `reconstruction_after_mraf.npy`
- `reconstruction_after_wgs_xy.npy`
- `wgs_weights_final.npy`
- `wgs_weights_2d_final.npy`
- `wgs_x_weights_final.npy`
- `metrics.csv`
- `report.txt`
- `target_full_intensity.png`
- `target_constraint_intensity.png`
- `target_profiles.png`
- `target_intensity.png`
- `target_amplitude.png`
- `masks.png`
- `initial_reconstruction_intensity.png`
- `refined_reconstruction_intensity.png`
- `center_profiles_compare.png`
- `center_profiles_compare_initial_mraf_wgs.png`
- `center_profiles_compare_mraf_xy_xonly.png`
- `wgs_weights_final.png`
- `wgs_weight_hist.png`
- `wgs_x_weights_final.png`
- `wgs_x_weights_profile.png`
- `phase0.png`
- `phase_refined.png`
- `convergence_metrics.png`
- `diagnostics_python/diagnostics_report.txt`
- `diagnostics_python/diagnostics_metrics.csv`
- `diagnostics_python/diagnostics_metrics.json`
- `diagnostics_python/center_profiles_diagnostics.png`
- `diagnostics_python/intensity_diagnostics.png`
- `diagnostics_python/derivative_sidelobe_diagnostics.png`

The lightweight Python diagnostics can also be run on an existing case:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_diagnostics_case.py "E:\program\Point2P\rtad_mraf_gs_python\artifacts\20260428-164830_rtad_mraf_gs"
```

## Recent Trial Notes

These notes summarize the parameter trials run on 2026-04-28. They are
empirical observations for the current `phase0_wrapped_rad` case, not general
proofs.

Direct flat-local WGS, no MRAF warmup:

```text
method = wgs
wgs_strategy = flat_local
iters = 200
mraf_factor = 0.4
bg_factor = 0.05
wgs_feedback_exponent = 0.3
wgs_update_every = 5
```

This direct WGS mode produced the best flat-core uniformity so far. With
`wgs_weight_min=0.5` and `wgs_weight_max=2.0`, the result reached about
`uniformity_rms_percent = 96.19%` with `size50_x/y = 330.48 / 121.08 um`
(`artifacts/20260428-221535_rtad_mraf_gs_truncI0135`). The `13.5%` size and
efficiency changes were small enough for the current platform-focused search.

WGS max sweep with direct WGS:

```text
wgs_weight_min = 0.5
wgs_weight_max = 1.0, 1.5, 2.0, 2.5, 3.0
```

Artifacts:

```text
max=1.0: artifacts/20260428-225058_rtad_mraf_gs_truncI0135
max=1.5: artifacts/20260428-225032_rtad_mraf_gs_truncI0135
max=2.0: artifacts/20260428-221535_rtad_mraf_gs_truncI0135
max=2.5: artifacts/20260428-224354_rtad_mraf_gs_truncI0135
max=3.0: artifacts/20260428-224415_rtad_mraf_gs_truncI0135
```

Observed trend: `max=1.0` is too restrictive and does not repair the platform.
Once `max > 1`, the x-direction center profile changes only mildly as max is
increased. In y, the profile shape also changes modestly, while `size50_y`
gradually grows and the shoulder becomes higher. For the current round we chose
to fix `wgs_weight_max = 1.5` as the less aggressive platform-WGS setting.

The combined center-profile figure is:

```text
artifacts/wgs_weight_max_summary/center_profiles_compare_wgs_weight_max_1p0_1p5_2p0_2p5_3p0.png
```

WGS min sweep with direct WGS and `max=1.5`:

```text
wgs_weight_max = 1.5
wgs_weight_min = 0.5, 0.7, 0.8
```

Artifacts:

```text
min=0.5: artifacts/20260428-225032_rtad_mraf_gs_truncI0135
min=0.7: artifacts/20260428-225817_rtad_mraf_gs_truncI0135
min=0.8: artifacts/20260428-225843_rtad_mraf_gs_truncI0135
```

Observed trend: `wgs_weight_min` has almost no visible effect on the y direction
for this case. `min=0.5` and `min=0.7` are nearly identical; `min=0.8` starts to
limit the correction and slightly worsens flat uniformity. The current practical
choice is to keep `wgs_weight_min = 0.5`.

The combined center-profile figure is:

```text
artifacts/wgs_weight_min_summary/center_profiles_compare_wgs_weight_min_0p5_0p7_0p8_max1p5.png
```

Earlier MRAF/WGS trials:

```text
MRAF only, truncated RTAD release=exp(-2):
  artifacts/20260428-173456_rtad_mraf_gs_truncI0135

150 MRAF + 50 flat-local WGS:
  artifacts/20260428-180526_rtad_mraf_gs_mrafwgs_truncI0135
  uniformity improved from about 80.51% after MRAF to about 93.71%.

150 MRAF + 20 XY-WGS + 30 X-only:
  artifacts/20260428-183223_rtad_mraf_gs_mrafwgs_xyx_truncI0135
  froze y too early; y size/profile was worse than the 50-step flat-local WGS.

150 MRAF + 40 XY-WGS + 10 X-only:
  artifacts/20260428-215516_rtad_mraf_gs_mrafwgs_xyx_truncI0135
  recovered most of y while keeping a small x-only correction.

150 MRAF + 50 XY-WGS + 20 X-only:
  artifacts/20260428-215858_rtad_mraf_gs_mrafwgs_xyx_truncI0135
  gave better x profile and high uniformity, but the current focus shifted back
  to simpler direct flat-local WGS without x-only correction.
```

## First Recommended Parameters

```text
method = mraf_then_wgs
mraf_iters = 150
wgs_iters = 50
wgs_strategy = xy_then_x
wgs_xy_iters = 20
wgs_xonly_iters = 30
mraf_factor = 0.4
constraint_mode = truncated_rtad
release_level = 0.1353352832366127
bg_factor = 0.05
wgs_update_mask = flat
wgs_feedback = amplitude
wgs_xy_feedback_exponent = 0.3
wgs_x_feedback_exponent = 0.45
wgs_xy_update_every = 5
wgs_x_update_every = 5
wgs_xy_weight_min = 0.5
wgs_xy_weight_max = 2.0
wgs_x_weight_min = 0.5
wgs_x_weight_max = 2.5
delta_x_um = 15
delta_y_um = 8
guard_x_um = 20
guard_y_um = 12
```

The internal metrics are quick trend checks only. Final acceptance should still
be done by the existing Point2P diagnostics once connected.
