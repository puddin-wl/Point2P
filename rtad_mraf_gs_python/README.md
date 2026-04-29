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

WGS iteration-count sweep with direct WGS, `min=0.5`, and `max=1.5`:

```text
num_iters = 100, 150, 200, 250, 300
```

Artifacts:

```text
artifacts/wgs_num_iters_20260429-151618
```

Observed trend: after `num_iters=100`, the center profiles change only slightly.
Increasing from 100 to 300 improves `uniformity_rms_percent` only from about
94.45% to 94.71%, while the x/y profile shapes are almost overlapping. The
practical choice is to keep `num_iters = 200` for the next parameter sweeps,
because later iterations mostly add runtime with very small profile changes.

The combined center-profile figure is:

```text
artifacts/wgs_num_iters_20260429-151618/center_profiles_compare_num_iters_100_150_200_250_300.png
```

WGS feedback-exponent sweep with direct WGS, `num_iters=200`, `min=0.5`,
and `max=1.5`:

```text
wgs_feedback_exponent = 0.2, 0.3, 0.4, 0.5
wgs_feedback_exponent = 0.5, 0.6, 0.7, 0.8, 0.9
```

Artifacts:

```text
artifacts/wgs_feedback_exponent_20260429-152232
artifacts/wgs_feedback_exponent_hi_20260429-152556
```

Observed trend: increasing the feedback exponent improves flat uniformity
slightly up to about `0.6-0.7`. Above that, the improvement saturates and then
begins to drift back. The current practical setting for later sweeps is
`wgs_feedback_exponent = 0.7`.

MRAF free-region factor sweep with direct WGS, `num_iters=200`,
`wgs_feedback_exponent=0.7`, `min=0.5`, and `max=1.5`:

```text
mraf_factor = 0.2, 0.4, 0.6, 0.8
mraf_factor = 1.0, 1.2
```

Artifacts:

```text
artifacts/mraf_factor_20260429-153322
artifacts/mraf_factor_overrelease_20260429-153812
```

Observed trend: larger `mraf_factor` gives the optimization more freedom in the
RTAD free region. This strongly improves flat uniformity, but it also broadens
the edge/halo structure and pushes `size50_x/y` outward. `mraf_factor=1.2`
over-releases the free region and makes the profile much too wide. The current
working conclusion is that `mraf_factor = 1.0` is a useful high-uniformity
baseline, while values above 1 are mainly diagnostic stress tests.

Release-level sweep at `mraf_factor=1.0`:

```text
release_level = 0.1353352832366127, 0.10, 0.05
release_level = 0.05, 0.025
```

Artifacts:

```text
artifacts/release_level_mraf1_20260429-155649
artifacts/release_level_mraf1_low_20260429-160039
```

Observed trend: `13.5%` and `10%` are nearly identical. Lowering the release
level to `5%` constrains more of the edge, pulls the profile size back, and
keeps most of the high uniformity from `mraf_factor=1.0`. Lowering further to
`2.5%` gives almost no additional benefit. The current practical setting is
`release_level = 0.05`.

Background attenuation sweep at `mraf_factor=1.0` and `release_level=0.05`:

```text
bg_mode = attenuate
bg_factor = 0, 0.02, 0.05, 0.1, 0.2
bg_factor = 0.2, 0.3, 0.4, 0.5
bg_factor = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
```

Artifacts are archived together under:

```text
artifacts/background_attenuation_sweeps_20260429/
  bg_factor_mraf1_release005_20260429-160705
  bg_factor_high_mraf1_release005_20260429-161055
  bg_factor_to_keep_mraf1_release005_20260429-161433
```

Observed trend: relaxing the far-background constraint has the largest positive
effect in this round. As `bg_factor` approaches `1.0`, the platform becomes
flatter and the x-direction `size50` moves close to the 330 um target. The
tradeoff is that energy is allowed to remain in the far background:
`background_fraction` rises and the diagnostics-side efficiency/side-lobe
scores worsen, especially at `bg_factor=1.0`.

Current no-background-attenuation baseline from the final background sweep
(`bg_factor_to_keep_mraf1_release005_20260429-161433/bg_factor_1p00`):

```text
method = wgs
wgs_strategy = flat_local
num_iters = 200
wgs_feedback_exponent = 0.7
wgs_weight_min = 0.5
wgs_weight_max = 1.5
mraf_factor = 1.0
release_level = 0.05
bg_mode = attenuate
bg_factor = 1.0
```

With `bg_mode="attenuate"` and `bg_factor=1.0`, the far-background projection is
effectively left unchanged. The working judgment after these sweeps is that not
attenuating the far background is the strongest improvement so far. Lowering
`release_level` to `0.05` and using `mraf_factor=1.0` also remain useful, but
the next parameter sweeps should be based on the no-background-attenuation
baseline unless the background policy itself is being retested.

After cutting off the old artifact history, a new no-background-attenuation
`mraf_factor` sweep tested `0.8`, `1.0`, and `1.2`:

```text
artifacts/mraf_factor_nobg_baseline_20260429-170901
```

Working judgment: `mraf_factor=0.8` and `1.0` are broadly similar, but `0.8`
has a better visual profile in the details. `mraf_factor=1.2` is clearly too
aggressive: it improves the nominal flatness metric, but damages the outside
shape and should not be used as the next baseline. The next sweeps use
`mraf_factor=0.8` unless this parameter is being retested.

The follow-up no-background release-level sweep used `mraf_factor=0.8` and
tested `0.05`, `0.10`, `0.135`, and `0.15`:

```text
artifacts/release_level_nobg_mraf0p8_20260429-171931
```

Working judgment: `release_level` does not strongly affect this result over the
tested range. `0.05`, `0.10`, and `0.135` are visually and metrically close, so
the next sweeps can return to the original 13.5% release setting. `0.15` is less
attractive and should not be used as the working point. Starting with the next
comparison figures, include RMS nonuniformity (`rms_nonuniformity_percent`) and
use `efficiency_e2_percent` as the standard diffraction-efficiency metric in the
on-figure summary table. Do not label `efficiency_flat` as diffraction
efficiency; it is only the power fraction inside the flat core.

The follow-up no-background WGS feedback-exponent sweep used
`mraf_factor=0.8`, `release_level=0.1353352832366127`, and tested `0.2` through
`0.8`:

```text
artifacts/wgs_feedback_exponent_nobg_mraf0p8_rel0135_20260429-172747
```

Working observation: the RMS metric improves as the exponent rises from `0.2`
to roughly `0.5`, then changes only mildly. Visually, `wgs_feedback_exponent=0.8`
has the lowest shoulder in this sweep, even though the global RMS metric is not
the absolute minimum at that point. The next sweeps use
`wgs_feedback_exponent=0.8` as the temporary visual-profile baseline.

The next no-background `wgs_weight_max` sweep used `wgs_feedback_exponent=0.8`
and tested `1.5`, `2.0`, `2.5`, and `3.0`:

```text
artifacts/wgs_weight_max_nobg_mraf0p8_rel0135_exp0p8_20260429-173633
```

Working judgment: `wgs_weight_max=1.5` through `2.5` are broadly similar in the
center profiles, with larger max improving RMS but also pushing the y-direction
size upward. `3.0` gives little extra benefit. The current working point is
`wgs_weight_max=2.0`; keep `wgs_weight_min=0.5`.

A fixed-baseline run for direct comparison against the initial phase was then
generated with the current working parameters:

```text
method = wgs
wgs_strategy = flat_local
num_iters = 200
wgs_feedback_exponent = 0.8
wgs_weight_min = 0.5
wgs_weight_max = 2.0
mraf_factor = 0.8
release_level = 0.1353352832366127
bg_mode = attenuate
bg_factor = 1.0

artifacts/fixed_baseline_initial_compare_20260429-174537
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
