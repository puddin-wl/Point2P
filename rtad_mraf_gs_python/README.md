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

Recommended MRAF followed by flat-core WGS:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --phase-mat "E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat" --phase-var phase0_wrapped_rad --method mraf_then_wgs --mraf-iters 150 --wgs-iters 50 --mraf-factor 0.4 --wgs-feedback-exponent 0.3 --wgs-update-every 5 --release-level 0.1353352832366127 --bg-factor 0.05
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
- `reconstruction_refined.npy`
- `reconstruction_after_mraf.npy`
- `wgs_weights_final.npy`
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
- `wgs_weights_final.png`
- `wgs_weight_hist.png`
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

## First Recommended Parameters

```text
method = mraf_then_wgs
mraf_iters = 150
wgs_iters = 50
mraf_factor = 0.4
constraint_mode = truncated_rtad
release_level = 0.1353352832366127
bg_factor = 0.05
wgs_update_mask = flat
wgs_feedback = amplitude
wgs_feedback_exponent = 0.3
wgs_update_every = 5
wgs_weight_min = 0.5
wgs_weight_max = 2.0
delta_x_um = 15
delta_y_um = 8
guard_x_um = 20
guard_y_um = 12
```

The internal metrics are quick trend checks only. Final acceptance should still
be done by the existing Point2P diagnostics once connected.
