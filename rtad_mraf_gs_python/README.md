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

The target is defined first as intensity. The amplitude target used by GS/MRAF
is then:

```text
A_target = sqrt(I_target)
```

The intensity size50 is fixed to 330 um x 120 um. With
`a50 = W50/2` and `b50 = H50/2`, the first version uses a separable
raised-cosine edge:

```text
a0 = a50 - delta_x_um
a1 = a50 + delta_x_um
b0 = b50 - delta_y_um
b1 = b50 + delta_y_um

C(u; u0, u1) = 1                                      for u <= u0
             = 0.5 * (1 + cos(pi*(u-u0)/(u1-u0)))     for u0 < u < u1
             = 0                                      for u >= u1

I_target(x,y) = C(abs(x); a0, a1) * C(abs(y); b0, b1)
```

Masks are separated into flat core, descending edge, support, free guard band,
and background. Flat uniformity metrics use only the flat core, not the edge.

## Refinement Formula

GS projection:

```text
E_far = FFT(A_input * exp(i phase))
E_far' = A_target * exp(i angle(E_far))
phase = angle(IFFT(E_far'))
```

MRAF projection in this implementation:

```text
signal/support:  E_far' = W * exp(i angle(E_far))
free guard band: E_far' = mraf_factor * E_far
background:      keep, attenuate, or zero according to bg_mode
```

Default `bg_mode="attenuate"` and `bg_factor=0.25` weakly constrain the far
background without hard-zeroing everything outside the support. This is
intentional for RTAD edge behavior.

WGS uses a simplified Leonardo-style update on `mask_flat` by default:

```text
W <- W * (T / F)^feedback_exponent
```

where `F` is the current farfield amplitude scaled to the target signal power.
Weights are clipped and L2 normalized after each update.

## How To Run

From this folder:

```powershell
cd E:\program\Point2P\rtad_mraf_gs_python
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_rtad_mraf_gs_case.py --phase-mat "E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat" --phase-var phase0_wrapped_rad --iters 200 --method mraf --mraf-factor 0.4
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
- `--delta-x`
- `--delta-y`
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
- `reconstruction_refined.npy`
- `metrics.csv`
- `report.txt`
- `target_intensity.png`
- `target_amplitude.png`
- `masks.png`
- `initial_reconstruction_intensity.png`
- `refined_reconstruction_intensity.png`
- `center_profiles_compare.png`
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
method = mraf
num_iters = 200
mraf_factor = 0.4 or 0.5
delta_x_um = 15
delta_y_um = 8
guard_x_um = 20
guard_y_um = 12
```

The internal metrics are quick trend checks only. Final acceptance should still
be done by the existing Point2P diagnostics once connected.
