# lab_test_f100mm — DOE refinement for f=100mm lab setup

Wavelength 532nm, focal length 100mm, input beam 5mm, target 330x120um.

The only parameter changed from the baseline (rtad_mraf_gs_python) is focal length: 429mm → 100mm.

## Quick start

### 1. Generate Romero-Dickey initial phase (MATLAB)

```matlab
cd E:\program\Point2P\lab_test_f100mm\matlab
run_initial_phase_generation
```

This produces `phase0.mat` under `matlab/artifacts/<timestamp>/`.

### 2. Run DOE refinement (Python)

```bash
cd lab_test_f100mm
python run_rtad_mraf_gs_case.py --phase-mat matlab/artifacts/<timestamp>/phase0.mat
```

### 3. (Optional) Smoke test without MATLAB

```bash
python run_rtad_mraf_gs_case.py --smoke-shape 2048 --no-cupy
```

## Key parameter impact

| Parameter | Baseline | This test |
|-----------|----------|-----------|
| focal_length_m | 429e-3 | 100e-3 |
| dx_doe | ~44.6um | ~10.4um |
| DOE grid extent | ~91.3mm | ~21.3mm |

All other parameters (wavelength=532nm, beam=5mm, target=330x120um, N=2048, focal_dx=2.5um) are identical to baseline.

## Results (WGS 200 iterations, xy_then_x strategy)

| Metric | Value |
|--------|-------|
| size50_x | 330.81 μm (target: 330) |
| size50_y | 116.94 μm (target: 120) |
| RMS nonuniformity | 0.27% |
| Efficiency (e^-2) | 99.14% |
| Peak overshoot | 1.5% |

## Why the initial phase is already nearly perfect at f=100mm

A key observation: at this short focal length, the Romero-Dickey analytical phase
is already excellent, and WGS refinement provides only marginal improvement.

### Theoretical explanation (Romero-Dickey, JOSA A 1996)

The Romero-Dickey method solves for the DOE phase via stationary-phase approximation.
Its accuracy is governed by a dimensionless parameter **β**:

```
β_x = 2π · r_{1/e} · R_{o,x} / (λ f)
β_y = 2π · r_{1/e} · R_{o,y} / (λ f)
```

where r_{1/e} = D_{1/e²} / (2√2), R_o = target_size / √π.

Larger β means the stationary-phase approximation is more accurate,
with residual error scaling as O(1/β²).

| Parameter | f = 429mm (baseline) | f = 100mm (this test) | Ratio |
|-----------|----------------------|-----------------------|-------|
| β_x | ~9.1 | 38.9 | 4.3× |
| β_y | ~3.3 | 14.1 | 4.3× |
| dx_doe | 44.6 μm | 10.4 μm | 1/4.3 |
| Fresnel number N_f | ~27 | ~117 | 4.3× |

When β > 10, the system enters the geometric-optics regime:
- Each ray from the input Gaussian maps uniquely to a focal-plane target position
- Diffraction corrections become negligible
- The analytical phase is already the exact solution of the energy-conservation equation

At f=100mm, β_x ≈ 39 and β_y ≈ 14 — deep in the geometric-optics regime.
The initial phase alone produces a near-perfect flat-top, and WGS refinement
reduces RMS nonuniformity from a few percent down to 0.27%, which is mostly
residual from the raised-cosine edge roll-off and discrete grid limits.

In contrast, at f=429mm (β_x≈9, β_y≈3), diffraction effects are still
significant, and WGS/MRAF refinement meaningfully reshapes the focal-plane profile.

### Quick start

1. Generate Romero-Dickey initial phase (MATLAB):
```matlab
cd E:\program\Point2P\lab_test_f100mm\matlab
run_initial_phase_generation
```

2. Run DOE refinement (Python, optional — initial phase alone may suffice):
```bash
cd lab_test_f100mm
python run_rtad_mraf_gs_case.py --phase-mat matlab/artifacts/<timestamp>/phase0.mat \
    --phase-var phase0_wrapped_rad --method wgs --wgs-xy-iters 80 --wgs-xonly-iters 120
```

## Output

Results go to `artifacts/<timestamp>_rtad_mraf_gs_.../`:
- `phase_refined.npy` / `phase_refined.mat` — final DOE phase
- `target.npz` — target masks and parameters
- `report.txt` — summary metrics
- `diagnostics_python/` — detailed focal-plane diagnostics
