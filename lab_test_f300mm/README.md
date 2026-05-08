# lab_test_f300mm — DOE refinement for f=300mm lab setup

Wavelength 532nm, focal length 300mm, input beam 5mm, target 330x120um.

The only parameter changed from the baseline (rtad_mraf_gs_python) is focal length: 429mm → 300mm.

## Quick start

### 1. Generate Romero-Dickey initial phase (MATLAB)

```matlab
cd E:\program\Point2P\lab_test_f300mm\matlab
run_initial_phase_generation
```

This produces `phase0.mat` under `matlab/artifacts/<timestamp>/`.

### 2. Run DOE refinement (Python, direct WGS, no MRAF warmup)

```bash
cd lab_test_f300mm
python run_rtad_mraf_gs_case.py --phase-mat matlab/artifacts/<timestamp>/phase0.mat \
    --phase-var phase0_wrapped_rad --method wgs --wgs-strategy flat_local \
    --iters 200 --wgs-feedback-exponent 0.8 --wgs-weight-min 0.5 \
    --wgs-weight-max 2.0 --bg-factor 0.9
```

### 3. (Optional) Smoke test without MATLAB

```bash
python run_rtad_mraf_gs_case.py --smoke-shape 2048 --no-cupy
```

## Key parameter impact

| Parameter | Baseline (f=429mm) | This test (f=300mm) |
|-----------|---------------------|----------------------|
| focal_length_m | 429e-3 | 300e-3 |
| dx_doe | ~44.6um | ~31.2um |
| DOE grid extent | ~91.3mm | ~63.8mm |

All other parameters (wavelength=532nm, beam=5mm, target=330x120um, N=2048, focal_dx=2.5um) are identical to baseline.

## Results (WGS 200 iterations, flat_local strategy, bg_factor=0.9)

| Metric | Value | Target |
|--------|-------|--------|
| size50_x | 329.3 μm | 330 μm |
| size50_y | 118.1 μm | 120 μm |
| size13.5_x | 343.3 μm | — |
| size13.5_y | 131.5 μm | — |
| transition_13.5_90_x | 15.3 μm | — |
| transition_13.5_90_y | 13.6 μm | — |
| RMS nonuniformity | 1.22% | < 2% |
| Efficiency (e^-2) | 94.91% | — |
| Peak overshoot | 3.5% | — |
| Derivative sidelobe x/y | ~1.8e-5 / ~1.2e-5 | — |

Both RMS nonuniformity and efficiency are better than the f=429mm baseline
(1.86% / 92.52%), as expected from the higher β values at this shorter focal length.

## Theoretical expectation: intermediate diffraction regime

The Romero-Dickey method accuracy is governed by the dimensionless parameter **β**:

```
β_x = 2π · r_{1/e} · R_{o,x} / (λ f)
β_y = 2π · r_{1/e} · R_{o,y} / (λ f)
```

where r_{1/e} = D_{1/e²} / (2√2), R_o = target_size / √π.

Larger β means the stationary-phase approximation is more accurate,
with residual error scaling as O(1/β²).

| Parameter | f = 429mm (baseline) | f = 300mm (this test) | f = 100mm | Ratio vs. 429mm |
|-----------|----------------------|-----------------------|-----------|-----------------|
| β_x | ~9.1 | 13.0 | 38.9 | 1.43× |
| β_y | ~3.3 | 4.7 | 14.1 | 1.43× |
| dx_doe | 44.6 μm | 31.2 μm | 10.4 μm | 1/1.43 |
| Fresnel number N_f | ~27 | ~39 | ~117 | 1.43× |

At f=300mm, the system sits in an **intermediate regime**:
- **x-direction**: β_x ≈ 13.0 > 10, close to the geometric-optics threshold.
  The analytical phase produces a good flat-top profile in x.
- **y-direction**: β_y ≈ 4.7 < 10, still in the diffraction-influenced regime.
  WGS refinement provides meaningful improvement over the initial phase,
  bringing RMS nonuniformity down to 1.22%.

This contrasts with f=100mm (both β_x, β_y >> 10, initial phase nearly perfect,
WGS only marginal) and f=429mm (both β_x, β_y near or below 10, WGS strongly
reshapes the profile, RMS ~1.86%).

## Output

Results go to `artifacts/<timestamp>_rtad_mraf_gs_.../`:
- `phase_refined.npy` / `phase_refined.mat` — final DOE phase
- `target.npz` — target masks and parameters
- `report.txt` — summary metrics
- `diagnostics_python/` — detailed focal-plane diagnostics
