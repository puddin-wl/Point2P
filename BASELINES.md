# Baselines

## Initial Phase Baseline

Status: verified and frozen as the current Romero-Dickey initial phase baseline.

### Phase Output

```text
E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942
```

Key file for downstream MRAF/WGS optimization:

```matlab
load('E:\program\Point2P\initial_phase_generation\artifacts\20260428-141942\phase0.mat')
```

Default parameters:

- `beta_x = 9.060975545`
- `beta_y = 3.294900198`
- `target size = 330 x 120 um`
- `lambda = 532 nm`
- `f = 429 mm`
- `input 1/e^2 diameter = 5 mm`
- `clear aperture = 15 mm`

### Diagnostic Output

```text
E:\program\Point2P\result_diagnostics\artifacts\20260428-141916
```

Key baseline metrics:

- `size50_x/y = 321.931 / 109.254 um`
- `size13p5_x/y = 370.088 / 153.315 um`
- `transition_13p5_90_x/y = 48.999 / 45.735 um`
- `core_rms = 0.163565466`
- `shoulder_peak_x/y = -0.002643 / -0.215922`
- `true_side_lobe_x/y = none / none`
- `energy_inside_signal = 0.863962`

### Policy

- Do not delete these two artifact folders.
- Do not change `initial_phase_generation` core logic unless creating a new named baseline.
- Future MRAF/WGS scripts should load this baseline `phase0.mat` rather than re-implementing initial phase generation.

