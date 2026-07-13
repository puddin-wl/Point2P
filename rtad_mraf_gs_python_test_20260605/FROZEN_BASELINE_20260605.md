# Frozen Baseline 20260605

This directory is a frozen reproducible baseline copied from
`E:\program\Point2P\rtad_mraf_gs_python`.

Reference result directory:
`E:\program\Point2P\rtad_mraf_gs_python\artifacts\20260605-144020_rtad_mraf_gs_truncI0135`

Key success conditions:
- `phase0 = artifacts/phase0_beam6p5mm/phase0.mat`
- `phase_var = phase0_wrapped_rad`
- `beam = 6.5 mm`
- `swap_phase_xy = False`
- `constraint_mode = truncated_rtad`
- `release_level = exp(-2)`
- `method = wgs`
- `wgs_strategy = flat_local`
- `mraf_factor = 0.8`
- `bg_factor = 0.9`

Recommended reproduction commands:

```powershell
python run_rtad_mraf_gs_case.py
python run_pipeline.py
```

This is a frozen baseline. Future experiments should be based on a copy of
this directory. Do not directly overwrite this baseline.
