# Session Notes 2026-07-13

## Final Working Method

- Final simulation/design artifact to treat as the working production baseline:
  `E:\program\Point2P\rtad_mraf_gs_python\artifacts\20260605-144020_rtad_mraf_gs_truncI0135`
- This is the run to treat as the final method unless a later, explicitly named replacement is created.
- This artifact also has downstream SLM-related outputs (`shift_sweep`, `slm_output`, blaze/grating-related exports), so it matches the branch of work used for SLM smile/shift compensation.

### Baseline configuration in that artifact

- Initial phase:
  `E:\program\Point2P\rtad_mraf_gs_python\artifacts\phase0_beam6p5mm\phase0.mat`
- Input phase variable: `phase0_wrapped_rad`
- Beam:
  `input_gaussian_1e2_diameter_m = 6.5 mm`
- Target:
  `W50 = 330 um`, `H50 = 120 um`
- RTAD edge parameters:
  `delta_x = 15 um`, `delta_y = 8 um`
- Constraint:
  `constraint_mode = truncated_rtad`
  `release_level = exp(-2) = 0.1353352832366127`
- Refinement:
  `method = wgs`
  `wgs_strategy = flat_local`
  `num_iters = 200`
  `wgs_feedback_exponent = 0.8`
  `wgs_weight_min = 0.5`
  `wgs_weight_max = 1.5`
  `mraf_factor = 0.8`
  `bg_factor = 0.9`

## Current Experimental Issue

- The measured experimental spot size is still smaller than the design target in both axes.
- There are currently two active explanations and both must remain on the table:
  1. Optical-path distortion / alignment / scan-lens or galvo related geometric distortion.
  2. Size-definition mismatch between simulation-side metrics and experiment/vendor-side metrology.

### Important interpretation already supported by existing real-test notes

- The experiment does not appear to support one single universal width definition that reproduces the vendor-style reported size in both axes at the same time.
- Existing records in `real_test` already show that a target width of `330 x 120 um` maps to different relative profile levels on the two axes:
  - X target `330 um` corresponds to about `11.56%` of `median(flat)`
  - Y target `120 um` corresponds to about `27.46%` of `median(flat)`
- Therefore, if the experimental readout is `310 x 100 um`, the first conclusion should not be "uniform scale error"; it is more likely a mix of geometry distortion and metrology-definition mismatch.

## Practical Decision Rule

1. First adjust the optical path and reduce obvious geometric distortion/skew/shear.
2. Re-measure using the same experiment-side metrology path.
3. If the size is still systematically too small after optical alignment is judged acceptable, then compensate in the design target.

## About Enlarging the Program Target

- This is not hard.
- The current Python target generator is parameterized by `W50` and `H50`, so the simplest compensation path is to increase those target values.
- The main target parameters live in:
  `E:\program\Point2P\rtad_mraf_gs_python\config_default.py`
  and are used by:
  `E:\program\Point2P\rtad_mraf_gs_python\src\rtad_target.py`

### Important caution

- If the observed mismatch is `330 / 310` in X and `120 / 100` in Y, those are not the same scale factor:
  - X scale factor: `1.0645`
  - Y scale factor: `1.2000`
- So this is not truly isotropic enlargement.
- If compensation is needed, X and Y should usually be calibrated separately unless a later optical fix removes the anisotropy.

## Definition Caution For Future Compensation

- The RTAD paper in `text/1-s2.0-S0030399225003640-main.pdf` discusses descending edges mainly in terms of **amplitude**.
- The current Point2P Python/MATLAB target generators define the target first in **intensity**, then convert by `A = sqrt(I)`.
- Therefore, the project metric `size90` means **90% intensity width**, not necessarily the same thing as the paper's "flat-top region no less than 90% of peak amplitude".
- These two must not be mixed during target redesign.

## If the New Requirement Is "size90 should equal 330 x 120 um"

- This is straightforward only after fixing the metric convention.
- In the current Point2P codebase, the natural interpretation is:
  `size90 = 90% intensity width`
- The current RTAD target is defined by `W50/H50`, but `size90` is already determined analytically by the same raised-cosine profile.
- For the current raised-cosine intensity definition, `size90` can be converted back to the required `W50/H50` analytically when `delta_x/delta_y` stay fixed.

### Conversion with current edge parameters

- For one axis with edge half-width `delta`:
  `W90 = W50 - 1.18066894 * delta`
- Therefore:
  `W50 = W90 + 1.18066894 * delta`

With the current baseline values:

- X axis:
  `delta_x = 15 um`
  If desired `size90_x = 330 um`, use
  `W50_x = 330 + 1.18066894 * 15 = 347.71 um`
- Y axis:
  `delta_y = 8 um`
  If desired `size90_y = 120 um`, use
  `H50_y = 120 + 1.18066894 * 8 = 129.45 um`

So, keeping the same edge-shape parameters, the equivalent target would be approximately:

- `W50 = 347.7 um`
- `H50 = 129.4 um`

This is a moderate enlargement, not a difficult change.

## Recommended Next Step

- Do not change `delta_x`, `delta_y`, or `release_level` in the first compensation step.
- First keep the edge-shape definition fixed and only change the target width parameters.
- If compensation is needed after optical cleanup, test a small grid such as:
  - `W50 = 345, 348, 351 um`
  - `H50 = 128, 130, 132 um`
- Compare those against the same experiment-side readout path before changing anything more complicated.
