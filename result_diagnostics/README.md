# Result Diagnostics Module

This folder contains standalone post-processing utilities for focal-plane results. It is for analyzing already-computed intensity arrays; it does not generate phase, run MRAF, or modify targets.

## Files

- `default_diagnostics_config.m` — target size and diagnostic-level defaults.
- `compute_focal_diagnostics.m` — computes center profiles, size crossings, transition widths, core RMS/PV, shoulder, side-lobe, and energy ratios.
- `plot_center_profiles_diagnostics.m` — saves x/y center profile plots with 90%, 50%, and 13.5% crossing markers.
- `write_diagnostics_report.m` — writes a text report of key metrics.
- `run_diagnostics_example.m` — example using the standalone initial phase generator result.

## Basic Usage

```matlab
addpath('E:\program\Point2P\result_diagnostics')
cfg = default_diagnostics_config();
diagnostics = compute_focal_diagnostics(focal_x_m, focal_y_m, intensity, 'cfg', cfg);
write_diagnostics_report('diagnostics_report.txt', diagnostics);
plot_center_profiles_diagnostics(diagnostics, 'center_profiles_diagnostics.png', 'cfg', cfg);
```

## Main Metrics

- `size90_x_um`, `size90_y_um`
- `size50_x_um`, `size50_y_um`
- `size13p5_x_um`, `size13p5_y_um`
- `transition_13p5_90_x_um`, `transition_13p5_90_y_um`
- `core_rms`, `core_pv`, `core_uniformity`
- `shoulder_peak_x`, `shoulder_peak_y`
- `side_lobe_peak_x_rel_to_core`, `side_lobe_peak_y_rel_to_core`
- `true_side_lobe_peak_x_rel_to_core`, `true_side_lobe_peak_y_rel_to_core`
- `energy_inside_signal`, `energy_inside_guard`, `noise_energy_ratio`, `null_energy_ratio`

## Side-Lobe Definition

Two side-lobe style metrics are exposed:

1. `side_lobe_peak_*_rel_to_core`: max tail level outside a configurable exclusion half-width.
2. `true_side_lobe_peak_*_rel_to_core`: derivative-based local maximum after the 13.5% crossing.

The derivative-based metric is intended to distinguish a separated lobe from a monotonic shoulder/tail.

