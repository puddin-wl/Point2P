# DOE Romero-Dickey MATLAB flat-top workflow

Current recommended workflow:

```matlab
cd DOE_ROMERO_DICKEY_MATLAB
run_flat_core_free_edge_mraf
```

This project now uses a simplified MRAF refinement path:

1. Generate the initial phase with the existing Romero-Dickey / point-to-point code.
2. Use that phase as the starting phase for WGS-MRAF.
3. Constrain only the central `flat_core` to a unit-amplitude flat top.
4. Leave the surrounding `free_edge` unconstrained so transition structure and side lobes can carry error.
5. Weakly suppress the far `outer_suppress` region so energy does not spread without bound.
6. Rank results by flat-core quality first, not by minimum side lobe or narrowest transition.

## Main target: `flat_core_free_edge`

The recommended target is not RD-derived, gamma-remapped, pivot-remapped, sigmoid, or hard rectangle. It is a simple three-region target:

- `flat_core`: central region with `target.amplitude = 1`; WGS weights update only here.
- `free_edge`: surrounding transition / side-lobe reservoir; `target.amplitude = NaN`; projection keeps current amplitude times `free_factor`.
- `outer_suppress`: far outside the free region; `target.amplitude = NaN`; projection keeps current amplitude times `outer_factor`.

## Key parameters

- `flat_fraction_x/y`: fraction of the nominal 330 x 120 um target half-size that is forced flat. Larger values flatten a larger area but are harder to optimize. Smaller values flatten more easily but reduce the guaranteed flat region.
- `free_fraction_x/y`: size of the free transition / side-lobe reservoir. Larger values give side lobes more room and may improve flatness. Smaller values compress side lobes and can create shoulder or worse flatness.
- `wgs_exponent`: WGS correction strength. `0.05` is gentle, `0.10` is recommended, `0.15` is stronger and may introduce noise.
- `outer_factor`: far-field outer suppression. `1.0` means no suppression, `0.5` is recommended, `0.35` is stronger. Avoid very small values because they push energy back toward the edge.
- `n_iter`: recommended 20 to 30. Start small and inspect the trend.

## Default cases

`run_flat_core_free_edge_mraf` runs three small cases:

- Case A: recommended baseline (`wgs_exponent=0.10`, `n_iter=20`, `outer_factor=0.50`).
- Case B: stronger flattening (`wgs_exponent=0.15`, `n_iter=30`).
- Case C: slightly stronger outer suppression (`free_fraction` smaller, `outer_factor=0.35`).

Outputs are written under:

```text
artifacts/flat_core_free_edge_mraf_<timestamp>/
```

The top-level `summary.csv` contains `flat_core_rms`, `flat_core_pv`, `flat_core_uniformity`, `side_lobe_energy_ratio`, `outer_energy_ratio`, `size50_x_um`, `size50_y_um`, and transition widths.

## Legacy workflows

Older RD-derived target, gamma, pivot50, hard-rectangle, free-threshold, y-scale, and WGS probe scripts are kept only for comparison. They are not the recommended main workflow.
