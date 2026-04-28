# Initial Phase Generation Module

This folder contains the standalone Romero-Dickey initial phase generator used by the MATLAB DOE project. It is self-contained: the files in this folder are enough to generate the initial phase and one ideal Fourier-lens focal-plane flat-top result.

## Files

- `default_initial_phase_config.m` — standalone copy of all configuration values needed for initial phase generation.
- `generate_initial_phase.m` — standalone generator containing the useful grid, Gaussian input, RD phase, wrapping, FFT, and save logic.
- `run_initial_phase_generation.m` — one-command demo script that saves phase and focal-plane plots.

## Main Entry

```matlab
phase_data = generate_initial_phase(cfg);
```

If no config is passed, the module uses `default_initial_phase_config.m`:

```matlab
phase_data = generate_initial_phase([]);
```

## One-Command Run

```matlab
cd E:\program\Point2P\initial_phase_generation
run_initial_phase_generation
```

Outputs are written to:

```text
E:\program\Point2P\initial_phase_generation\artifacts\<timestamp>\
```

Optional forward-propagation inspection and saving:

```matlab
phase_data = generate_initial_phase(cfg, ...
    'do_forward', true, ...
    'output_dir', fullfile(project_root, 'artifacts', 'initial_phase_only', 'test'));
```

## What It Generates

- `phase_data.phase0_unwrapped_rad`: unwrapped analytical Romero-Dickey phase in radians.
- `phase_data.phase0_wrapped_rad`: wrapped phase in `[0, 2*pi)`, with `NaN` outside aperture.
- `phase_data.input_amplitude`: DOE-plane Gaussian amplitude.
- `phase_data.aperture_mask`: clear aperture mask.
- `phase_data.focal_x_m/focal_y_m`: focal-plane coordinates.
- `phase_data.initial_intensity_norm`: optional ideal FFT intensity when `'do_forward'` is true.

## Design Rule

Downstream scripts should call `generate_initial_phase(cfg)` instead of manually calling:

```matlab
make_grid
gaussian_input_field
build_separable_phase_2d
```

The legacy functions remain in `DOE_ROMERO_DICKEY_MATLAB/src` for compatibility and documentation, but the clean source of phase0 is this folder.

## Current Physics

The generated phase is a separable Romero-Dickey analytical phase:

```text
phase0(x,y) = phase_x_RD(x) + phase_y_RD(y)
```

It uses the project configuration values for wavelength, focal length, 5 mm incident Gaussian beam, 15 mm clear aperture, and 330 x 120 um target size.
