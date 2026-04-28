% run_initial_phase_only Generate and inspect the Romero-Dickey initial phase only.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Run only the phase0 generation path and one ideal FFT forward propagation.
%   This script does not run MRAF, does not build any MRAF target, and does not
%   modify the phase.
% Inputs:
%   default_config.m supplies wavelength, focal length, DOE grid, 15 mm clear
%   aperture, 5 mm incident Gaussian beam, and 330 x 120 um target size.
% Outputs:
%   artifacts/initial_phase_only/<timestamp>/phase0.mat, phase0.png,
%   initial_intensity.png, initial_x_profile.png, initial_y_profile.png,
%   config_snapshot.mat.
% Physical meaning:
%   phase0 is the analytical separable Romero-Dickey DOE phase used as the
%   baseline and as the initial condition for MRAF.
% Used by:
%   Manual inspection of initial phase quality without running MRAF.
% Notes:
%   This script is intentionally minimal and should not be used for target or
%   MRAF parameter tuning.

clear; close all; clc;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));
addpath(fullfile(fileparts(project_root), 'initial_phase_generation'));

cfg = default_config(project_root);
stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
out_dir = fullfile(project_root, 'artifacts', 'initial_phase_only', char(stamp));

phase_data = generate_initial_phase(cfg, 'do_forward', true, 'output_dir', out_dir, 'figure_dpi', cfg.figure_dpi);

fprintf('\nInitial phase only artifacts written to:\n%s\n', out_dir);
fprintf('phase0 size: %d x %d, unwrapped rad min/max inside aperture: %.6g / %.6g\n', ...
    size(phase_data.phase0_unwrapped_rad, 1), size(phase_data.phase0_unwrapped_rad, 2), ...
    min(phase_data.phase0_unwrapped_rad(phase_data.aperture_mask), [], 'all'), max(phase_data.phase0_unwrapped_rad(phase_data.aperture_mask), [], 'all'));
fprintf('focal dx: %.6f um/pixel; target size: %.3f x %.3f um\n', cfg.focal_dx_m*1e6, cfg.target_size_x_m*1e6, cfg.target_size_y_m*1e6);
