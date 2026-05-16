% run_initial_phase_generation Standalone initial phase generation demo.
%
% Run this script from this folder or from MATLAB with this folder on path.
% It uses only files inside E:\program\Point2P\initial_phase_generation.

clear; close all; clc;

module_root = fileparts(mfilename('fullpath'));
addpath(module_root);

cfg = default_initial_phase_config(module_root);
stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
out_dir = fullfile(cfg.output_root, char(stamp));

phase_data = generate_initial_phase(cfg, 'do_forward', true, 'output_dir', out_dir, 'figure_dpi', cfg.figure_dpi);

fprintf('\nStandalone initial phase generated.\n');
fprintf('Output folder: %s\n', out_dir);
fprintf('phase0 size: %d x %d\n', size(phase_data.phase0_unwrapped_rad, 1), size(phase_data.phase0_unwrapped_rad, 2));
fprintf('lambda = %.3f nm, f = %.3f mm\n', cfg.lambda_m * 1e9, cfg.f_m * 1e3);
fprintf('input 1/e^2 diameter = %.3f mm, clear aperture = %.3f mm\n', cfg.input_1e2_diameter_m * 1e3, cfg.aperture_diameter_m * 1e3);
fprintf('target size = %.3f x %.3f um\n', cfg.target_size_x_m * 1e6, cfg.target_size_y_m * 1e6);
fprintf('beta_x/y = %.9f / %.9f\n', cfg.beta_x, cfg.beta_y);
fprintf('phase0 unwrapped min/max inside aperture = %.9g / %.9g rad\n', ...
    phase_data.phase_info.unwrapped_min_rad, phase_data.phase_info.unwrapped_max_rad);
