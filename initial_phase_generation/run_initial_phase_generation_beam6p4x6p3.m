% run_initial_phase_generation_beam6p4x6p3
%
% Generate an initial phase using the measured elliptical Gaussian beam:
%   X 1/e^2 intensity diameter = 6.4 mm
%   Y 1/e^2 intensity diameter = 6.3 mm

clear; close all; clc;

module_root = fileparts(mfilename('fullpath'));
addpath(module_root);

cfg = default_initial_phase_config(module_root);
cfg.input_1e2_diameter_x_m = 6.4e-3;
cfg.input_1e2_diameter_y_m = 6.3e-3;
cfg.input_1e2_radius_x_m = cfg.input_1e2_diameter_x_m / 2;
cfg.input_1e2_radius_y_m = cfg.input_1e2_diameter_y_m / 2;
cfg.input_1e_radius_x_m = cfg.input_1e2_radius_x_m / sqrt(2);
cfg.input_1e_radius_y_m = cfg.input_1e2_radius_y_m / sqrt(2);
cfg.input_1e2_diameter_m = mean([cfg.input_1e2_diameter_x_m, cfg.input_1e2_diameter_y_m]);
cfg.input_1e2_radius_m = cfg.input_1e2_diameter_m / 2;
cfg.input_1e_radius_m = cfg.input_1e2_radius_m / sqrt(2);
cfg.beta_x = 2 * pi * cfg.input_1e_radius_x_m * cfg.Ro_x_m / (cfg.lambda_m * cfg.f_m);
cfg.beta_y = 2 * pi * cfg.input_1e_radius_y_m * cfg.Ro_y_m / (cfg.lambda_m * cfg.f_m);

out_dir = fullfile(cfg.output_root, 'phase0_beamX6p4mm_Y6p3mm');
phase_data = generate_initial_phase(cfg, 'do_forward', true, 'output_dir', out_dir, 'figure_dpi', cfg.figure_dpi);

fprintf('\nElliptical-beam initial phase generated.\n');
fprintf('Output folder: %s\n', out_dir);
fprintf('phase0 size: %d x %d\n', size(phase_data.phase0_unwrapped_rad, 1), size(phase_data.phase0_unwrapped_rad, 2));
fprintf('input 1/e^2 diameter x/y = %.3f / %.3f mm\n', cfg.input_1e2_diameter_x_m * 1e3, cfg.input_1e2_diameter_y_m * 1e3);
fprintf('target size = %.3f x %.3f um\n', cfg.target_size_x_m * 1e6, cfg.target_size_y_m * 1e6);
fprintf('beta_x/y = %.9f / %.9f\n', cfg.beta_x, cfg.beta_y);
