function cfg = default_initial_phase_config(module_root)
% default_initial_phase_config Standalone config for RD initial phase generation.
%
% This is a copied/minimized configuration containing only parameters needed to
% generate the Romero-Dickey phase0 and one ideal Fourier-lens focal-plane FFT.
% It intentionally does not include MRAF target/projection/sweep parameters.

if nargin < 1 || isempty(module_root)
    module_root = fileparts(mfilename('fullpath'));
end

cfg = struct();
cfg.module_root = module_root;

% Optical constants.
cfg.lambda_m = 532e-9;
cfg.f_m = 429e-3;

% Clear aperture / pupil. This is not the illuminated Gaussian diameter.
cfg.aperture_diameter_m = 15e-3;
cfg.aperture_radius_m = cfg.aperture_diameter_m / 2;

% Illuminated input Gaussian beam. Diameter is at 1/e^2 intensity.
cfg.input_1e2_diameter_m = 5e-3;
cfg.input_1e2_radius_m = cfg.input_1e2_diameter_m / 2;
cfg.input_1e_radius_m = cfg.input_1e2_radius_m / sqrt(2);

% Desired rectangular flat-top size used by the RD initial phase.
cfg.target_size_x_m = 330e-6;
cfg.target_size_y_m = 120e-6;
cfg.target_half_x_m = cfg.target_size_x_m / 2;
cfg.target_half_y_m = cfg.target_size_y_m / 2;

% Computational grid and focal-plane sampling.
cfg.N = 2048;
cfg.requested_focal_dx_m = 2.5e-6;
cfg.doe_grid_extent_m = cfg.lambda_m * cfg.f_m / cfg.requested_focal_dx_m;
cfg.dx_doe_m = cfg.doe_grid_extent_m / cfg.N;
cfg.focal_dx_m = cfg.lambda_m * cfg.f_m / (cfg.N * cfg.dx_doe_m);

% Romero-Dickey separable phase settings.
cfg.phase_method = "romero_dickey_separable";
cfg.phase_sign = 1;
cfg.phase_scale_x = 1;
cfg.phase_scale_y = 1;

% RD output scale convention and beta values.
cfg.Ro_definition = "full_size_over_sqrt_pi";
cfg.Ro_x_m = cfg.target_size_x_m / sqrt(pi);
cfg.Ro_y_m = cfg.target_size_y_m / sqrt(pi);
cfg.beta_x = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_x_m / (cfg.lambda_m * cfg.f_m);
cfg.beta_y = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_y_m / (cfg.lambda_m * cfg.f_m);

% Plot/output defaults used by run_initial_phase_generation.m.
cfg.figure_dpi = 150;
cfg.output_root = fullfile(module_root, 'artifacts');
end
