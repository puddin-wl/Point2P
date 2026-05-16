function cfg = experiments_wgs_config()
% Experiment WGS configuration parameters.
% All physical parameters must match lab_test_f200mm/.

%% Physical parameters
cfg.lambda_m        = 532e-9;
cfg.f_m             = 200e-3;
cfg.aperture_diameter_m = 15e-3;
cfg.input_1e2_diameter_m = 7.0e-3;  % matches best f200mm result (RMS 0.18%)

%% Computation grid (matches simulation)
cfg.N               = 2048;
cfg.focal_dx_um     = 2.5;
cfg.focal_dy_um     = 2.5;

%% Derived grid parameters
cfg.doe_extent_m    = cfg.lambda_m * cfg.f_m / (cfg.focal_dx_um * 1e-6);
cfg.dx_doe_m        = cfg.doe_extent_m / cfg.N;
cfg.focal_dx_m      = cfg.lambda_m * cfg.f_m / (cfg.N * cfg.dx_doe_m);

%% Target geometry
cfg.W50_um          = 330;
cfg.H50_um          = 120;
cfg.delta_x_um      = 15;
cfg.delta_y_um      = 8;
cfg.guard_x_um      = 20;
cfg.guard_y_um      = 12;
cfg.release_level   = exp(-2);
cfg.constraint_mode = 'truncated_rtad';
cfg.target_mode     = 'separable';

%% SLM
cfg.slm_width       = 1920;
cfg.slm_height      = 1080;
cfg.slm_pitch_um    = 6.4;

%% Camera
cfg.camera_pixel_um = 3.45;

%% WGS parameters (more conservative than simulation)
cfg.max_iters       = 30;
cfg.wgs_feedback_exponent = 0.5;
cfg.wgs_weight_min  = 0.5;
cfg.wgs_weight_max  = 1.5;
cfg.wgs_update_every = 5;
cfg.bg_factor       = 0.9;
cfg.rms_target      = 0.02;
cfg.convergence_window = 5;

%% Paths
cfg.phase0_path     = 'phase_refined.mat';  % copy from lab_test_f200mm artifacts
cfg.output_root     = 'artifacts';
