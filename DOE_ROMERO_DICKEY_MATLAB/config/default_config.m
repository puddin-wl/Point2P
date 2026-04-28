function cfg = default_config(project_root)
% default_config 创建 Romero-Dickey DOE 项目的统一参数结构。
%
% 所有物理量使用 SI 单位。默认 DOE 数值窗口由期望焦平面采样反推，
% 因此计算窗口不是 5 mm 光斑，也不是 15 mm clear aperture。

if nargin < 1 || isempty(project_root)
    project_root = fileparts(fileparts(mfilename('fullpath')));
end

cfg = struct();
cfg.project_root = project_root;

cfg.lambda_m = 532e-9;
cfg.f_m = 429e-3;
cfg.doe_to_lens_m = 200e-3;

cfg.aperture_diameter_m = 15e-3;
cfg.aperture_radius_m = cfg.aperture_diameter_m / 2;
cfg.lens_aperture_diameter_m = cfg.aperture_diameter_m;
cfg.lens_aperture_radius_m = cfg.lens_aperture_diameter_m / 2;
cfg.apply_lens_aperture_in_full_model = true;

cfg.input_1e2_diameter_m = 5e-3;
cfg.input_1e2_radius_m = cfg.input_1e2_diameter_m / 2;
cfg.input_1e_radius_m = cfg.input_1e2_radius_m / sqrt(2);

cfg.target_size_x_m = 330e-6;
cfg.target_size_y_m = 120e-6;
cfg.target_half_x_m = cfg.target_size_x_m / 2;
cfg.target_half_y_m = cfg.target_size_y_m / 2;

cfg.N = 2048;
cfg.requested_focal_dx_m = 2.5e-6;
cfg.doe_grid_extent_m = cfg.lambda_m * cfg.f_m / cfg.requested_focal_dx_m;
cfg.dx_doe_m = cfg.doe_grid_extent_m / cfg.N;
cfg.focal_dx_m = cfg.lambda_m * cfg.f_m / (cfg.N * cfg.dx_doe_m);

cfg.phase_method = "romero_dickey_separable";
cfg.phase_sign = 1;
cfg.phase_scale_x = 1;
cfg.phase_scale_y = 1;
cfg.run_opposite_phase_sign = false;

cfg.target_edge = "hard";
cfg.transition_width_x_m = 10e-6;
cfg.transition_width_y_m = 8e-6;

cfg.Ro_definition = "full_size_over_sqrt_pi";
cfg.Ro_x_m = cfg.target_size_x_m / sqrt(pi);
cfg.Ro_y_m = cfg.target_size_y_m / sqrt(pi);
cfg.beta_x = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_x_m / (cfg.lambda_m * cfg.f_m);
cfg.beta_y = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_y_m / (cfg.lambda_m * cfg.f_m);

cfg.normalization = "max";
cfg.core_fraction = 0.45;
cfg.guard_fraction = 1.5;
cfg.side_lobe_exclusion_fraction = 1.15;
cfg.profile_transition_low = 0.13;
cfg.profile_transition_high = 0.90;
cfg.center_profile_xlim_factor = 1.8;
cfg.figure_dpi = 150;

cfg.mraf = struct();
cfg.mraf.enable = true;
cfg.mraf.enabled = true;
cfg.mraf.target_mode = "flat_core_free_edge";
cfg.mraf.n_iter = 20;
cfg.mraf.method = "WGS-MRAF";
cfg.mraf.mraf_factor = 1.0;
cfg.mraf.use_three_region_projection = true;
cfg.mraf.free_factor = 1.0;
cfg.mraf.outer_factor = 0.5;
cfg.mraf.transition_cap_blend = 0.35;
cfg.mraf.phase_blend = 0.20;
cfg.mraf.use_wgs = true;
cfg.mraf.wgs_exponent = 0.10;
cfg.mraf.flat_fraction_x = 0.80;
cfg.mraf.flat_fraction_y = 0.80;
cfg.mraf.free_fraction_x = 1.35;
cfg.mraf.free_fraction_y = 1.50;
cfg.mraf.guard_fraction_x = 1.02;
cfg.mraf.guard_fraction_y = 1.04;
cfg.mraf.lobe_fraction_x = 1.25;
cfg.mraf.lobe_fraction_y = 1.35;
cfg.mraf.guard_cap_outer_intensity = 0.12;
cfg.mraf.guard_blend = 0.35;
cfg.mraf.output_dir = fullfile(project_root, "artifacts", "mraf_flat_core_free_edge");
cfg.mraf.save_every = 1;
cfg.mraf.normalization = "flat_core_mean";
cfg.mraf.preset_name = "flat_core_free_edge";
cfg.mraf.variant_name = "flat_core_caseA_baseline";
cfg.mraf.artifact_group = "flat_core_free_edge_mraf";

% LEGACY: Parameters below are kept only so older comparison scripts can run.
% The recommended workflow is run_flat_core_free_edge_mraf.m and does not use
% RD-derived gamma/pivot/sigmoid target remapping or free-threshold sweeps.
cfg.mraf.rd_target_gamma = 1.0;
cfg.mraf.rd_target_sigmoid_k = 8;
cfg.mraf.rd_target_sigmoid_t = 0.5;
cfg.mraf.rd_target_edge_mode = "gamma";
cfg.mraf.target_edge_mode = "identity";
cfg.mraf.pivot50_inner_power = 0.95;
cfg.mraf.pivot50_outer_power = 1.05;
cfg.mraf.target_x_scale = 1.0;
cfg.mraf.target_y_scale = 1.0;
cfg.mraf.core_threshold = 0.85;
cfg.mraf.signal_threshold = 0.50;
cfg.mraf.edge_low_threshold = 0.10;
cfg.mraf.noise_threshold = 0.05;
cfg.mraf.free_threshold = 0.05;
cfg.mraf.y_free_reservoir_enable = false;
cfg.mraf.y_free_reservoir_half_y_um = 72;
cfg.mraf.y_free_reservoir_max_intensity = 0.50;
cfg.mraf.wgs_y_edge_damping_enable = false;
cfg.mraf.wgs_y_edge_damping_factor = 0.50;
cfg.mraf.wgs_y_edge_damping_low = 0.50;
cfg.mraf.wgs_y_edge_damping_high = 0.90;
cfg.mraf.flat_low = 0.75;
cfg.mraf.flat_high = 0.95;


stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
cfg.save_root = fullfile(project_root, 'artifacts', char(stamp));
end
