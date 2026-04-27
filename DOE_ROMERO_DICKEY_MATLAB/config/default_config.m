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

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
cfg.save_root = fullfile(project_root, 'artifacts', char(stamp));
end
