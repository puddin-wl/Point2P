% run_one 运行 Romero-Dickey 矩形平顶 DOE 相位设计与诊断。
%
% 本脚本生成默认配置、DOE 面输入场、解析 separable phase-only 相位，
% 使用 Fourier lens 焦平面 FFT 模型传播，并保存诊断图、剖面和指标。

clear; close all; clc;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(fileparts(project_root), 'initial_phase_generation'));

cfg = default_config(project_root);
if ~exist(cfg.save_root, 'dir')
    mkdir(cfg.save_root);
end

fprintf('\n=== Romero-Dickey MATLAB DOE design ===\n');
fprintf('DOE grid size: N = %d, extent = %.6f mm\n', cfg.N, cfg.doe_grid_extent_m * 1e3);
fprintf('DOE dx: %.6f um\n', cfg.dx_doe_m * 1e6);
fprintf('Focal plane dx: %.6f um/pixel\n', cfg.focal_dx_m * 1e6);
fprintf('Clear aperture: %.6f mm diameter\n', cfg.aperture_diameter_m * 1e3);
fprintf('Input 1/e^2 intensity diameter: %.6f mm\n', cfg.input_1e2_diameter_m * 1e3);
fprintf('Converted input 1/e intensity radius ri: %.6f mm\n', cfg.input_1e_radius_m * 1e3);
fprintf('Target size: %.6f um x %.6f um at 50%% intensity boundary\n', cfg.target_size_x_m * 1e6, cfg.target_size_y_m * 1e6);
fprintf('Ro convention: %s\n', cfg.Ro_definition);
fprintf('beta_x = %.6f, beta_y = %.6f\n', cfg.beta_x, cfg.beta_y);
fprintf('Beta note: smaller beta means a harder stationary-phase / finite-aperture design; beta_y is the limiting direction here.\n');
fprintf('Save root: %s\n\n', cfg.save_root);

phase_data = generate_initial_phase(cfg);
X_m = phase_data.X_m;
Y_m = phase_data.Y_m;
x_m = phase_data.x_m;
y_m = phase_data.y_m;
grid = phase_data.grid;
input_field = phase_data.input_field;
input_amplitude = phase_data.input_amplitude;
input_intensity = phase_data.input_intensity;
aperture_mask = phase_data.aperture_mask;
phase_unwrapped_rad = phase_data.phase0_unwrapped_rad;
phase_wrapped_rad = phase_data.phase0_wrapped_rad;
phase_info = phase_data.phase_info;
doe_field = input_field .* exp(1i * phase_unwrapped_rad);

fprintf('Step 1/2: ideal Fourier-lens FFT baseline...\n');
[ideal_focal_x_m, ideal_focal_y_m, ideal_intensity, ideal_focal] = fresnel_to_focal_fft(doe_field, cfg.dx_doe_m, cfg.lambda_m, cfg.f_m);
[ideal_target_intensity, ideal_target_info] = rectangular_target(ideal_focal_x_m, ideal_focal_y_m, cfg);
ideal_metrics = compute_metrics(ideal_focal_x_m, ideal_focal_y_m, ideal_intensity, ideal_target_intensity, cfg);

fprintf('Step 2/2: full propagation model with %.3f mm DOE-to-lens distance...\n', cfg.doe_to_lens_m * 1e3);
[full_focal_x_m, full_focal_y_m, full_intensity, full_focal] = full_lens_propagate(doe_field, X_m, Y_m, cfg.dx_doe_m, cfg);
[full_target_intensity, full_target_info] = rectangular_target(full_focal_x_m, full_focal_y_m, cfg);
full_metrics = compute_metrics(full_focal_x_m, full_focal_y_m, full_intensity, full_target_intensity, cfg);

save_artifacts(cfg, grid, phase_info, ideal_focal, ideal_target_info, ideal_metrics, ...
    x_m, y_m, ideal_focal_x_m, ideal_focal_y_m, input_amplitude, input_intensity, ...
    aperture_mask, phase_unwrapped_rad, phase_wrapped_rad, ideal_intensity, ideal_target_intensity, 'ideal_');
save_artifacts(cfg, grid, phase_info, full_focal, full_target_info, full_metrics, ...
    x_m, y_m, full_focal_x_m, full_focal_y_m, input_amplitude, input_intensity, ...
    aperture_mask, phase_unwrapped_rad, phase_wrapped_rad, full_intensity, full_target_intensity, 'full_');
save_comparison_report(cfg, ideal_metrics, full_metrics);

plot_diagnostics(cfg, ideal_focal_x_m, ideal_focal_y_m, ideal_intensity, ideal_target_intensity, ...
    input_amplitude, input_intensity, aperture_mask, phase_wrapped_rad, ideal_metrics, 'ideal_');
plot_diagnostics(cfg, full_focal_x_m, full_focal_y_m, full_intensity, full_target_intensity, ...
    input_amplitude, input_intensity, aperture_mask, phase_wrapped_rad, full_metrics, 'full_');

fprintf('\nDone. Artifacts written to:\n%s\n', cfg.save_root);
