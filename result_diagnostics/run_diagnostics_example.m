% run_diagnostics_example Demonstrate standalone post-processing diagnostics.
%
% This example generates an initial RD result via initial_phase_generation, then
% computes center profiles, sizes, transition widths, shoulder and side-lobe
% diagnostics using only result_diagnostics utilities.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
diag_root = fullfile(repo_root, 'result_diagnostics');
phase_root = fullfile(repo_root, 'initial_phase_generation');
addpath(diag_root);
addpath(phase_root);

phase_cfg = default_initial_phase_config(phase_root);
phase_data = generate_initial_phase(phase_cfg, 'do_forward', true, 'save_png', false);

diag_cfg = default_diagnostics_config();
diagnostics = compute_focal_diagnostics(phase_data.focal_x_m, phase_data.focal_y_m, phase_data.initial_intensity_norm, 'cfg', diag_cfg, 'normalization', 'max');

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
out_dir = fullfile(diag_root, 'artifacts', char(stamp));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save(fullfile(out_dir, 'diagnostics.mat'), 'diagnostics', 'diag_cfg', '-v7.3');
write_diagnostics_report(fullfile(out_dir, 'diagnostics_report.txt'), diagnostics);
plot_center_profiles_diagnostics(diagnostics, fullfile(out_dir, 'center_profiles_diagnostics.png'), 'cfg', diag_cfg);

fprintf('\nDiagnostics example written to:\n%s\n', out_dir);
fprintf('size50 x/y = %.3f / %.3f um\n', diagnostics.size50_x_um, diagnostics.size50_y_um);
fprintf('transition 13.5-90 x/y = %.3f / %.3f um\n', diagnostics.transition_13p5_90_x_um, diagnostics.transition_13p5_90_y_um);
fprintf('shoulder x/y = %.6g / %.6g\n', diagnostics.shoulder_peak_x, diagnostics.shoulder_peak_y);
