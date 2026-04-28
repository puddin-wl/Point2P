% export_phase_for_python exports the Romero-Dickey initial phase for Python MRAF/WGS refinement.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Generate the same phase0 as run_one.m/run_mraf_one.m and save it to a simple
%   MAT file for external Python refinement or inspection.
% Inputs:
%   default_config.m, make_grid.m, gaussian_input_field.m, and
%   build_separable_phase_2d.m provide the same physical parameters and RD phase
%   path used by the MATLAB workflows.
% Outputs:
%   exports/phase_rd.mat containing phase_rd_wrapped, phase_rd_unwrapped,
%   phase_wrapped, phase_unwrapped, phase_rad, and export_meta.
% Physical meaning:
%   This exports the analytical separable Romero-Dickey phase0, not an MRAF
%   refined phase.
% Used by:
%   External Python code that wants the MATLAB RD initial phase.
% Notes:
%   NaN values outside the aperture are converted to zero for easier external
%   consumption.
clear; close all; clc;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(fileparts(project_root), 'initial_phase_generation'));

cfg = default_config(project_root);
export_dir = fullfile(project_root, 'exports');
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
end

phase_data = generate_initial_phase(cfg);
phase_rd_unwrapped = phase_data.phase0_unwrapped_rad;
phase_rd_wrapped = phase_data.phase0_wrapped_rad;
phase_info = phase_data.phase_info;
phase_rd_unwrapped(~isfinite(phase_rd_unwrapped)) = 0;
phase_rd_wrapped(~isfinite(phase_rd_wrapped)) = 0;

phase_unwrapped = phase_rd_unwrapped;
phase_wrapped = phase_rd_wrapped;
phase_rad = phase_rd_wrapped;

export_meta = struct();
export_meta.source = 'DOE_ROMERO_DICKEY_MATLAB export_phase_for_python';
export_meta.generated_at = char(datetime('now'));
export_meta.N = cfg.N;
export_meta.lambda_m = cfg.lambda_m;
export_meta.f_m = cfg.f_m;
export_meta.doe_to_lens_m = cfg.doe_to_lens_m;
export_meta.aperture_diameter_m = cfg.aperture_diameter_m;
export_meta.input_1e2_diameter_m = cfg.input_1e2_diameter_m;
export_meta.target_size_x_m = cfg.target_size_x_m;
export_meta.target_size_y_m = cfg.target_size_y_m;
export_meta.requested_focal_dx_m = cfg.requested_focal_dx_m;
export_meta.focal_dx_m = cfg.focal_dx_m;
export_meta.doe_grid_extent_m = cfg.doe_grid_extent_m;
export_meta.dx_doe_m = cfg.dx_doe_m;
export_meta.phase_method = cfg.phase_method;
export_meta.phase_sign = cfg.phase_sign;
export_meta.phase_scale_x = cfg.phase_scale_x;
export_meta.phase_scale_y = cfg.phase_scale_y;
export_meta.phase_info = phase_info;

out_path = fullfile(export_dir, 'phase_rd.mat');
save(out_path, 'phase_rd_wrapped', 'phase_rd_unwrapped', 'phase_wrapped', 'phase_unwrapped', 'phase_rad', 'export_meta', '-v7');

fprintf('\nExported Romero-Dickey phase for Python refinement:\n');
fprintf('  path: %s\n', out_path);
fprintf('  variables: phase_rd_wrapped, phase_rd_unwrapped, phase_wrapped, phase_unwrapped, phase_rad, export_meta\n');
fprintf('  shape: %d x %d\n', size(phase_rd_wrapped, 1), size(phase_rd_wrapped, 2));
fprintf('  wrapped min/max: %.9g / %.9g rad\n', min(phase_rd_wrapped(:)), max(phase_rd_wrapped(:)));
fprintf('  unwrapped min/max: %.9g / %.9g rad\n', min(phase_rd_unwrapped(:)), max(phase_rd_unwrapped(:)));
fprintf('  focal dx: %.9g um/pixel\n', cfg.focal_dx_m * 1e6);
fprintf('  target 50%% size: %.9g x %.9g um\n\n', cfg.target_size_x_m * 1e6, cfg.target_size_y_m * 1e6);
