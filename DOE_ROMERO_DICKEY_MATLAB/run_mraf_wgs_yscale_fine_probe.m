% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_wgs_yscale_fine_probe()
% run_mraf_wgs_yscale_fine_probe Fine y-scale check around target_y_scale ~= 1.02.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

yscales = [1.020, 1.022];
rows = [];
for k = 1:numel(yscales)
    ys = yscales(k);
    variant_name = sprintf('wgs_exp020_yscale%04d', round(ys*1000));
    variant = struct('variant_name', string(variant_name), 'preset_name', "wgs_exp020_yscale_fine_probe", 'artifact_group', "wgs_yscale_probe", ...
        'target_edge_mode', "pivot50", 'pivot50_inner_power', 0.75, 'pivot50_outer_power', 1.25, ...
        'target_x_scale', 1.0, 'target_y_scale', ys, ...
        'free_threshold', 0.135, 'edge_low_threshold', 0.135, ...
        'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
        'phase_blend', 0.25, 'method', "WGS-MRAF", 'use_wgs', true, 'wgs_exponent', 0.20);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {k, string(variant_name), ys, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, m.rms_90, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        fullprop.mraf_metrics.size50_x_um, fullprop.mraf_metrics.size50_y_um, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

t = cell2table(rows, 'VariableNames', {'case_index','variant_name','target_y_scale','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um','shoulder_peak_x','shoulder_peak_y','core_rms','rms_90', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core','efficiency_inside_signal','efficiency_inside_guard', ...
    'fullprop_size50_x_um','fullprop_size50_y_um','fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'wgs_yscale_probe', char(stamp + "_wgs_yscale_fine_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(t, fullfile(summary_root, 'wgs_yscale_fine_metrics.csv'));
fid=fopen(fullfile(summary_root, 'wgs_yscale_fine_summary.txt'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'WGS y-scale fine probe around 1.02\n====================================\n\n');
for i=1:height(t)
    fprintf(fid,'%s target_y_scale=%.4f\n',t.variant_name(i),t.target_y_scale(i));
    fprintf(fid,'  artifact: %s\n',string(t.artifact_root(i)));
    fprintf(fid,'  size50 %.3f x %.3f um; size13.5 %.3f x %.3f um\n',t.size50_x_um(i),t.size50_y_um(i),t.size13p5_x_um(i),t.size13p5_y_um(i));
    fprintf(fid,'  TW13.5-90 %.3f x %.3f um\n',t.transition_13p5_90_x_um(i),t.transition_13p5_90_y_um(i));
    fprintf(fid,'  shoulder %.6g x %.6g; core_rms %.6g; side_lobe %.6g x %.6g; efficiency %.6g / %.6g\n\n',t.shoulder_peak_x(i),t.shoulder_peak_y(i),t.core_rms(i),t.side_lobe_peak_x_rel_to_core(i),t.side_lobe_peak_y_rel_to_core(i),t.efficiency_inside_signal(i),t.efficiency_inside_guard(i));
end
fprintf('\nWGS y-scale fine summary written to:\n%s\n', summary_root);
end

