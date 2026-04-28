% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_iter_sweep()
% run_mraf_iter_sweep Small single-variable n_iter sweep for caseC_stable.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

iters = [5 10 15];
rows = [];
artifact_roots = strings(numel(iters), 1);
for k = 1:numel(iters)
    variant = struct('variant_name', "caseC_iter" + string(iters(k)), ...
        'preset_name', "caseC_stable", 'artifact_group', "caseC_iter_sweep", 'n_iter', iters(k), ...
        'rd_target_gamma', 1.0, 'mraf_factor', 1.0, 'phase_blend', 0.25, ...
        'method', "GS-MRAF", 'use_wgs', false, 'wgs_exponent', 0);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    artifact_roots(k) = string(cfg.save_root);
    m = result.mraf.metrics;
    rd = result.rd.metrics;
    fm = fullprop.mraf_metrics;
    rows = [rows; {iters(k), artifact_roots(k), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        rd.shoulder_peak_x, rd.shoulder_peak_y, rd.core_rms, ...
        fm.shoulder_peak_x, fm.shoulder_peak_y, fm.core_rms}]; %#ok<AGROW>
end

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'caseC_iter_sweep', char(stamp + "_caseC_iter_sweep_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end

sweep_table = cell2table(rows, 'VariableNames', {'n_iter','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um', ...
    'shoulder_peak_x','shoulder_peak_y','core_rms', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core', ...
    'efficiency_inside_signal','efficiency_inside_guard', ...
    'rd_shoulder_peak_x','rd_shoulder_peak_y','rd_core_rms', ...
    'fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
writetable(sweep_table, fullfile(summary_root, 'sweep_metrics.csv'));
write_sweep_summary(fullfile(summary_root, 'sweep_summary.txt'), sweep_table);
plot_sweep_summary(fullfile(summary_root, 'sweep_summary_plot.png'), sweep_table);

fprintf('\nSweep summary written to:\n%s\n', summary_root);
end

function write_sweep_summary(path, t)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'caseC_stable n_iter sweep summary\n=================================\n\n');
fprintf(fid, 'Fixed params: gamma=1.0, mraf_factor=1.0, phase_blend=0.25, method=GS-MRAF, WGS off.\n');
fprintf(fid, 'Priority: shoulder not increasing, core_rms improves, size50 near 330 x 120 um, side lobe stable.\n\n');
for k = 1:height(t)
    fprintf(fid, 'n_iter=%d\n', t.n_iter(k));
    fprintf(fid, '  artifact: %s\n', string(t.artifact_root(k)));
    fprintf(fid, '  size50 = %.3f x %.3f um, size13.5 = %.3f x %.3f um\n', t.size50_x_um(k), t.size50_y_um(k), t.size13p5_x_um(k), t.size13p5_y_um(k));
    fprintf(fid, '  TW13.5-90 = %.3f x %.3f um\n', t.transition_13p5_90_x_um(k), t.transition_13p5_90_y_um(k));
    fprintf(fid, '  shoulder = %.6g x %.6g, core_rms = %.6g\n', t.shoulder_peak_x(k), t.shoulder_peak_y(k), t.core_rms(k));
    fprintf(fid, '  side_lobe = %.6g x %.6g, efficiency = %.6g / %.6g\n\n', t.side_lobe_peak_x_rel_to_core(k), t.side_lobe_peak_y_rel_to_core(k), t.efficiency_inside_signal(k), t.efficiency_inside_guard(k));
end
[~, best_idx] = min(t.shoulder_peak_x + 0.5*t.core_rms + 0.001*abs(t.size50_x_um - 330) + 0.001*abs(t.size50_y_um - 120));
fprintf(fid, 'Suggested next baseline by conservative score: n_iter=%d\n', t.n_iter(best_idx));
end

function plot_sweep_summary(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 800]);
tiledlayout(2,2,'TileSpacing','compact');
nexttile; plot(t.n_iter, t.shoulder_peak_x, 'ro-', t.n_iter, t.shoulder_peak_y, 'bo-', 'LineWidth', 1.2); grid on; xlabel('n iter'); ylabel('shoulder peak'); legend('x','y'); title('n_iter vs shoulder');
nexttile; plot(t.n_iter, t.core_rms, 'ko-', 'LineWidth', 1.2); grid on; xlabel('n iter'); ylabel('core RMS'); title('n_iter vs core RMS');
nexttile; plot(t.n_iter, t.size50_x_um, 'ro-', t.n_iter, t.size50_y_um, 'bo-', 'LineWidth', 1.2); grid on; xlabel('n iter'); ylabel('Size50 (um)'); legend('x','y'); title('n_iter vs Size50');
nexttile; plot(t.n_iter, t.transition_13p5_90_x_um, 'ro-', t.n_iter, t.transition_13p5_90_y_um, 'bo-', 'LineWidth', 1.2); grid on; xlabel('n iter'); ylabel('TW13.5-90 (um)'); legend('x','y'); title('n_iter vs transition');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end
