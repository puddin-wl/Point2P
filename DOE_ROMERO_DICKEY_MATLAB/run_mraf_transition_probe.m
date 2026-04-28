% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_transition_probe()
% run_mraf_transition_probe Conservative gamma-only probe around caseC_iter10_baseline.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

gammas = [1.00 1.02 1.04 1.06 1.08];
rows = [];
for k = 1:numel(gammas)
    gamma_value = gammas(k);
    gamma_tag = sprintf('gamma%03d', round(gamma_value * 100));
    variant = struct('variant_name', "transition_probe_" + string(gamma_tag), ...
        'preset_name', "caseC_stable", 'artifact_group', "transition_gamma_probe", 'n_iter', 10, ...
        'rd_target_gamma', gamma_value, 'mraf_factor', 1.0, 'phase_blend', 0.25, ...
        'method', "GS-MRAF", 'use_wgs', false, 'wgs_exponent', 0);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rd = result.rd.metrics;
    rows = [rows; {gamma_value, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        rd.size50_x_um, rd.size50_y_um, rd.size13p5_x_um, rd.size13p5_y_um, ...
        rd.transition_13p5_90_x_um, rd.transition_13p5_90_y_um, ...
        rd.shoulder_peak_x, rd.shoulder_peak_y, rd.core_rms, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

probe_table = cell2table(rows, 'VariableNames', {'gamma','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um', ...
    'shoulder_peak_x','shoulder_peak_y','core_rms', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core', ...
    'efficiency_inside_signal','efficiency_inside_guard', ...
    'rd_size50_x_um','rd_size50_y_um','rd_size13p5_x_um','rd_size13p5_y_um', ...
    'rd_transition_13p5_90_x_um','rd_transition_13p5_90_y_um', ...
    'rd_shoulder_peak_x','rd_shoulder_peak_y','rd_core_rms', ...
    'fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
probe_table = add_pass_fail(probe_table);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'transition_gamma_probe', char(stamp + "_transition_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(probe_table, fullfile(summary_root, 'transition_probe_metrics.csv'));
write_transition_probe_summary(fullfile(summary_root, 'transition_probe_summary.txt'), probe_table);
plot_transition_probe_summary(fullfile(summary_root, 'transition_probe_summary_plot.png'), probe_table);
plot_transition_probe_profile_overlay(fullfile(summary_root, 'transition_probe_profile_overlay.png'), probe_table);

fprintf('\nTransition probe summary written to:\n%s\n', summary_root);
end

function t = add_pass_fail(t)
b = t(1, :);
n = height(t);
pass = false(n, 1);
reason = strings(n, 1);
for i = 1:n
    checks = true(8, 1);
    messages = strings(8, 1);
    checks(1) = t.shoulder_peak_x(i) <= b.shoulder_peak_x + 0.01;
    messages(1) = "shoulder_x";
    checks(2) = t.shoulder_peak_y(i) <= b.shoulder_peak_y + 0.01;
    messages(2) = "shoulder_y";
    checks(3) = t.core_rms(i) <= b.core_rms * 1.10;
    messages(3) = "core_rms";
    checks(4) = abs(t.size50_x_um(i) - 330) <= 5;
    messages(4) = "size50_x";
    checks(5) = abs(t.size50_y_um(i) - 120) <= 5;
    messages(5) = "size50_y";
    checks(6) = t.side_lobe_peak_x_rel_to_core(i) <= b.side_lobe_peak_x_rel_to_core + 0.02 && ...
        t.side_lobe_peak_y_rel_to_core(i) <= b.side_lobe_peak_y_rel_to_core + 0.02;
    messages(6) = "side_lobe";
    checks(7) = t.efficiency_inside_signal(i) >= b.efficiency_inside_signal * 0.98;
    messages(7) = "efficiency";
    checks(8) = isfinite(t.transition_13p5_90_x_um(i)) && isfinite(t.transition_13p5_90_y_um(i));
    messages(8) = "finite_transition";
    pass(i) = all(checks);
    failed = messages(~checks);
    if isempty(failed)
        reason(i) = "PASS";
    else
        reason(i) = "FAIL: " + strjoin(failed, ", ");
    end
end
t.pass = pass;
t.pass_fail_reason = reason;
t.delta_transition_x_um = t.transition_13p5_90_x_um - b.transition_13p5_90_x_um;
t.delta_transition_y_um = t.transition_13p5_90_y_um - b.transition_13p5_90_y_um;
t.delta_shoulder_x = t.shoulder_peak_x - b.shoulder_peak_x;
t.delta_shoulder_y = t.shoulder_peak_y - b.shoulder_peak_y;
t.delta_core_rms = t.core_rms - b.core_rms;
end

function write_transition_probe_summary(path, t)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Transition probe summary\n========================\n\n');
fprintf(fid, 'Baseline: gamma=1.00, n_iter=10, mraf_factor=1.0, phase_blend=0.25, GS-MRAF, WGS off.\n');
fprintf(fid, 'Only rd_target_gamma is varied: 1.00, 1.02, 1.04, 1.06, 1.08.\n\n');
fprintf(fid, 'Pass criteria relative to gamma=1.00 baseline:\n');
fprintf(fid, '- shoulder_peak_x/y increase <= 0.01\n');
fprintf(fid, '- core_rms increase <= 10%%\n');
fprintf(fid, '- size50_x within 330 +/- 5 um, size50_y within 120 +/- 5 um\n');
fprintf(fid, '- side_lobe x/y increase <= 0.02\n');
fprintf(fid, '- efficiency_inside_signal drop <= 2%%\n\n');
for i = 1:height(t)
    fprintf(fid, 'gamma = %.2f  %s\n', t.gamma(i), t.pass_fail_reason(i));
    fprintf(fid, '  artifact: %s\n', string(t.artifact_root(i)));
    fprintf(fid, '  size50 = %.3f x %.3f um; size13.5 = %.3f x %.3f um\n', t.size50_x_um(i), t.size50_y_um(i), t.size13p5_x_um(i), t.size13p5_y_um(i));
    fprintf(fid, '  TW13.5-90 = %.3f x %.3f um; delta = %.3f x %.3f um\n', t.transition_13p5_90_x_um(i), t.transition_13p5_90_y_um(i), t.delta_transition_x_um(i), t.delta_transition_y_um(i));
    fprintf(fid, '  shoulder = %.6g x %.6g; delta = %.6g x %.6g\n', t.shoulder_peak_x(i), t.shoulder_peak_y(i), t.delta_shoulder_x(i), t.delta_shoulder_y(i));
    fprintf(fid, '  core_rms = %.6g; side_lobe = %.6g x %.6g; efficiency_signal = %.6g\n\n', t.core_rms(i), t.side_lobe_peak_x_rel_to_core(i), t.side_lobe_peak_y_rel_to_core(i), t.efficiency_inside_signal(i));
end
acceptable = t(t.pass, :);
if isempty(acceptable)
    fprintf(fid, 'No gamma probe passed. Do not continue increasing gamma. Suggested next directions only: phase_blend micro-scan, flat_low micro-scan, or separately designed outer-edge-only sharpen.\n');
else
    score = acceptable.delta_transition_x_um + acceptable.delta_transition_y_um + 20*max(acceptable.delta_shoulder_x, 0) + 20*max(acceptable.delta_shoulder_y, 0) + 5*max(acceptable.delta_core_rms, 0);
    [~, idx] = min(score);
    fprintf(fid, 'Best conservative candidate among PASS rows: gamma=%.2f\n', acceptable.gamma(idx));
end
end

function plot_transition_probe_summary(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1350 900]);
tiledlayout(3,2,'TileSpacing','compact');
base_gamma = 1.00;
nexttile; plot_pair(t.gamma, t.transition_13p5_90_x_um, t.transition_13p5_90_y_um, 'TW13.5-90 (um)', 'gamma vs transition'); xline(base_gamma, 'k--');
nexttile; plot_pair(t.gamma, t.shoulder_peak_x, t.shoulder_peak_y, 'shoulder peak', 'gamma vs shoulder'); xline(base_gamma, 'k--');
nexttile; plot(t.gamma, t.core_rms, 'ko-', 'LineWidth', 1.2); grid on; xlabel('gamma'); ylabel('core RMS'); title('gamma vs core RMS'); xline(base_gamma, 'k--');
nexttile; plot_pair(t.gamma, t.size50_x_um, t.size50_y_um, 'Size50 (um)', 'gamma vs Size50'); xline(base_gamma, 'k--');
nexttile; plot_pair(t.gamma, t.side_lobe_peak_x_rel_to_core, t.side_lobe_peak_y_rel_to_core, 'side lobe rel core', 'gamma vs side lobe'); xline(base_gamma, 'k--');
nexttile; plot(t.gamma, t.efficiency_inside_signal, 'ko-', 'LineWidth', 1.2); grid on; xlabel('gamma'); ylabel('efficiency inside signal'); title('gamma vs efficiency'); xline(base_gamma, 'k--');
sgtitle('Transition probe: gamma-only, caseC_iter10_baseline fixed');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end

function plot_pair(x, yx, yy, ylabel_text, title_text)
plot(x, yx, 'ro-', x, yy, 'bo-', 'LineWidth', 1.2); grid on;
xlabel('gamma'); ylabel(ylabel_text); title(title_text); legend('x','y', 'Location', 'best');
end

function plot_transition_probe_profile_overlay(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 800]);
tiledlayout(1,2,'TileSpacing','compact');
colors = lines(height(t));
for i = 1:height(t)
    x_csv = fullfile(string(t.artifact_root(i)), 'ideal_fft', 'x_profile_rd_target_mraf.csv');
    y_csv = fullfile(string(t.artifact_root(i)), 'ideal_fft', 'y_profile_rd_target_mraf.csv');
    xt = readtable(x_csv);
    yt = readtable(y_csv);
    ax = nexttile(1);
    if i == 1
        plot(ax, xt.x_m*1e6, xt.rd_intensity_core_norm, 'k--', 'LineWidth', 1.3, 'DisplayName', 'RD baseline'); hold(ax, 'on');
    end
    if abs(t.gamma(i) - 1.0) < 1e-12
        name = 'caseC baseline gamma=1.00'; width = 1.8;
    else
        name = sprintf('gamma=%.2f', t.gamma(i)); width = 1.0;
    end
    plot(ax, xt.x_m*1e6, xt.mraf_intensity_core_norm, '-', 'Color', colors(i,:), 'LineWidth', width, 'DisplayName', name);
    ax = nexttile(2);
    if i == 1
        plot(ax, yt.y_m*1e6, yt.rd_intensity_core_norm, 'k--', 'LineWidth', 1.3, 'DisplayName', 'RD baseline'); hold(ax, 'on');
    end
    plot(ax, yt.y_m*1e6, yt.mraf_intensity_core_norm, '-', 'Color', colors(i,:), 'LineWidth', width, 'DisplayName', name);
end
ax = nexttile(1); grid(ax, 'on'); xlim(ax, [-360 360]); ylim(ax, [0 1.35]); yline(ax, 0.90, 'g:'); yline(ax, 0.50, 'r:'); yline(ax, 0.135, 'm:'); title(ax, 'X center profiles: RD, caseC baseline, gamma probes'); xlabel(ax, 'x (um)'); ylabel(ax, 'normalized intensity'); legend(ax, 'Location', 'northeastoutside');
ax = nexttile(2); grid(ax, 'on'); xlim(ax, [-220 220]); ylim(ax, [0 1.35]); yline(ax, 0.90, 'g:'); yline(ax, 0.50, 'r:'); yline(ax, 0.135, 'm:'); title(ax, 'Y center profiles: RD, caseC baseline, gamma probes'); xlabel(ax, 'y (um)'); ylabel(ax, 'normalized intensity'); legend(ax, 'Location', 'northeastoutside');
sgtitle('Transition probe profile overlay: check shoulder before accepting narrower transition');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end

