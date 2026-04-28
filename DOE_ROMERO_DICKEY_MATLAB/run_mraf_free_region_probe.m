% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function summary_root = run_mraf_free_region_probe()
% run_mraf_free_region_probe Probe free/noise threshold around caseC baseline.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

thresholds = [0.05 0.08 0.10 0.135 0.15];
rows = [];
for k = 1:numel(thresholds)
    threshold = thresholds(k);
    tag = threshold_tag(threshold);
    variant = struct('variant_name', "free_probe_" + string(tag), ...
        'preset_name', "caseC_stable", 'artifact_group', "free_region_probe", 'n_iter', 10, ...
        'rd_target_gamma', 1.0, 'mraf_factor', 1.0, 'phase_blend', 0.25, ...
        'method', "GS-MRAF", 'use_wgs', false, 'wgs_exponent', 0, ...
        'free_threshold', threshold, 'edge_low_threshold', threshold);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    m = result.mraf.metrics;
    rows = [rows; {threshold, string(cfg.save_root), ...
        m.size50_x_um, m.size50_y_um, m.size13p5_x_um, m.size13p5_y_um, ...
        m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, ...
        m.shoulder_peak_x, m.shoulder_peak_y, m.core_rms, ...
        m.side_lobe_peak_x_rel_to_core, m.side_lobe_peak_y_rel_to_core, ...
        m.efficiency_inside_signal, m.efficiency_inside_guard, ...
        result.rd.metrics.shoulder_peak_x, result.rd.metrics.shoulder_peak_y, result.rd.metrics.core_rms, ...
        fullprop.mraf_metrics.shoulder_peak_x, fullprop.mraf_metrics.shoulder_peak_y, fullprop.mraf_metrics.core_rms}]; %#ok<AGROW>
end

probe_table = cell2table(rows, 'VariableNames', {'free_threshold','artifact_root', ...
    'size50_x_um','size50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um', ...
    'shoulder_peak_x','shoulder_peak_y','core_rms', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core', ...
    'efficiency_inside_signal','efficiency_inside_guard', ...
    'rd_shoulder_peak_x','rd_shoulder_peak_y','rd_core_rms', ...
    'fullprop_shoulder_peak_x','fullprop_shoulder_peak_y','fullprop_core_rms'});
probe_table = add_free_probe_status(probe_table);

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
summary_root = fullfile(project_root, 'artifacts', 'free_region_probe', char(stamp + "_free_probe_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(probe_table, fullfile(summary_root, 'free_probe_metrics.csv'));
write_free_probe_summary(fullfile(summary_root, 'free_probe_summary.txt'), probe_table);
plot_free_probe_summary(fullfile(summary_root, 'free_probe_summary_plot.png'), probe_table);
plot_free_probe_profile_overlay(fullfile(summary_root, 'free_probe_profile_overlay.png'), probe_table);

fprintf('\nFree region probe summary written to:\n%s\n', summary_root);
end

function tag = threshold_tag(threshold)
if abs(threshold - 0.135) < 1e-12
    tag = 'thr0135';
else
    tag = sprintf('thr%03d', round(threshold * 100));
end
end

function t = add_free_probe_status(t)
b = t(1, :);
n = height(t);
status = strings(n, 1);
reason = strings(n, 1);
for i = 1:n
    hard_fail = strings(0, 1);
    warning = strings(0, 1);
    if t.shoulder_peak_x(i) > b.shoulder_peak_x + 0.005, hard_fail(end+1) = "shoulder_x"; end %#ok<AGROW>
    if t.shoulder_peak_y(i) > b.shoulder_peak_y + 0.005, hard_fail(end+1) = "shoulder_y"; end %#ok<AGROW>
    if t.core_rms(i) > b.core_rms * 1.10, hard_fail(end+1) = "core_rms"; end %#ok<AGROW>
    if abs(t.size50_x_um(i) - 330) > 5, hard_fail(end+1) = "size50_x"; end %#ok<AGROW>
    if abs(t.size50_y_um(i) - 120) > 5, hard_fail(end+1) = "size50_y"; end %#ok<AGROW>
    if t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core + 0.05 || t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core + 0.05
        hard_fail(end+1) = "side_lobe_large"; %#ok<AGROW>
    elseif t.side_lobe_peak_x_rel_to_core(i) > b.side_lobe_peak_x_rel_to_core || t.side_lobe_peak_y_rel_to_core(i) > b.side_lobe_peak_y_rel_to_core
        warning(end+1) = "side_lobe"; %#ok<AGROW>
    end
    if t.efficiency_inside_signal(i) < b.efficiency_inside_signal
        warning(end+1) = "efficiency_signal"; %#ok<AGROW>
    end
    if t.size13p5_x_um(i) > b.size13p5_x_um + 3 || t.size13p5_y_um(i) > b.size13p5_y_um + 3
        warning(end+1) = "size13p5_expansion"; %#ok<AGROW>
    end
    if ~isempty(hard_fail)
        status(i) = "FAIL";
        reason(i) = "FAIL: " + strjoin(hard_fail, ", ");
    elseif ~isempty(warning)
        status(i) = "PASS_WITH_SIDELOBE_WARNING";
        reason(i) = "WARN: " + strjoin(warning, ", ");
    else
        status(i) = "PASS";
        reason(i) = "PASS";
    end
end
t.status = status;
t.pass_fail_reason = reason;
t.delta_shoulder_x = t.shoulder_peak_x - b.shoulder_peak_x;
t.delta_shoulder_y = t.shoulder_peak_y - b.shoulder_peak_y;
t.delta_transition_x_um = t.transition_13p5_90_x_um - b.transition_13p5_90_x_um;
t.delta_transition_y_um = t.transition_13p5_90_y_um - b.transition_13p5_90_y_um;
t.delta_core_rms = t.core_rms - b.core_rms;
t.delta_side_lobe_x = t.side_lobe_peak_x_rel_to_core - b.side_lobe_peak_x_rel_to_core;
t.delta_side_lobe_y = t.side_lobe_peak_y_rel_to_core - b.side_lobe_peak_y_rel_to_core;
end

function write_free_probe_summary(path, t)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Free/noise region probe summary\n===============================\n\n');
fprintf(fid, 'Fixed params: n_iter=10, gamma=1.0, mraf_factor=1.0, phase_blend=0.25, GS-MRAF, WGS off.\n');
fprintf(fid, 'Only free_threshold is varied. free/noise = I_rd_norm < free_threshold or outside guard.\n');
fprintf(fid, 'Baseline free_threshold = 0.05. Side lobe is warning unless +0.05 over baseline.\n\n');
for i = 1:height(t)
    fprintf(fid, 'free_threshold = %.3f  %s\n', t.free_threshold(i), t.pass_fail_reason(i));
    fprintf(fid, '  artifact: %s\n', string(t.artifact_root(i)));
    fprintf(fid, '  size50 = %.3f x %.3f um; size13.5 = %.3f x %.3f um\n', t.size50_x_um(i), t.size50_y_um(i), t.size13p5_x_um(i), t.size13p5_y_um(i));
    fprintf(fid, '  TW13.5-90 = %.3f x %.3f um; delta = %.3f x %.3f um\n', t.transition_13p5_90_x_um(i), t.transition_13p5_90_y_um(i), t.delta_transition_x_um(i), t.delta_transition_y_um(i));
    fprintf(fid, '  shoulder = %.6g x %.6g; delta = %.6g x %.6g\n', t.shoulder_peak_x(i), t.shoulder_peak_y(i), t.delta_shoulder_x(i), t.delta_shoulder_y(i));
    fprintf(fid, '  core_rms = %.6g; side_lobe = %.6g x %.6g; delta side_lobe = %.6g x %.6g\n', t.core_rms(i), t.side_lobe_peak_x_rel_to_core(i), t.side_lobe_peak_y_rel_to_core(i), t.delta_side_lobe_x(i), t.delta_side_lobe_y(i));
    fprintf(fid, '  efficiency signal/guard = %.6g / %.6g\n\n', t.efficiency_inside_signal(i), t.efficiency_inside_guard(i));
end
ok = t(t.status ~= "FAIL", :);
if isempty(ok)
    fprintf(fid, 'No free_threshold candidate passed hard constraints. Keep baseline free_threshold=0.05.\n');
else
    score = ok.delta_shoulder_x + ok.delta_shoulder_y + 0.2*ok.core_rms + 0.01*abs(ok.delta_transition_x_um) + 0.01*abs(ok.delta_transition_y_um) + 0.2*max(ok.delta_side_lobe_x, 0) + 0.2*max(ok.delta_side_lobe_y, 0);
    [~, idx] = min(score);
    fprintf(fid, 'Suggested conservative candidate: free_threshold=%.3f (%s)\n', ok.free_threshold(idx), ok.status(idx));
end
end

function plot_free_probe_summary(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 1050]);
tiledlayout(4,2,'TileSpacing','compact');
base = 0.05;
nexttile; plot_pair(t.free_threshold, t.shoulder_peak_x, t.shoulder_peak_y, 'shoulder peak', 'free threshold vs shoulder', base);
nexttile; plot_pair(t.free_threshold, t.transition_13p5_90_x_um, t.transition_13p5_90_y_um, 'TW13.5-90 (um)', 'free threshold vs transition', base);
nexttile; plot(t.free_threshold, t.core_rms, 'ko-', 'LineWidth', 1.2); grid on; xlabel('free threshold'); ylabel('core RMS'); title('free threshold vs core RMS'); xline(base, 'k--');
nexttile; plot_pair(t.free_threshold, t.size50_x_um, t.size50_y_um, 'Size50 (um)', 'free threshold vs Size50', base);
nexttile; plot_pair(t.free_threshold, t.size13p5_x_um, t.size13p5_y_um, 'Size13.5 (um)', 'free threshold vs Size13.5', base);
nexttile; plot_pair(t.free_threshold, t.side_lobe_peak_x_rel_to_core, t.side_lobe_peak_y_rel_to_core, 'side lobe rel core', 'free threshold vs side lobe', base);
nexttile; plot(t.free_threshold, t.efficiency_inside_signal, 'ro-', t.free_threshold, t.efficiency_inside_guard, 'bo-', 'LineWidth', 1.2); grid on; xlabel('free threshold'); ylabel('efficiency'); legend('signal','guard'); title('free threshold vs efficiency'); xline(base, 'k--');
nexttile; axis off; text(0, 1, status_text(t), 'VerticalAlignment', 'top', 'FontName', 'Consolas');
sgtitle('Free/noise region probe: allow outer side lobe, protect shoulder');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end

function plot_pair(x, yx, yy, ylabel_text, title_text, base)
plot(x, yx, 'ro-', x, yy, 'bo-', 'LineWidth', 1.2); grid on;
xlabel('free threshold'); ylabel(ylabel_text); title(title_text); legend('x','y', 'Location', 'best'); xline(base, 'k--');
end

function txt = status_text(t)
txt = "threshold  status" + newline;
for i = 1:height(t)
    txt = txt + sprintf('%.3f      %s\n', t.free_threshold(i), t.status(i));
end
end

function plot_free_probe_profile_overlay(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1350 800]);
tiledlayout(1,2,'TileSpacing','compact');
ax1 = nexttile; ax2 = nexttile;
colors = lines(height(t));
for row_idx = 1:height(t)
    xt = readtable(fullfile(string(t.artifact_root(row_idx)), 'ideal_fft', 'x_profile_rd_target_mraf.csv'));
    yt = readtable(fullfile(string(t.artifact_root(row_idx)), 'ideal_fft', 'y_profile_rd_target_mraf.csv'));
    if row_idx == 1
        plot(ax1, xt.x_m*1e6, xt.rd_intensity_core_norm, 'k--', 'LineWidth', 1.3, 'DisplayName', 'RD baseline'); hold(ax1, 'on');
        plot(ax2, yt.y_m*1e6, yt.rd_intensity_core_norm, 'k--', 'LineWidth', 1.3, 'DisplayName', 'RD baseline'); hold(ax2, 'on');
    end
    if abs(t.free_threshold(row_idx) - 0.05) < 1e-12
        display_name = 'caseC baseline free=0.05'; line_width = 1.8;
    else
        display_name = sprintf('free=%.3f', t.free_threshold(row_idx)); line_width = 1.0;
    end
    plot(ax1, xt.x_m*1e6, xt.mraf_intensity_core_norm, '-', 'Color', colors(row_idx,:), 'LineWidth', line_width, 'DisplayName', display_name);
    plot(ax2, yt.y_m*1e6, yt.mraf_intensity_core_norm, '-', 'Color', colors(row_idx,:), 'LineWidth', line_width, 'DisplayName', display_name);
end
format_profile_axis(ax1, 'X center profiles: RD, caseC baseline, free-threshold probes', 'x (um)', [-380 380]);
format_profile_axis(ax2, 'Y center profiles: RD, caseC baseline, free-threshold probes', 'y (um)', [-230 230]);
sgtitle('Free/noise region probe profile overlay: inspect shoulder and acceptable side lobe');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end

function format_profile_axis(ax, title_text, xlabel_text, x_limits)
grid(ax, 'on'); xlim(ax, x_limits); ylim(ax, [0 1.35]);
yline(ax, 0.90, 'g:'); yline(ax, 0.50, 'r:'); yline(ax, 0.135, 'm:');
title(ax, title_text); xlabel(ax, xlabel_text); ylabel(ax, 'normalized intensity'); legend(ax, 'Location', 'northeastoutside');
end
