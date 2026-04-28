function results = run_slmsuite_like_mraf_caseE()
% run_slmsuite_like_mraf_caseE Sweep slmsuite-like MRAF mraf_factor cases.
%
% CaseE_slmsuite_like_mraf uses amplitude target semantics:
%   signal rectangle: target amplitude = 1
%   MRAF noise ring: target amplitude = NaN, keep mraf_factor * current amp
%   null/background: target amplitude = 0

close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

cases = struct( ...
    'name', {"CaseE_mraf04", "CaseE_mraf06", "CaseE_mraf08"}, ...
    'mraf_factor', {0.4, 0.6, 0.8});

rows = cell(numel(cases), 1);
results = cell(numel(cases), 1);
for k = 1:numel(cases)
    variant = base_variant(cases(k).name, cases(k).mraf_factor);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    results{k} = struct('cfg', cfg, 'result', result, 'fullprop', fullprop);
    rows{k} = metrics_row(cases(k).name, cfg.save_root, cases(k).mraf_factor, result.mraf.metrics, fullprop.mraf_metrics);
end

summary = vertcat(rows{:});
summary_root = fullfile(project_root, 'artifacts', 'slmsuite_like_mraf', char(string(datetime('now', 'Format', 'yyyyMMdd-HHmmss')) + "_CaseE_summary"));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
writetable(summary, fullfile(summary_root, 'slmsuite_like_mraf_caseE_summary.csv'));
write_summary_text(fullfile(summary_root, 'slmsuite_like_mraf_caseE_summary.txt'), summary);
plot_summary(fullfile(summary_root, 'slmsuite_like_mraf_caseE_summary_plot.png'), summary);
fprintf('\nCaseE slmsuite-like summary folder:\n%s\n', summary_root);
fprintf('Summary CSV:\n%s\n', fullfile(summary_root, 'slmsuite_like_mraf_caseE_summary.csv'));
end

function variant = base_variant(name, mraf_factor)
variant = struct('variant_name', string(name), ...
    'preset_name', "CaseE_slmsuite_like_mraf", ...
    'artifact_group', "slmsuite_like_mraf", ...
    'target_mode', "slmsuite_like_mraf", ...
    'method', "WGS-MRAF", ...
    'use_wgs', true, ...
    'wgs_exponent', 0.08, ...
    'n_iter', 20, ...
    'phase_blend', 0.20, ...
    'use_three_region_projection', true, ...
    'mraf_factor', mraf_factor);
end

function row = metrics_row(name, artifact_path, mraf_factor, m, fm)
row = table(string(name), string(artifact_path), mraf_factor, ...
    m.signal_rms, m.signal_pv, m.shoulder_peak, m.noise_energy_ratio, m.null_energy_ratio, ...
    m.measured_size50_x_um, m.measured_size50_y_um, m.measured_size135_x_um, m.measured_size135_y_um, m.e2_efficiency, ...
    fm.signal_rms, fm.signal_pv, fm.shoulder_peak, fm.noise_energy_ratio, fm.null_energy_ratio, fm.e2_efficiency, ...
    'VariableNames', {'case_name','artifact_path','mraf_factor', ...
    'ideal_signal_rms','ideal_signal_pv','ideal_shoulder_peak','ideal_noise_energy_ratio','ideal_null_energy_ratio', ...
    'ideal_measured_size50_x_um','ideal_measured_size50_y_um','ideal_measured_size135_x_um','ideal_measured_size135_y_um','ideal_e2_efficiency', ...
    'full_signal_rms','full_signal_pv','full_shoulder_peak','full_noise_energy_ratio','full_null_energy_ratio','full_e2_efficiency'});
end

function write_summary_text(path, summary)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'CaseE slmsuite-like MRAF mraf_factor sweep\n');
fprintf(fid, '==========================================\n\n');
for k = 1:height(summary)
    fprintf(fid, '%s\n', summary.case_name(k));
    fprintf(fid, '  artifact = %s\n', summary.artifact_path(k));
    fprintf(fid, '  mraf_factor = %.3f\n', summary.mraf_factor(k));
    fprintf(fid, '  signal_rms = %.8g, signal_pv = %.8g\n', summary.ideal_signal_rms(k), summary.ideal_signal_pv(k));
    fprintf(fid, '  shoulder_peak = %.8g\n', summary.ideal_shoulder_peak(k));
    fprintf(fid, '  noise_energy_ratio = %.8g, null_energy_ratio = %.8g\n', summary.ideal_noise_energy_ratio(k), summary.ideal_null_energy_ratio(k));
    fprintf(fid, '  size50 = %.4f x %.4f um, size13.5 = %.4f x %.4f um\n', summary.ideal_measured_size50_x_um(k), summary.ideal_measured_size50_y_um(k), summary.ideal_measured_size135_x_um(k), summary.ideal_measured_size135_y_um(k));
    fprintf(fid, '  e2_efficiency = %.8g\n\n', summary.ideal_e2_efficiency(k));
end
[~, best_idx] = min(summary.ideal_signal_rms + 0.25 .* summary.ideal_shoulder_peak);
fprintf(fid, 'Heuristic candidate = %s (balances signal_rms and shoulder_peak; inspect images before accepting).\n', summary.case_name(best_idx));
end

function plot_summary(path, summary)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 780]);
tiledlayout(2,3,'TileSpacing','compact');
x = summary.mraf_factor;
nexttile; plot(x, summary.ideal_signal_rms, '-o', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('signal RMS');
nexttile; plot(x, summary.ideal_signal_pv, '-o', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('signal PV');
nexttile; plot(x, summary.ideal_shoulder_peak, '-o', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('shoulder peak in noise');
nexttile; plot(x, summary.ideal_noise_energy_ratio, '-o', x, summary.ideal_null_energy_ratio, '-s', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('energy / signal'); legend('noise','null');
nexttile; plot(x, summary.ideal_measured_size50_x_um, '-o', x, summary.ideal_measured_size50_y_um, '-s', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('size50 (um)'); legend('x','y');
nexttile; plot(x, summary.ideal_e2_efficiency, '-o', 'LineWidth', 1.3); grid on; xlabel('mraf\_factor'); ylabel('e2 efficiency');
sgtitle('CaseE slmsuite-like MRAF sweep');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end
