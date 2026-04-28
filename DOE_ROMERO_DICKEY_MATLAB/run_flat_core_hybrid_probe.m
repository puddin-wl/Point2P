function summary_root = run_flat_core_hybrid_probe()
% run_flat_core_hybrid_probe Hybrid of B flattening and C y-side control.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
artifact_group = "flat_core_hybrid_probe_" + stamp;

cases = {
    struct('name', "hybrid_Bwgs_Cy_outer050", 'free_x', 1.35, 'free_y', 1.35, 'outer_factor', 0.50)
    struct('name', "hybrid_Bwgs_Cy_outer035", 'free_x', 1.35, 'free_y', 1.35, 'outer_factor', 0.35)
};
rows = [];
for k = 1:numel(cases)
    c = cases{k};
    variant = struct('variant_name', c.name, 'preset_name', "flat_core_hybrid", 'artifact_group', artifact_group, ...
        'target_mode', "flat_core_free_edge", 'method', "WGS-MRAF", 'use_wgs', true, ...
        'use_three_region_projection', true, ...
        'flat_fraction_x', 0.80, 'flat_fraction_y', 0.80, ...
        'free_fraction_x', c.free_x, 'free_fraction_y', c.free_y, ...
        'wgs_exponent', 0.15, 'n_iter', 30, 'phase_blend', 0.18, ...
        'free_factor', 1.00, 'outer_factor', c.outer_factor);
    [cfg, result, fullprop] = run_mraf_one(variant); %#ok<ASGLU>
    copy_core_case_artifacts(cfg, result);
    m = result.mraf.metrics; rd = result.rd.metrics;
    rows = [rows; {k, string(c.name), string(cfg.save_root), c.free_x, c.free_y, c.outer_factor, ...
        rd.flat_core_rms, m.flat_core_rms, m.flat_core_pv, m.flat_core_uniformity, ...
        m.side_lobe_energy_ratio, m.outer_energy_ratio, m.size50_x_um, m.size50_y_um, ...
        m.transition_width_x_um, m.transition_width_y_um, m.core_rms, m.outer_tail_peak_y_rel_to_core, ...
        m.shoulder_peak_x, m.shoulder_peak_y}]; %#ok<AGROW>
end

summary_root = fullfile(project_root, 'artifacts', char(artifact_group));
if ~exist(summary_root, 'dir'), mkdir(summary_root); end
t = cell2table(rows, 'VariableNames', {'case_index','case_name','artifact_root','free_fraction_x','free_fraction_y','outer_factor', ...
    'rd_flat_core_rms','flat_core_rms','flat_core_pv','flat_core_uniformity','side_lobe_energy_ratio','outer_energy_ratio', ...
    'size50_x_um','size50_y_um','transition_width_x_um','transition_width_y_um','legacy_core_rms','outer_tail_y','shoulder_peak_x','shoulder_peak_y'});
t.flat_core_rms_improvement = t.rd_flat_core_rms - t.flat_core_rms;
t = sortrows(t, {'flat_core_rms','flat_core_pv'});
writetable(t, fullfile(summary_root, 'summary.csv'));
write_summary_text(fullfile(summary_root, 'summary.txt'), t);
plot_summary(fullfile(summary_root, 'summary_plot.png'), t);
fprintf('\nFlat-core hybrid probe summary written to:\n%s\n', summary_root);
end

function copy_core_case_artifacts(cfg, result)
metrics = result.mraf.metrics; %#ok<NASGU>
rd_metrics = result.rd.metrics; %#ok<NASGU>
save(fullfile(cfg.save_root, 'ideal_fft', 'metrics_mraf.mat'), 'metrics');
save(fullfile(cfg.save_root, 'ideal_fft', 'metrics_rd.mat'), 'rd_metrics');
copy_if_exists(fullfile(cfg.save_root, 'common', 'phase_mraf_wrapped.png'), fullfile(cfg.save_root, 'phase_mraf_wrapped.png'));
copy_if_exists(fullfile(cfg.save_root, 'ideal_fft', 'focal_mraf_intensity.png'), fullfile(cfg.save_root, 'focal_mraf_intensity.png'));
copy_if_exists(fullfile(cfg.save_root, 'ideal_fft', 'center_profiles_rd_vs_target_vs_mraf.png'), fullfile(cfg.save_root, 'center_profiles.png'));
copy_if_exists(fullfile(cfg.save_root, 'target_masks', 'flat_core_mask.png'), fullfile(cfg.save_root, 'flat_core_mask.png'));
copy_if_exists(fullfile(cfg.save_root, 'target_masks', 'free_edge_mask.png'), fullfile(cfg.save_root, 'free_edge_mask.png'));
end

function copy_if_exists(src, dst)
if exist(src, 'file'), copyfile(src, dst); end
end

function write_summary_text(path, t)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Flat-core hybrid probe summary\n==============================\n\n');
fprintf(fid, 'Hybrid = Case B WGS/n_iter/phase_blend + C-style tighter y free region. Rows sorted by flat_core_rms.\n\n');
for i = 1:height(t)
    fprintf(fid, '%s\n', t.case_name(i));
    fprintf(fid, '  artifact: %s\n', string(t.artifact_root(i)));
    fprintf(fid, '  free_fraction %.2f x %.2f, outer_factor %.2f\n', t.free_fraction_x(i), t.free_fraction_y(i), t.outer_factor(i));
    fprintf(fid, '  flat_core_rms %.6g (RD %.6g, improvement %.6g), pv %.6g, uniformity %.6g\n', t.flat_core_rms(i), t.rd_flat_core_rms(i), t.flat_core_rms_improvement(i), t.flat_core_pv(i), t.flat_core_uniformity(i));
    fprintf(fid, '  side_lobe_energy %.6g, outer_energy %.6g, size50 %.3f x %.3f um, transition %.3f x %.3f um\n', t.side_lobe_energy_ratio(i), t.outer_energy_ratio(i), t.size50_x_um(i), t.size50_y_um(i), t.transition_width_x_um(i), t.transition_width_y_um(i));
    fprintf(fid, '  shoulder %.6g x %.6g, outer_tail_y %.6g\n\n', t.shoulder_peak_x(i), t.shoulder_peak_y(i), t.outer_tail_y(i));
end
fprintf(fid, 'Recommended hybrid: %s\n', t.case_name(1));
end

function plot_summary(path, t)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 760]); tiledlayout(2,3,'TileSpacing','compact');
x = 1:height(t); labels = cellstr(t.case_name);
nexttile; bar(x, t.flat_core_rms); title('flat core RMS'); grid on;
nexttile; bar(x, t.flat_core_pv); title('flat core PV'); grid on;
nexttile; bar(x, [t.side_lobe_energy_ratio t.outer_energy_ratio]); title('energy ratios'); legend('free','outer'); grid on;
nexttile; plot(x, t.size50_x_um, 'ro-', x, t.size50_y_um, 'bo-', 'LineWidth', 1.2); title('size50'); legend('x','y'); grid on;
nexttile; plot(x, t.transition_width_x_um, 'ro-', x, t.transition_width_y_um, 'bo-', 'LineWidth', 1.2); title('transition width'); legend('x','y'); grid on;
nexttile; plot(x, t.shoulder_peak_x, 'ro-', x, t.shoulder_peak_y, 'bo-', 'LineWidth', 1.2); title('shoulder'); legend('x','y'); grid on;
for ax = findall(fig, 'Type', 'axes')'
    set(ax, 'XTick', x, 'XTickLabel', labels); xtickangle(ax, 25);
end
sgtitle('Flat-core hybrid probe');
exportgraphics(fig, path, 'Resolution', 150); close(fig);
end
