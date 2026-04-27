function plot_center_profiles_only(cfg, focal_x_m, focal_y_m, metrics, normalize_core, filename)
% plot_center_profiles_only 仅重新绘制并保存中心剖面诊断图。
%
% 该函数用于不重跑传播仿真的情况下，从已保存 metrics 重新生成
% center_profiles.png / center_profiles_norm.png。

if ~isfield(cfg, 'center_profile_xlim_factor')
    cfg.center_profile_xlim_factor = 1.8;
end

fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'compact');
plot_profile(nexttile, focal_x_m, metrics.x_profile, metrics.target_x_profile, cfg.target_half_x_m, metrics.crossing_x_50_m, metrics.crossing_x_13_m, metrics.crossing_x_90_m, 'x center profile', normalize_core, cfg.center_profile_xlim_factor);
plot_profile(nexttile, focal_y_m, metrics.y_profile, metrics.target_y_profile, cfg.target_half_y_m, metrics.crossing_y_50_m, metrics.crossing_y_13_m, metrics.crossing_y_90_m, 'y center profile', normalize_core, cfg.center_profile_xlim_factor);
exportgraphics(fig, fullfile(cfg.save_root, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end

function plot_profile(ax, axis_m, profile, target_profile, half_target_m, crossing_50_m, crossing_13_m, crossing_90_m, title_text, normalize_core, xlim_factor)
if normalize_core
    core = abs(axis_m) <= 0.5 * half_target_m;
    profile = profile ./ max(mean(profile(core), 'omitnan'), eps);
end
plot(ax, axis_m * 1e6, profile, 'b-', 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, axis_m * 1e6, target_profile, 'k--', 'LineWidth', 1.0);
yline(ax, 0.5, 'r:'); yline(ax, 0.13, 'm:'); yline(ax, 0.90, 'g:');
xline(ax, -half_target_m * 1e6, 'k:'); xline(ax, half_target_m * 1e6, 'k:');
draw_crossings(ax, crossing_50_m, 'r');
draw_crossings(ax, crossing_13_m, 'm');
draw_crossings(ax, crossing_90_m, 'g');
xlim(ax, [-xlim_factor * half_target_m, xlim_factor * half_target_m] * 1e6);
ylim(ax, [0, max(1.2, min(2, max(profile) * 1.05))]);
xlabel(ax, 'Focal coordinate (um)'); ylabel(ax, 'Intensity'); title(ax, title_text);
legend(ax, {'simulation', 'target', '50%', '13%', '90%'}, 'Location', 'best');
end

function draw_crossings(ax, crossings_m, color)
for k = 1:numel(crossings_m)
    xline(ax, crossings_m(k) * 1e6, '-', 'Color', color, 'Alpha', 0.35);
end
end

