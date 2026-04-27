function plot_diagnostics(cfg, focal_x_m, focal_y_m, focal_intensity, target_intensity, ...
    input_amplitude, input_intensity, aperture_mask, phase_wrapped_rad, metrics, prefix)
% plot_diagnostics 保存输入、相位、焦平面、ROI 和中心剖面诊断图。

if nargin < 10
    prefix = '';
end

save_image(cfg, input_amplitude, [prefix 'input_amplitude.png'], 'Input amplitude', 'DOE x (mm)', 'DOE y (mm)', cfg.doe_grid_extent_m * 1e3);
save_image(cfg, input_intensity, [prefix 'input_intensity.png'], 'Input intensity', 'DOE x (mm)', 'DOE y (mm)', cfg.doe_grid_extent_m * 1e3);
save_image(cfg, double(aperture_mask), [prefix 'aperture_mask.png'], '15 mm clear aperture mask', 'DOE x (mm)', 'DOE y (mm)', cfg.doe_grid_extent_m * 1e3);
save_image(cfg, phase_wrapped_rad, [prefix 'phase_wrapped.png'], 'Wrapped phase (rad)', 'DOE x (mm)', 'DOE y (mm)', cfg.doe_grid_extent_m * 1e3);
save_phase_overlay(cfg, input_amplitude, phase_wrapped_rad, prefix);

save_focal_image(cfg, focal_x_m, focal_y_m, focal_intensity, [prefix 'focal_intensity.png'], false);
save_focal_image(cfg, focal_x_m, focal_y_m, focal_intensity, [prefix 'focal_intensity_log.png'], true);
save_focal_image(cfg, focal_x_m, focal_y_m, target_intensity, [prefix 'target.png'], false);
save_roi_image(cfg, focal_x_m, focal_y_m, focal_intensity, target_intensity, prefix);
plot_center_profiles_only(cfg, focal_x_m, focal_y_m, metrics, false, [prefix 'center_profiles.png']);
plot_center_profiles_only(cfg, focal_x_m, focal_y_m, metrics, true, [prefix 'center_profiles_norm.png']);
save_edge_profiles(cfg, focal_x_m, focal_y_m, metrics, prefix);
end

function save_image(cfg, data, filename, title_text, xlabel_text, ylabel_text, extent_mm)
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(linspace(-extent_mm/2, extent_mm/2, size(data, 2)), linspace(-extent_mm/2, extent_mm/2, size(data, 1)), data);
axis image; colormap(gca, 'turbo'); colorbar;
xlabel(xlabel_text); ylabel(ylabel_text); title(title_text);
exportgraphics(fig, fullfile(cfg.save_root, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end

function save_phase_overlay(cfg, input_amplitude, phase_wrapped_rad, prefix)
fig = figure('Visible', 'off', 'Color', 'w');
extent_mm = cfg.doe_grid_extent_m * 1e3;
imagesc(linspace(-extent_mm/2, extent_mm/2, size(phase_wrapped_rad, 2)), linspace(-extent_mm/2, extent_mm/2, size(phase_wrapped_rad, 1)), phase_wrapped_rad);
axis image; colormap(gca, 'hsv'); colorbar; hold on;
contour(linspace(-extent_mm/2, extent_mm/2, size(input_amplitude, 2)), linspace(-extent_mm/2, extent_mm/2, size(input_amplitude, 1)), input_amplitude, [exp(-1) exp(-1)], 'w', 'LineWidth', 1.2);
title('Wrapped phase with 1/e amplitude contour');
xlabel('DOE x (mm)'); ylabel('DOE y (mm)');
exportgraphics(fig, fullfile(cfg.save_root, [prefix 'phase_with_beam_overlay.png']), 'Resolution', cfg.figure_dpi);
close(fig);
end

function save_focal_image(cfg, focal_x_m, focal_y_m, intensity, filename, use_log)
fig = figure('Visible', 'off', 'Color', 'w');
plot_data = intensity;
if use_log
    plot_data = log10(max(intensity, 1e-8));
end
imagesc(focal_x_m * 1e6, focal_y_m * 1e6, plot_data);
axis image; colormap(gca, 'turbo'); colorbar; hold on;
rectangle('Position', [-cfg.target_half_x_m, -cfg.target_half_y_m, cfg.target_size_x_m, cfg.target_size_y_m] * 1e6, 'EdgeColor', 'w', 'LineWidth', 1.2);
xlim([-600 600]); ylim([-400 400]);
xlabel('Focal x (um)'); ylabel('Focal y (um)');
if use_log
    title('Focal intensity log10 normalized');
else
    title('Focal intensity normalized');
end
exportgraphics(fig, fullfile(cfg.save_root, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end

function save_roi_image(cfg, focal_x_m, focal_y_m, intensity, target_intensity, prefix)
fig = figure('Visible', 'off', 'Color', 'w');
roi = abs(focal_x_m) <= cfg.guard_fraction * cfg.target_half_x_m & true;
roi_y = abs(focal_y_m) <= cfg.guard_fraction * cfg.target_half_y_m;
imagesc(focal_x_m(roi) * 1e6, focal_y_m(roi_y) * 1e6, intensity(roi_y, roi));
axis image; colormap(gca, 'turbo'); colorbar; hold on;
contour(focal_x_m(roi) * 1e6, focal_y_m(roi_y) * 1e6, target_intensity(roi_y, roi), [0.5 0.5], 'w', 'LineWidth', 1.2);
title('ROI intensity with target 50% contour');
xlabel('Focal x (um)'); ylabel('Focal y (um)');
exportgraphics(fig, fullfile(cfg.save_root, [prefix 'roi_intensity.png']), 'Resolution', cfg.figure_dpi);
close(fig);
end

function plot_profile(ax, axis_m, profile, target_profile, half_target_m, crossing_50_m, crossing_13_m, crossing_90_m, title_text, normalize_core)
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
xlim_factor = 1.8;
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

function save_edge_profiles(cfg, focal_x_m, focal_y_m, metrics, prefix)
fig = figure('Visible', 'off', 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact');
ax = nexttile;
plot(ax, focal_x_m * 1e6, metrics.x_profile, 'b-', 'LineWidth', 1.2); hold(ax, 'on');
xline(ax, cfg.target_half_x_m * 1e6, 'k:'); xline(ax, -cfg.target_half_x_m * 1e6, 'k:');
yline(ax, 0.13, 'm:'); yline(ax, 0.90, 'g:'); yline(ax, 0.5, 'r:');
xlim(ax, [-2.2 2.2] * cfg.target_half_x_m * 1e6); title(ax, 'x edge diagnostics');
xlabel(ax, 'x (um)'); ylabel(ax, 'Intensity');
ax = nexttile;
plot(ax, focal_y_m * 1e6, metrics.y_profile, 'b-', 'LineWidth', 1.2); hold(ax, 'on');
xline(ax, cfg.target_half_y_m * 1e6, 'k:'); xline(ax, -cfg.target_half_y_m * 1e6, 'k:');
yline(ax, 0.13, 'm:'); yline(ax, 0.90, 'g:'); yline(ax, 0.5, 'r:');
xlim(ax, [-2.2 2.2] * cfg.target_half_y_m * 1e6); title(ax, 'y edge diagnostics');
xlabel(ax, 'y (um)'); ylabel(ax, 'Intensity');
exportgraphics(fig, fullfile(cfg.save_root, [prefix 'edge_diagnostic_profiles.png']), 'Resolution', cfg.figure_dpi);
close(fig);
end
