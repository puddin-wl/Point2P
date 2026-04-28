function plot_mraf_diagnostics(cfg, focal_x_m, focal_y_m, rd, target, mraf, fullprop)
% plot_mraf_diagnostics Save MRAF diagnostics into artifact subfolders.
save_phase(cfg, rd.phase_wrapped, 'common', 'phase_rd_wrapped.png', 'RD wrapped phase');
save_phase(cfg, mraf.phase_wrapped, 'common', 'phase_mraf_wrapped.png', 'MRAF wrapped phase');
save_image(cfg, angle(exp(1i * (mraf.phase - rd.phase))), 'common', 'phase_delta_mraf_minus_rd.png', 'MRAF minus RD phase', 'hsv');

save_focal(cfg, focal_x_m, focal_y_m, rd.intensity_norm, rd.metrics, 'target_masks', 'rd_baseline_intensity.png', false);
if strcmpi(string(target.mode), "hard_rectangle")
    save_focal(cfg, focal_x_m, focal_y_m, target.amplitude, [], 'target_masks', 'hard_rectangle_target_amplitude.png', false);
    save_focal(cfg, focal_x_m, focal_y_m, target.intensity_plot, [], 'target_masks', 'hard_rectangle_target_intensity.png', false);
    save_mask(cfg, target.masks.signal, 'target_masks', 'hard_rectangle_signal_mask.png', 'hard rectangle signal mask');
    save_mask(cfg, target.masks.free, 'target_masks', 'hard_rectangle_free_mask.png', 'hard rectangle free mask');
else
    save_focal(cfg, focal_x_m, focal_y_m, target.amplitude, [], 'target_masks', 'rd_derived_target_amplitude.png', false);
    save_focal(cfg, focal_x_m, focal_y_m, target.intensity_plot, [], 'target_masks', 'rd_derived_target_intensity.png', false);
end
save_focal(cfg, focal_x_m, focal_y_m, target.blend_map, [], 'target_masks', 'blend_map.png', false);
save_mask(cfg, target.masks.core, 'target_masks', 'core_mask.png', 'core mask');
save_mask(cfg, target.masks.signal, 'target_masks', 'signal_mask.png', 'signal mask');
save_mask(cfg, target.masks.edge, 'target_masks', 'edge_mask.png', 'edge mask');
save_mask(cfg, target.masks.noise, 'target_masks', 'noise_mask.png', 'noise/free mask');
save_mask(cfg, target.masks.free, 'target_masks', 'free_mask.png', 'free mask');
if isfield(target.masks, 'target')
    save_mask(cfg, target.masks.target, 'target_masks', 'target_mask.png', 'three-region target mask');
end
if isfield(target.masks, 'free_noise')
    save_mask(cfg, target.masks.free_noise, 'target_masks', 'free_noise_mask.png', 'three-region free/noise mask');
end
if isfield(target.masks, 'outer_suppress')
    save_mask(cfg, target.masks.outer_suppress, 'target_masks', 'outer_suppress_mask.png', 'three-region outer suppress mask');
end
if isfield(target.masks, 'y_free_reservoir')
    save_mask(cfg, target.masks.y_free_reservoir, 'target_masks', 'y_free_reservoir_mask.png', 'y free reservoir mask');
end
if isfield(target.masks, 'wgs_y_edge_damping_band')
    save_mask(cfg, target.masks.wgs_y_edge_damping_band, 'target_masks', 'wgs_y_edge_damping_band.png', 'WGS y edge damping band');
end
save_remap_curve(cfg, target, 'target_masks', 'pivot50_remap_curve.png');

save_focal(cfg, focal_x_m, focal_y_m, rd.intensity_norm, rd.metrics, 'ideal_fft', 'focal_rd_intensity.png', false);
save_focal(cfg, focal_x_m, focal_y_m, mraf.intensity_norm, mraf.metrics, 'ideal_fft', 'focal_mraf_intensity.png', false);
save_focal(cfg, focal_x_m, focal_y_m, mraf.intensity_norm, mraf.metrics, 'ideal_fft', 'focal_mraf_log.png', true);
save_same_scale(cfg, focal_x_m, focal_y_m, rd.intensity_norm, mraf.intensity_norm, 'ideal_fft', 'focal_rd_vs_mraf_same_scale.png', 'Full focal RD vs MRAF');
save_roi(cfg, focal_x_m, focal_y_m, rd.intensity_norm, rd.metrics, 'ideal_fft', 'roi_rd_intensity.png');
save_roi(cfg, focal_x_m, focal_y_m, mraf.intensity_norm, mraf.metrics, 'ideal_fft', 'roi_mraf_intensity.png');
save_roi_same_scale(cfg, focal_x_m, focal_y_m, rd.intensity_norm, mraf.intensity_norm, rd.metrics, mraf.metrics, 'ideal_fft', 'roi_rd_vs_mraf_same_scale.png');
save_profiles(cfg, focal_x_m, focal_y_m, rd.metrics, target, mraf.metrics, 'ideal_fft', 'center_profiles_rd_vs_target_vs_mraf.png');
save_edge_profiles(cfg, focal_x_m, focal_y_m, rd.metrics, target, mraf.metrics, 'ideal_fft', 'edge_profiles_rd_vs_target_vs_mraf.png');
save_summary(cfg, focal_x_m, focal_y_m, rd.intensity_norm, mraf.intensity_norm, rd.metrics, mraf.metrics, target, 'ideal_fft', 'summary_rd_vs_mraf.png', 'Ideal FFT');

if nargin >= 7 && ~isempty(fullprop)
    save_focal(cfg, fullprop.focal_x_m, fullprop.focal_y_m, fullprop.rd_intensity_norm, fullprop.rd_metrics, 'full_propagation', 'fullprop_rd_intensity.png', false);
    save_focal(cfg, fullprop.focal_x_m, fullprop.focal_y_m, fullprop.mraf_intensity_norm, fullprop.mraf_metrics, 'full_propagation', 'fullprop_mraf_intensity.png', false);
    save_full_profiles(fullprop, 'full_propagation', 'fullprop_center_profiles_rd_vs_mraf.png');
    save_summary(cfg, fullprop.focal_x_m, fullprop.focal_y_m, fullprop.rd_intensity_norm, fullprop.mraf_intensity_norm, fullprop.rd_metrics, fullprop.mraf_metrics, [], 'full_propagation', 'summary_fullprop_rd_vs_mraf.png', 'Full propagation');
end
end

function path = out_path(cfg, folder_key, filename)
path = fullfile(cfg.artifact_dirs.(folder_key), filename);
end

function save_phase(cfg, phase_wrapped, folder_key, filename, title_text)
plot_data = phase_wrapped; plot_data(~isfinite(plot_data)) = NaN;
save_image(cfg, plot_data, folder_key, filename, title_text, 'hsv');
end

function save_image(cfg, data, folder_key, filename, title_text, cmap)
if nargin < 6, cmap = 'turbo'; end
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(data); axis image off; colormap(gca, cmap); colorbar; title(title_text, 'Interpreter', 'none');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end


function save_remap_curve(cfg, target, folder_key, filename)
fig = figure('Visible', 'off', 'Color', 'w');
plot(target.remap_curve_x, target.remap_curve_y, 'b-', 'LineWidth', 1.5); hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.0);
plot(0.5, 0.5, 'ro', 'MarkerFaceColor', 'r');
grid on; axis([0 1 0 1]); xlabel('I_{rd,norm}'); ylabel('I_{env}');
title(sprintf('target edge remap: %s, inner=%.3f, outer=%.3f, remap(0.5)=%.3f', string(target.edge_mode), target.pivot50_inner_power, target.pivot50_outer_power, target.remap_anchor_50), 'Interpreter', 'none');
legend('remap', 'identity', '0.5 anchor', 'Location', 'best');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function save_mask(cfg, mask, folder_key, filename, title_text)
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(double(mask)); axis image off; colormap(gca, 'gray'); colorbar; title(title_text);
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end

function save_focal(cfg, focal_x_m, focal_y_m, intensity, metrics, folder_key, filename, use_log)
plot_data = intensity; plot_data(~isfinite(plot_data)) = 0;
if use_log, plot_data = log10(max(plot_data, 1e-8)); end
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(focal_x_m * 1e6, focal_y_m * 1e6, plot_data); axis image; colormap(gca, 'turbo'); colorbar; hold on;
if ~isempty(metrics), draw_size_boxes(metrics); title(metric_title(filename, metrics), 'Interpreter', 'none'); else, title(filename, 'Interpreter', 'none'); end
xlim([-600 600]); ylim([-400 400]); xlabel('x (um)'); ylabel('y (um)');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi);
close(fig);
end

function draw_size_boxes(metrics)
if isfinite(metrics.size50_x_um) && isfinite(metrics.size50_y_um)
    rectangle('Position', [-metrics.size50_x_um/2, -metrics.size50_y_um/2, metrics.size50_x_um, metrics.size50_y_um], 'EdgeColor', 'w', 'LineWidth', 1.3, 'LineStyle', '-');
end
if isfinite(metrics.size13p5_x_um) && isfinite(metrics.size13p5_y_um)
    rectangle('Position', [-metrics.size13p5_x_um/2, -metrics.size13p5_y_um/2, metrics.size13p5_x_um, metrics.size13p5_y_um], 'EdgeColor', 'm', 'LineWidth', 1.2, 'LineStyle', '--');
end
end

function txt = metric_title(filename, metrics)
txt = sprintf('%s | Size50 %.1f x %.1f um | Size13.5 %.1f x %.1f um | TW13.5-90 %.1f x %.1f um', filename, metrics.size50_x_um, metrics.size50_y_um, metrics.size13p5_x_um, metrics.size13p5_y_um, metrics.transition_13p5_90_x_um, metrics.transition_13p5_90_y_um);
end

function save_same_scale(cfg, focal_x_m, focal_y_m, I1, I2, folder_key, filename, title_text)
fig = figure('Visible', 'off', 'Color', 'w'); tiledlayout(1,2,'TileSpacing','compact');
maxv = max([I1(:); I2(:)]);
nexttile; imagesc(focal_x_m*1e6, focal_y_m*1e6, I1, [0 maxv]); axis image; xlim([-600 600]); ylim([-400 400]); title('RD'); colorbar;
nexttile; imagesc(focal_x_m*1e6, focal_y_m*1e6, I2, [0 maxv]); axis image; xlim([-600 600]); ylim([-400 400]); title('MRAF'); colorbar; sgtitle(title_text);
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function save_roi(cfg, focal_x_m, focal_y_m, I, metrics, folder_key, filename)
roi_x = abs(focal_x_m) <= cfg.guard_fraction * cfg.target_half_x_m;
roi_y = abs(focal_y_m) <= cfg.guard_fraction * cfg.target_half_y_m;
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, I(roi_y, roi_x)); axis image; colormap(gca, 'turbo'); colorbar; hold on;
draw_size_boxes(metrics);
title(metric_title(filename, metrics), 'Interpreter', 'none'); xlabel('x (um)'); ylabel('y (um)');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function save_roi_same_scale(cfg, focal_x_m, focal_y_m, I1, I2, rd_metrics, mraf_metrics, folder_key, filename)
roi_x = abs(focal_x_m) <= cfg.guard_fraction * cfg.target_half_x_m;
roi_y = abs(focal_y_m) <= cfg.guard_fraction * cfg.target_half_y_m;
fig = figure('Visible', 'off', 'Color', 'w'); tiledlayout(1,2,'TileSpacing','compact');
maxv = max([I1(roi_y, roi_x); I2(roi_y, roi_x)], [], 'all');
nexttile; imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, I1(roi_y, roi_x), [0 maxv]); axis image; colormap(gca, 'turbo'); colorbar; hold on; draw_size_boxes(rd_metrics); title('RD ROI');
nexttile; imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, I2(roi_y, roi_x), [0 maxv]); axis image; colormap(gca, 'turbo'); colorbar; hold on; draw_size_boxes(mraf_metrics); title(metric_title('MRAF ROI', mraf_metrics), 'Interpreter', 'none');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function save_profiles(cfg, focal_x_m, focal_y_m, rd_metrics, target, mraf_metrics, folder_key, filename)
center_x = round(numel(focal_x_m)/2)+1; center_y = round(numel(focal_y_m)/2)+1;
tx = target.intensity(center_y, :); ty = target.intensity(:, center_x).';
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 760]); tiledlayout(2,1,'TileSpacing','compact');
ax = nexttile; plot_profile_with_annotations(ax, focal_x_m, rd_metrics.x_profile, mraf_metrics.x_profile, tx, rd_metrics, mraf_metrics, 'x');
ax = nexttile; plot_profile_with_annotations(ax, focal_y_m, rd_metrics.y_profile, mraf_metrics.y_profile, ty, rd_metrics, mraf_metrics, 'y');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function plot_profile_with_annotations(ax, axis_m, rd_profile, mraf_profile, target_profile, rd_metrics, mraf_metrics, direction)
plot(ax, axis_m*1e6, rd_profile, 'b-', axis_m*1e6, mraf_profile, 'r-', 'LineWidth', 1.2); hold(ax, 'on');
if ~isempty(target_profile), plot(ax, axis_m*1e6, target_profile, 'k--', 'LineWidth', 1.0); end
yline(ax, 0.90, 'g:'); yline(ax, 0.50, 'r:'); yline(ax, 0.135, 'm:');
draw_crossing_lines(ax, mraf_metrics, direction);
grid(ax, 'on'); xlabel(ax, sprintf('%s (um)', direction)); ylabel(ax, 'Normalized intensity');
if direction == 'x', xlim(ax, [-400 400]); shoulder = mraf_metrics.shoulder_peak_x; else, xlim(ax, [-220 220]); shoulder = mraf_metrics.shoulder_peak_y; end
draw_true_side_lobe_marker(ax, mraf_metrics, direction);
ylim(ax, [0 1.45]); title(ax, sprintf('%s center profile | shoulder %.4f | core_rms %.4f', upper(direction), shoulder, mraf_metrics.core_rms));
legend(ax, {'RD baseline','MRAF result','target','90%','50%','13.5%'}, 'Location', 'northeastoutside');
text(ax, 0.02, 0.95, profile_text(rd_metrics, mraf_metrics, direction), 'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', [0.4 0.4 0.4], 'FontName', 'Microsoft YaHei', 'FontSize', 9);
end

function draw_crossing_lines(ax, metrics, direction)
if direction == 'x'
    vals = [metrics.crossing50_x_left_um metrics.crossing50_x_right_um metrics.crossing13p5_x_left_um metrics.crossing13p5_x_right_um metrics.crossing90_x_left_um metrics.crossing90_x_right_um];
else
    vals = [metrics.crossing50_y_left_um metrics.crossing50_y_right_um metrics.crossing13p5_y_left_um metrics.crossing13p5_y_right_um metrics.crossing90_y_left_um metrics.crossing90_y_right_um];
end
colors = {'r','r','m','m','g','g'};
for k = 1:numel(vals)
    if isfinite(vals(k)), xline(ax, vals(k), '-', 'Color', colors{k}, 'Alpha', 0.35); end
end
end

function draw_true_side_lobe_marker(ax, metrics, direction)
if direction == 'x'
    has_peak = metrics.has_true_side_lobe_x > 0;
    peak_pos = metrics.true_side_lobe_pos_x_um;
    peak_val = metrics.true_side_lobe_peak_x_rel_to_core;
else
    has_peak = metrics.has_true_side_lobe_y > 0;
    peak_pos = metrics.true_side_lobe_pos_y_um;
    peak_val = metrics.true_side_lobe_peak_y_rel_to_core;
end
if has_peak && isfinite(peak_pos) && isfinite(peak_val)
    plot(ax, peak_pos, peak_val, 'vp', 'MarkerSize', 9, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'true side lobe');
end
end

function txt = profile_text(rd, mraf, direction)
if direction == 'x'
    if mraf.has_true_side_lobe_x > 0
        true_side_lobe_text = sprintf('true_side_lobe_x = %.4f @ %.1f um', mraf.true_side_lobe_peak_x_rel_to_core, mraf.true_side_lobe_pos_x_um);
    else
        true_side_lobe_text = 'true_side_lobe_x = none';
    end
    txt = sprintf(['X direction\n', ...
        '匀化光斑尺寸（50%%）= %.1f um\n', ...
        '匀化光斑尺寸（13.5%%）= %.1f um\n', ...
        '传输区宽度（13.5%%-90%%）= %.1f um\n', ...
        'shoulder_peak_x = %.4f (RD %.4f)\n', ...
        'outer_tail_x = %.4f\n%s\ncore_rms = %.4f (RD %.4f)'], ...
        mraf.size50_x_um, mraf.size13p5_x_um, mraf.transition_13p5_90_x_um, mraf.shoulder_peak_x, rd.shoulder_peak_x, ...
        mraf.outer_tail_peak_x_rel_to_core, true_side_lobe_text, mraf.core_rms, rd.core_rms);
else
    if mraf.has_true_side_lobe_y > 0
        true_side_lobe_text = sprintf('true_side_lobe_y = %.4f @ %.1f um', mraf.true_side_lobe_peak_y_rel_to_core, mraf.true_side_lobe_pos_y_um);
    else
        true_side_lobe_text = 'true_side_lobe_y = none';
    end
    txt = sprintf(['Y direction\n', ...
        '匀化光斑尺寸（50%%）= %.1f um\n', ...
        '匀化光斑尺寸（13.5%%）= %.1f um\n', ...
        '传输区宽度（13.5%%-90%%）= %.1f um\n', ...
        'shoulder_peak_y = %.4f (RD %.4f)\n', ...
        'outer_tail_y = %.4f\n%s\ncore_rms = %.4f (RD %.4f)'], ...
        mraf.size50_y_um, mraf.size13p5_y_um, mraf.transition_13p5_90_y_um, mraf.shoulder_peak_y, rd.shoulder_peak_y, ...
        mraf.outer_tail_peak_y_rel_to_core, true_side_lobe_text, mraf.core_rms, rd.core_rms);
end
end

function save_edge_profiles(cfg, focal_x_m, focal_y_m, rd_metrics, target, mraf_metrics, folder_key, filename)
save_profiles(cfg, focal_x_m, focal_y_m, rd_metrics, target, mraf_metrics, folder_key, filename);
end

function save_full_profiles(fullprop, folder_key, filename)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 760]); tiledlayout(2,1,'TileSpacing','compact');
ax = nexttile; plot_profile_with_annotations(ax, fullprop.focal_x_m, fullprop.rd_metrics.x_profile, fullprop.mraf_metrics.x_profile, [], fullprop.rd_metrics, fullprop.mraf_metrics, 'x');
ax = nexttile; plot_profile_with_annotations(ax, fullprop.focal_y_m, fullprop.rd_metrics.y_profile, fullprop.mraf_metrics.y_profile, [], fullprop.rd_metrics, fullprop.mraf_metrics, 'y');
exportgraphics(fig, out_path(fullprop.cfg, folder_key, filename), 'Resolution', fullprop.cfg.figure_dpi); close(fig);
end

function save_summary(cfg, focal_x_m, focal_y_m, Ird, Imraf, rd_metrics, mraf_metrics, target, folder_key, filename, title_prefix)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1500 950]); tiledlayout(2,3,'TileSpacing','compact');
roi_x = abs(focal_x_m) <= cfg.guard_fraction * cfg.target_half_x_m;
roi_y = abs(focal_y_m) <= cfg.guard_fraction * cfg.target_half_y_m;
maxv = max([Ird(roi_y, roi_x); Imraf(roi_y, roi_x)], [], 'all');
nexttile; imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, Ird(roi_y, roi_x), [0 maxv]); axis image; colorbar; title('RD ROI'); hold on; draw_size_boxes(rd_metrics);
nexttile; imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, Imraf(roi_y, roi_x), [0 maxv]); axis image; colorbar; title('MRAF ROI'); hold on; draw_size_boxes(mraf_metrics);
ax = nexttile; axis(ax, 'off'); text(ax, 0, 1, summary_text(rd_metrics, mraf_metrics), 'VerticalAlignment', 'top', 'FontName', 'Consolas', 'FontSize', 9);
center_x = round(numel(focal_x_m)/2)+1; center_y = round(numel(focal_y_m)/2)+1;
if isempty(target), tx = []; ty = []; else, tx = target.intensity(center_y, :); ty = target.intensity(:, center_x).'; end
ax = nexttile; plot_profile_with_annotations(ax, focal_x_m, rd_metrics.x_profile, mraf_metrics.x_profile, tx, rd_metrics, mraf_metrics, 'x');
ax = nexttile; plot_profile_with_annotations(ax, focal_y_m, rd_metrics.y_profile, mraf_metrics.y_profile, ty, rd_metrics, mraf_metrics, 'y');
ax = nexttile; plot(ax, focal_x_m*1e6, rd_metrics.x_profile - mraf_metrics.x_profile, 'k-', focal_y_m*1e6, rd_metrics.y_profile - mraf_metrics.y_profile, 'c-', 'LineWidth', 1.1); grid(ax, 'on'); xlim(ax, [-350 350]); title(ax, 'RD - MRAF profile delta'); legend(ax, {'x','y'});
sgtitle(sprintf('%s RD vs MRAF summary', title_prefix), 'Interpreter', 'none');
exportgraphics(fig, out_path(cfg, folder_key, filename), 'Resolution', cfg.figure_dpi); close(fig);
end

function txt = summary_text(rd, m)
txt = sprintf(['Metric                     RD          MRAF        Delta\n', ...
    'size50 x/y       %7.1f/%-7.1f %7.1f/%-7.1f %7.1f/%-7.1f\n', ...
    'size13.5 x/y     %7.1f/%-7.1f %7.1f/%-7.1f %7.1f/%-7.1f\n', ...
    'TW13.5-90 x/y    %7.1f/%-7.1f %7.1f/%-7.1f %7.1f/%-7.1f\n', ...
    'core_rms         %12.5g %12.5g %12.5g\n', ...
    'shoulder x/y     %7.4f/%-7.4f %7.4f/%-7.4f %7.4f/%-7.4f\n', ...
    'outer tail x/y   %7.4f/%-7.4f %7.4f/%-7.4f %7.4f/%-7.4f\n', ...
    'true SL x/y      %7.4f/%-7.4f %7.4f/%-7.4f %7.4f/%-7.4f\n', ...
    'eff signal/guard %7.4f/%-7.4f %7.4f/%-7.4f %7.4f/%-7.4f'], ...
    rd.size50_x_um, rd.size50_y_um, m.size50_x_um, m.size50_y_um, m.size50_x_um-rd.size50_x_um, m.size50_y_um-rd.size50_y_um, ...
    rd.size13p5_x_um, rd.size13p5_y_um, m.size13p5_x_um, m.size13p5_y_um, m.size13p5_x_um-rd.size13p5_x_um, m.size13p5_y_um-rd.size13p5_y_um, ...
    rd.transition_13p5_90_x_um, rd.transition_13p5_90_y_um, m.transition_13p5_90_x_um, m.transition_13p5_90_y_um, m.transition_13p5_90_x_um-rd.transition_13p5_90_x_um, m.transition_13p5_90_y_um-rd.transition_13p5_90_y_um, ...
    rd.core_rms, m.core_rms, m.core_rms-rd.core_rms, ...
    rd.shoulder_peak_x, rd.shoulder_peak_y, m.shoulder_peak_x, m.shoulder_peak_y, m.shoulder_peak_x-rd.shoulder_peak_x, m.shoulder_peak_y-rd.shoulder_peak_y, ...
    rd.outer_tail_peak_x_rel_to_core, rd.outer_tail_peak_y_rel_to_core, m.outer_tail_peak_x_rel_to_core, m.outer_tail_peak_y_rel_to_core, m.outer_tail_peak_x_rel_to_core-rd.outer_tail_peak_x_rel_to_core, m.outer_tail_peak_y_rel_to_core-rd.outer_tail_peak_y_rel_to_core, ...
    rd.true_side_lobe_peak_x_rel_to_core, rd.true_side_lobe_peak_y_rel_to_core, m.true_side_lobe_peak_x_rel_to_core, m.true_side_lobe_peak_y_rel_to_core, m.true_side_lobe_peak_x_rel_to_core-rd.true_side_lobe_peak_x_rel_to_core, m.true_side_lobe_peak_y_rel_to_core-rd.true_side_lobe_peak_y_rel_to_core, ...
    rd.efficiency_inside_signal, rd.efficiency_inside_guard, m.efficiency_inside_signal, m.efficiency_inside_guard, m.efficiency_inside_signal-rd.efficiency_inside_signal, m.efficiency_inside_guard-rd.efficiency_inside_guard);
end
