function plot_center_profiles_diagnostics(diagnostics, output_path, varargin)
% plot_center_profiles_diagnostics Save x/y center profile diagnostic plot.

opts = parse_options(varargin{:});
cfg = opts.cfg;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 760]);
tiledlayout(2, 1, 'TileSpacing', 'compact');
plot_one(nexttile, diagnostics.focal_x_m, diagnostics.x_profile, diagnostics.target_x_profile, diagnostics, 'x', cfg);
plot_one(nexttile, diagnostics.focal_y_m, diagnostics.y_profile, diagnostics.target_y_profile, diagnostics, 'y', cfg);
exportgraphics(fig, output_path, 'Resolution', cfg.figure_dpi);
close(fig);
end

function opts = parse_options(varargin)
opts = struct('cfg', default_diagnostics_config());
if mod(numel(varargin), 2) ~= 0, error('Options must be name/value pairs.'); end
for k = 1:2:numel(varargin)
    opts.(char(lower(string(varargin{k})))) = varargin{k + 1};
end
end

function plot_one(ax, axis_m, profile, target_profile, diagnostics, direction, cfg)
plot(ax, axis_m * 1e6, profile, 'b-', 'LineWidth', 1.2); hold(ax, 'on'); grid(ax, 'on');
if ~isempty(target_profile)
    plot(ax, axis_m * 1e6, target_profile, 'k--', 'LineWidth', 1.0);
end
yline(ax, 0.90, 'g:'); yline(ax, 0.50, 'r:'); yline(ax, 0.135, 'm:');
if direction == 'x'
    half_target_m = cfg.target_half_x_m;
    crossings = [diagnostics.crossing90_x_left_um diagnostics.crossing90_x_right_um diagnostics.crossing50_x_left_um diagnostics.crossing50_x_right_um diagnostics.crossing13p5_x_left_um diagnostics.crossing13p5_x_right_um];
    title_text = sprintf('X center profile | size50 %.2f um | shoulder %.4f | true side-lobe %.4f', diagnostics.size50_x_um, diagnostics.shoulder_peak_x, diagnostics.true_side_lobe_peak_x_rel_to_core);
else
    half_target_m = cfg.target_half_y_m;
    crossings = [diagnostics.crossing90_y_bottom_um diagnostics.crossing90_y_top_um diagnostics.crossing50_y_bottom_um diagnostics.crossing50_y_top_um diagnostics.crossing13p5_y_bottom_um diagnostics.crossing13p5_y_top_um];
    title_text = sprintf('Y center profile | size50 %.2f um | shoulder %.4f | true side-lobe %.4f', diagnostics.size50_y_um, diagnostics.shoulder_peak_y, diagnostics.true_side_lobe_peak_y_rel_to_core);
end
xline(ax, -half_target_m * 1e6, 'k:'); xline(ax, half_target_m * 1e6, 'k:');
colors = {'g','g','r','r','m','m'};
for k = 1:numel(crossings)
    if isfinite(crossings(k)), xline(ax, crossings(k), '-', 'Color', colors{k}, 'Alpha', 0.35); end
end
xlim(ax, [-cfg.center_profile_xlim_factor * half_target_m, cfg.center_profile_xlim_factor * half_target_m] * 1e6);
ylim(ax, [0, max(1.2, min(2.0, max(profile) * 1.05))]);
xlabel(ax, sprintf('%s focal coordinate (um)', direction)); ylabel(ax, 'normalized intensity'); title(ax, title_text, 'Interpreter', 'none');
if isempty(target_profile)
    legend(ax, {'result','90%','50%','13.5%'}, 'Location', 'northeastoutside');
else
    legend(ax, {'result','target','90%','50%','13.5%'}, 'Location', 'northeastoutside');
end
end
