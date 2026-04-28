function diagnostics = compute_focal_diagnostics(focal_x_m, focal_y_m, intensity, varargin)
% compute_focal_diagnostics Compute common focal-plane result diagnostics.
%
% Inputs:
%   focal_x_m, focal_y_m: focal-plane coordinate vectors in meters.
%   intensity: 2D focal-plane intensity array.
%
% Name/value options:
%   'cfg': diagnostics config, default_diagnostics_config() if omitted.
%   'normalization': 'core_mean' (default) or 'max'.
%   'target_intensity': optional target map for profile plotting.
%   'masks': optional struct with fields core/signal/noise/guard/etc.
%
% Output:
%   diagnostics struct containing center profiles, size crossings,
%   transition width, core RMS/PV, shoulder, side-lobe, and energy ratios.

opts = parse_options(varargin{:});
cfg = opts.cfg;
intensity = double(intensity);
intensity(~isfinite(intensity)) = 0;

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
masks = normalize_masks(opts.masks, XF_m, YF_m, cfg);

center_x_index = round(numel(focal_x_m) / 2) + 1;
center_y_index = round(numel(focal_y_m) / 2) + 1;

raw_core = intensity(masks.core);
core_mean_raw = mean(raw_core, 'omitnan');
if ~(isfinite(core_mean_raw) && core_mean_raw > 0)
    core_mean_raw = max(intensity(:));
end

switch lower(string(opts.normalization))
    case "max"
        norm_den = max(intensity(:));
        normalization_text = "global max";
    otherwise
        norm_den = core_mean_raw;
        normalization_text = "core mean";
end
I = intensity ./ max(norm_den, eps);

x_profile = I(center_y_index, :);
y_profile = I(:, center_x_index).';

x90 = crossing_pair(focal_x_m, x_profile, 0.90);
y90 = crossing_pair(focal_y_m, y_profile, 0.90);
x50 = crossing_pair(focal_x_m, x_profile, 0.50);
y50 = crossing_pair(focal_y_m, y_profile, 0.50);
x135 = crossing_pair(focal_x_m, x_profile, 0.135);
y135 = crossing_pair(focal_y_m, y_profile, 0.135);
[tw_x, tw_x_left, tw_x_right] = transition_from_pairs(x135, x90);
[tw_y, tw_y_bottom, tw_y_top] = transition_from_pairs(y135, y90);

core_values = I(masks.core);
core_mean = mean(core_values, 'omitnan');
core_norm = core_values ./ max(core_mean, eps);
signal_power = sum(intensity(masks.signal), 'all');
total_power = sum(intensity, 'all');

half50_x_m = pair_width_um(x50) * 1e-6 / 2;
half50_y_m = pair_width_um(y50) * 1e-6 / 2;

diagnostics = struct();
diagnostics.normalization = normalization_text;
diagnostics.core_mean_raw = core_mean_raw;
diagnostics.core_rms = sqrt(mean((core_norm - 1).^2, 'omitnan'));
diagnostics.core_pv = max(core_norm) - min(core_norm);
diagnostics.core_uniformity = 1 - diagnostics.core_pv;
diagnostics.size90_x_um = pair_width_um(x90);
diagnostics.size90_y_um = pair_width_um(y90);
diagnostics.size50_x_um = pair_width_um(x50);
diagnostics.size50_y_um = pair_width_um(y50);
diagnostics.size13p5_x_um = pair_width_um(x135);
diagnostics.size13p5_y_um = pair_width_um(y135);
diagnostics.transition_13p5_90_x_um = tw_x;
diagnostics.transition_13p5_90_y_um = tw_y;
diagnostics.transition_13p5_90_x_left_um = tw_x_left;
diagnostics.transition_13p5_90_x_right_um = tw_x_right;
diagnostics.transition_13p5_90_y_bottom_um = tw_y_bottom;
diagnostics.transition_13p5_90_y_top_um = tw_y_top;
diagnostics.crossing90_x_left_um = x90.left * 1e6;
diagnostics.crossing90_x_right_um = x90.right * 1e6;
diagnostics.crossing90_y_bottom_um = y90.left * 1e6;
diagnostics.crossing90_y_top_um = y90.right * 1e6;
diagnostics.crossing50_x_left_um = x50.left * 1e6;
diagnostics.crossing50_x_right_um = x50.right * 1e6;
diagnostics.crossing50_y_bottom_um = y50.left * 1e6;
diagnostics.crossing50_y_top_um = y50.right * 1e6;
diagnostics.crossing13p5_x_left_um = x135.left * 1e6;
diagnostics.crossing13p5_x_right_um = x135.right * 1e6;
diagnostics.crossing13p5_y_bottom_um = y135.left * 1e6;
diagnostics.crossing13p5_y_top_um = y135.right * 1e6;
[diagnostics.shoulder_peak_x, diagnostics.shoulder_peak_left, diagnostics.shoulder_peak_right] = shoulder_peak(focal_x_m, x_profile, half50_x_m, cfg);
[diagnostics.shoulder_peak_y, diagnostics.shoulder_peak_bottom, diagnostics.shoulder_peak_top] = shoulder_peak(focal_y_m, y_profile, half50_y_m, cfg);
diagnostics.side_lobe_peak_x_rel_to_core = outer_tail_peak(focal_x_m, x_profile, cfg.side_lobe_exclusion_fraction * cfg.target_half_x_m, core_mean);
diagnostics.side_lobe_peak_y_rel_to_core = outer_tail_peak(focal_y_m, y_profile, cfg.side_lobe_exclusion_fraction * cfg.target_half_y_m, core_mean);
[diagnostics.true_side_lobe_peak_x_rel_to_core, diagnostics.true_side_lobe_pos_x_um, diagnostics.has_true_side_lobe_x] = true_side_lobe_peak(focal_x_m, x_profile, x135, core_mean, cfg);
[diagnostics.true_side_lobe_peak_y_rel_to_core, diagnostics.true_side_lobe_pos_y_um, diagnostics.has_true_side_lobe_y] = true_side_lobe_peak(focal_y_m, y_profile, y135, core_mean, cfg);
diagnostics.energy_inside_signal = signal_power / max(total_power, eps);
diagnostics.energy_inside_guard = sum(intensity(masks.guard), 'all') / max(total_power, eps);
diagnostics.noise_energy_ratio = sum(intensity(masks.noise), 'all') / max(signal_power, eps);
diagnostics.null_energy_ratio = sum(intensity(masks.null), 'all') / max(signal_power, eps);
diagnostics.total_power = total_power;
diagnostics.x_profile = x_profile;
diagnostics.y_profile = y_profile;
diagnostics.focal_x_m = focal_x_m;
diagnostics.focal_y_m = focal_y_m;
diagnostics.masks = masks;
diagnostics.target_x_profile = target_profile(opts.target_intensity, center_y_index, 2);
diagnostics.target_y_profile = target_profile(opts.target_intensity, center_x_index, 1);
end

function opts = parse_options(varargin)
opts = struct('cfg', default_diagnostics_config(), 'normalization', "core_mean", 'target_intensity', [], 'masks', struct());
if mod(numel(varargin), 2) ~= 0, error('Options must be name/value pairs.'); end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    opts.(char(name)) = varargin{k + 1};
end
end

function masks = normalize_masks(input_masks, XF_m, YF_m, cfg)
base_core = abs(XF_m) <= cfg.core_fraction * cfg.target_size_x_m & abs(YF_m) <= cfg.core_fraction * cfg.target_size_y_m;
base_signal = abs(XF_m) <= cfg.target_half_x_m & abs(YF_m) <= cfg.target_half_y_m;
base_guard = abs(XF_m) <= cfg.guard_fraction * cfg.target_half_x_m & abs(YF_m) <= cfg.guard_fraction * cfg.target_half_y_m;
masks = input_masks;
if ~isfield(masks, 'core'), masks.core = base_core; end
if ~isfield(masks, 'signal'), masks.signal = base_signal; end
if ~isfield(masks, 'guard'), masks.guard = base_guard; end
if ~isfield(masks, 'noise'), masks.noise = masks.guard & ~masks.signal; end
if ~isfield(masks, 'null'), masks.null = ~masks.guard; end
end

function profile = target_profile(target_intensity, index, dim)
if isempty(target_intensity)
    profile = [];
elseif dim == 2
    profile = target_intensity(index, :);
else
    profile = target_intensity(:, index).';
end
end

function pair = crossing_pair(axis_m, profile, level)
positions = crossing_positions(axis_m, profile, level);
neg = positions(positions <= 0);
pos = positions(positions >= 0);
pair = struct('left', NaN, 'right', NaN, 'level', level);
if ~isempty(neg), pair.left = max(neg); end
if ~isempty(pos), pair.right = min(pos); end
end

function positions_m = crossing_positions(axis_m, profile, level)
profile = profile(:).'; axis_m = axis_m(:).';
above = profile >= level;
idx = find(diff(above) ~= 0);
positions_m = zeros(1, numel(idx));
for k = 1:numel(idx)
    i = idx(k);
    positions_m(k) = interp_cross(axis_m(i), axis_m(i + 1), profile(i), profile(i + 1), level);
end
end

function x_m = interp_cross(x1, x2, y1, y2, level)
if abs(y2 - y1) < eps
    x_m = (x1 + x2) / 2;
else
    x_m = x1 + (level - y1) * (x2 - x1) / (y2 - y1);
end
end

function width_um = pair_width_um(pair)
if isfinite(pair.left) && isfinite(pair.right)
    width_um = (pair.right - pair.left) * 1e6;
else
    width_um = NaN;
end
end

function [tw_um, left_um, right_um] = transition_from_pairs(low_pair, high_pair)
if isfinite(low_pair.left) && isfinite(high_pair.left)
    left_um = (high_pair.left - low_pair.left) * 1e6;
else
    left_um = NaN;
end
if isfinite(low_pair.right) && isfinite(high_pair.right)
    right_um = (low_pair.right - high_pair.right) * 1e6;
else
    right_um = NaN;
end
tw_um = mean([left_um, right_um], 'omitnan');
if isnan(left_um) && isnan(right_um), tw_um = NaN; end
end

function [peak_both, peak_neg, peak_pos] = shoulder_peak(axis_m, profile, half50_m, cfg)
if ~(isfinite(half50_m) && half50_m > 0)
    peak_both = NaN; peak_neg = NaN; peak_pos = NaN; return;
end
inner = cfg.shoulder_inner_fraction;
outer = cfg.shoulder_outer_fraction;
neg_band = axis_m <= -inner * half50_m & axis_m >= -outer * half50_m;
pos_band = axis_m >= inner * half50_m & axis_m <= outer * half50_m;
if any(neg_band), peak_neg = max(profile(neg_band)) - 1; else, peak_neg = NaN; end
if any(pos_band), peak_pos = max(profile(pos_band)) - 1; else, peak_pos = NaN; end
peak_both = max([peak_neg, peak_pos], [], 'omitnan');
end

function value = outer_tail_peak(axis_m, profile, exclusion_half_m, core_mean)
outside = abs(axis_m) > exclusion_half_m;
if any(outside), value = max(profile(outside)) / max(core_mean, eps); else, value = NaN; end
end

function [peak_value, peak_pos_um, has_peak] = true_side_lobe_peak(axis_m, profile, low_pair, core_mean, cfg)
profile = profile(:).'; axis_m = axis_m(:).';
smooth_profile = moving_average_1d(profile, cfg.true_side_lobe_smooth_window);
margin_m = 2.5e-6;
best_peak = 0;
best_pos = NaN;
if isfinite(low_pair.right)
    [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, smooth_profile, low_pair.right + margin_m, 1, cfg);
    if peak_height > best_peak, best_peak = peak_height; best_pos = peak_pos; end
end
if isfinite(low_pair.left)
    [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, smooth_profile, low_pair.left - margin_m, -1, cfg);
    if peak_height > best_peak, best_peak = peak_height; best_pos = peak_pos; end
end
peak_value = best_peak / max(core_mean, eps);
peak_pos_um = best_pos * 1e6;
has_peak = double(isfinite(best_pos) && peak_value > 0);
end

function [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, profile, start_m, side_sign, cfg)
if side_sign > 0
    idx = find(axis_m >= start_m);
else
    idx = find(axis_m <= start_m);
    idx = fliplr(idx);
end
peak_pos = NaN; peak_height = 0;
if numel(idx) < 5, return; end
axis_side = axis_m(idx); profile_side = profile(idx);
d = diff(profile_side);
for k = 2:numel(profile_side)-1
    if d(k-1) > 0 && d(k) <= 0 && profile_side(k) >= cfg.true_side_lobe_min_height
        left_min = min(profile_side(1:k));
        right_stop = min(numel(profile_side), k + 20);
        right_min = min(profile_side(k:right_stop));
        prominence = profile_side(k) - max(left_min, right_min);
        if prominence >= cfg.true_side_lobe_min_prominence && profile_side(k) > peak_height
            peak_height = profile_side(k);
            peak_pos = axis_side(k);
        end
    end
end
end

function y = moving_average_1d(x, window)
if numel(x) < window
    y = x;
else
    kernel = ones(1, window) ./ window;
    y = conv(x, kernel, 'same');
end
end
