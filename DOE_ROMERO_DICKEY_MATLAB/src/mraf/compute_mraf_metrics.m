function metrics = compute_mraf_metrics(focal_x_m, focal_y_m, intensity_raw, masks, input_power, cfg)
% compute_mraf_metrics Metrics normalized by core mean for RD/MRAF comparison.
intensity_raw = double(intensity_raw);
intensity_raw(~isfinite(intensity_raw)) = 0;
total_power_raw = sum(intensity_raw, 'all');
core_values_raw = intensity_raw(masks.core);
core_mean_raw = mean(core_values_raw, 'omitnan');
if ~(isfinite(core_mean_raw) && core_mean_raw > 0)
    core_mean_raw = max(intensity_raw(:));
end
I = intensity_raw ./ max(core_mean_raw, eps);
center_x_index = round(numel(focal_x_m) / 2) + 1;
center_y_index = round(numel(focal_y_m) / 2) + 1;
x_profile = I(center_y_index, :);
y_profile = I(:, center_x_index).';
core_values = I(masks.core);
core_mean = mean(core_values, 'omitnan');

x50 = crossing_pair(focal_x_m, x_profile, 0.50);
y50 = crossing_pair(focal_y_m, y_profile, 0.50);
x13p5 = crossing_pair(focal_x_m, x_profile, 0.135);
y13p5 = crossing_pair(focal_y_m, y_profile, 0.135);
x90 = crossing_pair(focal_x_m, x_profile, 0.90);
y90 = crossing_pair(focal_y_m, y_profile, 0.90);
[x_tw, x_tw_left, x_tw_right] = transition_from_pairs(x13p5, x90);
[y_tw, y_tw_left, y_tw_right] = transition_from_pairs(y13p5, y90);

metrics = struct();
metrics.normalization = "focal intensity normalized by RD-derived core mean";
metrics.size50_x_um = pair_width_um(x50);
metrics.size50_y_um = pair_width_um(y50);
metrics.size_50_x_um = metrics.size50_x_um;
metrics.size_50_y_um = metrics.size50_y_um;
metrics.size13p5_x_um = pair_width_um(x13p5);
metrics.size13p5_y_um = pair_width_um(y13p5);
metrics.transition_13p5_90_x_um = x_tw;
metrics.transition_13p5_90_y_um = y_tw;
metrics.transition_13p5_90_x_left_um = x_tw_left;
metrics.transition_13p5_90_x_right_um = x_tw_right;
metrics.transition_13p5_90_y_left_um = y_tw_left;
metrics.transition_13p5_90_y_right_um = y_tw_right;
metrics.transition_width_13_90_x_um = transition_width(focal_x_m, x_profile, 0.13, 0.90, cfg.target_half_x_m) * 1e6;
metrics.transition_width_13_90_y_um = transition_width(focal_y_m, y_profile, 0.13, 0.90, cfg.target_half_y_m) * 1e6;
metrics.transition_width_10_90_x_um = transition_width(focal_x_m, x_profile, 0.10, 0.90, cfg.target_half_x_m) * 1e6;
metrics.transition_width_10_90_y_um = transition_width(focal_y_m, y_profile, 0.10, 0.90, cfg.target_half_y_m) * 1e6;
metrics.crossing50_x_left_um = x50.left * 1e6;
metrics.crossing50_x_right_um = x50.right * 1e6;
metrics.crossing50_y_left_um = y50.left * 1e6;
metrics.crossing50_y_right_um = y50.right * 1e6;
metrics.crossing13p5_x_left_um = x13p5.left * 1e6;
metrics.crossing13p5_x_right_um = x13p5.right * 1e6;
metrics.crossing13p5_y_left_um = y13p5.left * 1e6;
metrics.crossing13p5_y_right_um = y13p5.right * 1e6;
metrics.crossing90_x_left_um = x90.left * 1e6;
metrics.crossing90_x_right_um = x90.right * 1e6;
metrics.crossing90_y_left_um = y90.left * 1e6;
metrics.crossing90_y_right_um = y90.right * 1e6;
metrics.has_crossing_warning = any(isnan([metrics.size50_x_um, metrics.size50_y_um, metrics.size13p5_x_um, metrics.size13p5_y_um, metrics.transition_13p5_90_x_um, metrics.transition_13p5_90_y_um]));
metrics.core_rms = sqrt(mean((core_values ./ max(core_mean, eps) - 1).^2, 'omitnan'));
roi90 = masks.core & I >= 0.90;
if nnz(roi90) < 4
    roi90 = masks.core;
end
values90 = I(roi90);
mean90 = mean(values90, 'omitnan');
metrics.rms_90 = sqrt(mean((values90 ./ max(mean90, eps) - 1).^2, 'omitnan'));
metrics.peak_to_valley = (max(core_values) - min(core_values)) / max(core_mean, eps);
half50_x_m = metrics.size50_x_um * 1e-6 / 2;
half50_y_m = metrics.size50_y_um * 1e-6 / 2;
[metrics.shoulder_peak_x, metrics.shoulder_peak_left, metrics.shoulder_peak_right] = shoulder_peak(focal_x_m, x_profile, half50_x_m);
[metrics.shoulder_peak_y, metrics.shoulder_peak_bottom, metrics.shoulder_peak_top] = shoulder_peak(focal_y_m, y_profile, half50_y_m);
metrics.overshoot_x = near_edge_extreme(focal_x_m, x_profile, cfg.target_half_x_m, 1, core_mean);
metrics.overshoot_y = near_edge_extreme(focal_y_m, y_profile, cfg.target_half_y_m, 1, core_mean);
metrics.undershoot_x = near_edge_extreme(focal_x_m, x_profile, cfg.target_half_x_m, -1, core_mean);
metrics.undershoot_y = near_edge_extreme(focal_y_m, y_profile, cfg.target_half_y_m, -1, core_mean);
tail_exclusion_x_m = cfg.side_lobe_exclusion_fraction * cfg.target_half_x_m;
tail_exclusion_y_m = cfg.side_lobe_exclusion_fraction * cfg.target_half_y_m;
metrics.outer_tail_peak_x_rel_to_core = outer_tail_peak(focal_x_m, x_profile, tail_exclusion_x_m, core_mean);
metrics.outer_tail_peak_y_rel_to_core = outer_tail_peak(focal_y_m, y_profile, tail_exclusion_y_m, core_mean);
metrics.side_lobe_peak_x_rel_to_core = metrics.outer_tail_peak_x_rel_to_core;
metrics.side_lobe_peak_y_rel_to_core = metrics.outer_tail_peak_y_rel_to_core;
[metrics.true_side_lobe_peak_x_rel_to_core, metrics.true_side_lobe_pos_x_um, metrics.has_true_side_lobe_x] = true_side_lobe_peak(focal_x_m, x_profile, x13p5, core_mean);
[metrics.true_side_lobe_peak_y_rel_to_core, metrics.true_side_lobe_pos_y_um, metrics.has_true_side_lobe_y] = true_side_lobe_peak(focal_y_m, y_profile, y13p5, core_mean);
metrics.efficiency_inside_signal = sum(intensity_raw(masks.signal), 'all') / max(total_power_raw, eps);
metrics.efficiency_inside_guard = sum(intensity_raw(masks.guard), 'all') / max(total_power_raw, eps);
metrics.total_power = total_power_raw;
metrics.input_power = input_power;
metrics.total_power_conservation_check = total_power_raw / max(input_power, eps);
metrics.x_profile = x_profile;
metrics.y_profile = y_profile;
end

function pair = crossing_pair(axis_m, profile, level)
positions = crossing_positions(axis_m, profile, level);
neg = positions(positions <= 0);
pos = positions(positions >= 0);
pair = struct('left', NaN, 'right', NaN, 'level', level);
if ~isempty(neg), pair.left = max(neg); end
if ~isempty(pos), pair.right = min(pos); end
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

function [peak_both, peak_neg, peak_pos] = shoulder_peak(axis_m, profile, half50_m)
if ~(isfinite(half50_m) && half50_m > 0)
    peak_both = NaN; peak_neg = NaN; peak_pos = NaN; return;
end
neg_band = axis_m <= -0.70 * half50_m & axis_m >= -1.05 * half50_m;
pos_band = axis_m >= 0.70 * half50_m & axis_m <= 1.05 * half50_m;
if any(neg_band), peak_neg = max(profile(neg_band)) - 1; else, peak_neg = NaN; end
if any(pos_band), peak_pos = max(profile(pos_band)) - 1; else, peak_pos = NaN; end
peak_both = max([peak_neg, peak_pos], [], 'omitnan');
end

function width_m = crossing_width(axis_m, profile, level)
positions = crossing_positions(axis_m, profile, level);
if numel(positions) >= 2
    width_m = max(positions) - min(positions);
else
    width_m = NaN;
end
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

function width_m = transition_width(axis_m, profile, low, high, half_target_m)
positive = axis_m >= 0;
axis_pos = axis_m(positive); profile_pos = profile(positive);
edge_band = axis_pos >= 0.2 * half_target_m & axis_pos <= 3.0 * half_target_m;
axis_pos = axis_pos(edge_band); profile_pos = profile_pos(edge_band);
if isempty(axis_pos), width_m = NaN; return; end
[~, edge_idx] = min(abs(axis_pos - half_target_m));
x_high = first_cross_after_reverse(axis_pos(1:edge_idx), profile_pos(1:edge_idx), high);
x_low = first_cross_forward(axis_pos(edge_idx:end), profile_pos(edge_idx:end), low);
if isnan(x_high) || isnan(x_low), width_m = NaN; else, width_m = abs(x_low - x_high); end
end

function x_m = first_cross_forward(axis_m, profile, level)
x_m = NaN;
for i = 1:(numel(profile) - 1)
    if (profile(i) - level) * (profile(i + 1) - level) <= 0
        x_m = interp_cross(axis_m(i), axis_m(i + 1), profile(i), profile(i + 1), level); return;
    end
end
end

function x_m = first_cross_after_reverse(axis_m, profile, level)
x_m = NaN;
for i = (numel(profile) - 1):-1:1
    if (profile(i) - level) * (profile(i + 1) - level) <= 0
        x_m = interp_cross(axis_m(i), axis_m(i + 1), profile(i), profile(i + 1), level); return;
    end
end
end

function x_m = interp_cross(x1, x2, y1, y2, level)
if abs(y2 - y1) < eps, x_m = (x1 + x2) / 2; else, x_m = x1 + (level - y1) * (x2 - x1) / (y2 - y1); end
end

function value = outer_tail_peak(axis_m, profile, exclusion_half_m, core_mean)
outside = abs(axis_m) > exclusion_half_m;
if any(outside), value = max(profile(outside)) / max(core_mean, eps); else, value = NaN; end
end

function [peak_value, peak_pos_um, has_peak] = true_side_lobe_peak(axis_m, profile, low_pair, core_mean)
profile = profile(:).'; axis_m = axis_m(:).';
smooth_profile = moving_average_1d(profile, 5);
min_height = 0.02;
min_prominence = 0.005;
margin_m = 2.5e-6;
best_peak = 0;
best_pos = NaN;
if isfinite(low_pair.right)
    [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, smooth_profile, low_pair.right + margin_m, 1, min_height, min_prominence);
    if peak_height > best_peak
        best_peak = peak_height; best_pos = peak_pos;
    end
end
if isfinite(low_pair.left)
    [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, smooth_profile, low_pair.left - margin_m, -1, min_height, min_prominence);
    if peak_height > best_peak
        best_peak = peak_height; best_pos = peak_pos;
    end
end
peak_value = best_peak / max(core_mean, eps);
peak_pos_um = best_pos * 1e6;
has_peak = double(isfinite(best_pos) && peak_value > 0);
end

function [peak_pos, peak_height] = true_side_lobe_one_side(axis_m, profile, start_m, side_sign, min_height, min_prominence)
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
    if d(k-1) > 0 && d(k) <= 0 && profile_side(k) >= min_height
        left_min = min(profile_side(1:k));
        right_stop = min(numel(profile_side), k + 20);
        right_min = min(profile_side(k:right_stop));
        prominence = profile_side(k) - max(left_min, right_min);
        if prominence >= min_prominence && profile_side(k) > peak_height
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

function value = near_edge_extreme(axis_m, profile, half_target_m, mode, core_mean)
band = abs(abs(axis_m) - half_target_m) <= 0.5 * half_target_m;
if ~any(band), value = NaN; return; end
if mode > 0
    value = max(profile(band)) / max(core_mean, eps) - 1;
else
    inside = abs(axis_m) <= half_target_m; use = band & inside;
    if any(use), value = 1 - min(profile(use)) / max(core_mean, eps); else, value = NaN; end
end
end
