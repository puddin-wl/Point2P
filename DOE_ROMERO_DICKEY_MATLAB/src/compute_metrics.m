function metrics = compute_metrics(focal_x_m, focal_y_m, intensity, target_intensity, cfg)
% compute_metrics 计算焦平面矩形平顶诊断指标。

intensity = intensity ./ max(intensity(:));
center_x_index = round(numel(focal_x_m) / 2) + 1;
center_y_index = round(numel(focal_y_m) / 2) + 1;
x_profile = intensity(center_y_index, :);
y_profile = intensity(:, center_x_index).';
target_x_profile = target_intensity(center_y_index, :);
target_y_profile = target_intensity(:, center_x_index).';

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
target_roi = abs(XF_m) <= cfg.target_half_x_m & abs(YF_m) <= cfg.target_half_y_m;
core_roi = abs(XF_m) <= cfg.core_fraction * cfg.target_size_x_m & abs(YF_m) <= cfg.core_fraction * cfg.target_size_y_m;
guard_roi = abs(XF_m) <= cfg.guard_fraction * cfg.target_half_x_m & abs(YF_m) <= cfg.guard_fraction * cfg.target_half_y_m;

core_values = intensity(core_roi);
core_mean = mean(core_values, 'omitnan');
core_rms = sqrt(mean((core_values ./ max(core_mean, eps) - 1).^2, 'omitnan'));
peak_to_valley = (max(core_values) - min(core_values)) / max(core_mean, eps);

metrics = struct();
metrics.normalization = "focal intensity normalized by global max";
metrics.fwhm_x_m = crossing_width(focal_x_m, x_profile, 0.5);
metrics.fwhm_y_m = crossing_width(focal_y_m, y_profile, 0.5);
metrics.transition_width_13_90_x_m = transition_width(focal_x_m, x_profile, cfg.profile_transition_low, cfg.profile_transition_high, cfg.target_half_x_m);
metrics.transition_width_13_90_y_m = transition_width(focal_y_m, y_profile, cfg.profile_transition_low, cfg.profile_transition_high, cfg.target_half_y_m);
metrics.core_mean = core_mean;
metrics.core_rms = core_rms;
metrics.peak_to_valley = peak_to_valley;
metrics.center_profile_std_x = std(x_profile(abs(focal_x_m) <= cfg.core_fraction * cfg.target_size_x_m), 0, 'omitnan') / max(mean(x_profile(abs(focal_x_m) <= cfg.core_fraction * cfg.target_size_x_m), 'omitnan'), eps);
metrics.center_profile_std_y = std(y_profile(abs(focal_y_m) <= cfg.core_fraction * cfg.target_size_y_m), 0, 'omitnan') / max(mean(y_profile(abs(focal_y_m) <= cfg.core_fraction * cfg.target_size_y_m), 'omitnan'), eps);
metrics.power_inside_target_roi = sum(intensity(target_roi), 'all');
metrics.power_inside_guard_roi = sum(intensity(guard_roi), 'all');
metrics.power_total = sum(intensity, 'all');
metrics.efficiency_target = metrics.power_inside_target_roi / max(metrics.power_total, eps);
metrics.efficiency_guard = metrics.power_inside_guard_roi / max(metrics.power_total, eps);
metrics.side_lobe_peak_x_rel_to_core = side_lobe_peak(focal_x_m, x_profile, cfg.side_lobe_exclusion_fraction * cfg.target_half_x_m, core_mean);
metrics.side_lobe_peak_y_rel_to_core = side_lobe_peak(focal_y_m, y_profile, cfg.side_lobe_exclusion_fraction * cfg.target_half_y_m, core_mean);
metrics.overshoot_near_edge_x = near_edge_extreme(focal_x_m, x_profile, cfg.target_half_x_m, 1, core_mean);
metrics.overshoot_near_edge_y = near_edge_extreme(focal_y_m, y_profile, cfg.target_half_y_m, 1, core_mean);
metrics.undershoot_inside_x = near_edge_extreme(focal_x_m, x_profile, cfg.target_half_x_m, -1, core_mean);
metrics.undershoot_inside_y = near_edge_extreme(focal_y_m, y_profile, cfg.target_half_y_m, -1, core_mean);
metrics.x_profile = x_profile;
metrics.y_profile = y_profile;
metrics.target_x_profile = target_x_profile;
metrics.target_y_profile = target_y_profile;
metrics.crossing_x_50_m = crossing_positions(focal_x_m, x_profile, 0.5);
metrics.crossing_y_50_m = crossing_positions(focal_y_m, y_profile, 0.5);
metrics.crossing_x_13_m = crossing_positions(focal_x_m, x_profile, cfg.profile_transition_low);
metrics.crossing_x_90_m = crossing_positions(focal_x_m, x_profile, cfg.profile_transition_high);
metrics.crossing_y_13_m = crossing_positions(focal_y_m, y_profile, cfg.profile_transition_low);
metrics.crossing_y_90_m = crossing_positions(focal_y_m, y_profile, cfg.profile_transition_high);
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
profile = profile(:).';
axis_m = axis_m(:).';
above = profile >= level;
idx = find(diff(above) ~= 0);
positions_m = zeros(1, numel(idx));
for k = 1:numel(idx)
    i = idx(k);
    x1 = axis_m(i); x2 = axis_m(i + 1);
    y1 = profile(i); y2 = profile(i + 1);
    if abs(y2 - y1) < eps
        positions_m(k) = (x1 + x2) / 2;
    else
        positions_m(k) = x1 + (level - y1) * (x2 - x1) / (y2 - y1);
    end
end
end

function width_m = transition_width(axis_m, profile, low, high, half_target_m)
positive = axis_m >= 0;
axis_pos = axis_m(positive);
profile_pos = profile(positive);
edge_band = axis_pos >= 0.2 * half_target_m & axis_pos <= 3.0 * half_target_m;
axis_pos = axis_pos(edge_band);
profile_pos = profile_pos(edge_band);
if isempty(axis_pos)
    width_m = NaN;
    return;
end
[~, edge_idx] = min(abs(axis_pos - half_target_m));
left_axis = axis_pos(1:edge_idx);
left_profile = profile_pos(1:edge_idx);
right_axis = axis_pos(edge_idx:end);
right_profile = profile_pos(edge_idx:end);
x_high = first_cross_after_reverse(left_axis, left_profile, high);
x_low = first_cross_forward(right_axis, right_profile, low);
if isnan(x_high) || isnan(x_low)
    width_m = NaN;
else
    width_m = abs(x_low - x_high);
end
end

function x_m = first_cross_forward(axis_m, profile, level)
x_m = NaN;
for i = 1:(numel(profile) - 1)
    if (profile(i) - level) * (profile(i + 1) - level) <= 0
        x_m = interp_cross(axis_m(i), axis_m(i + 1), profile(i), profile(i + 1), level);
        return;
    end
end
end

function x_m = first_cross_after_reverse(axis_m, profile, level)
x_m = NaN;
for i = (numel(profile) - 1):-1:1
    if (profile(i) - level) * (profile(i + 1) - level) <= 0
        x_m = interp_cross(axis_m(i), axis_m(i + 1), profile(i), profile(i + 1), level);
        return;
    end
end
end

function x_m = interp_cross(x1, x2, y1, y2, level)
if abs(y2 - y1) < eps
    x_m = (x1 + x2) / 2;
else
    x_m = x1 + (level - y1) * (x2 - x1) / (y2 - y1);
end
end

function value = side_lobe_peak(axis_m, profile, exclusion_half_m, core_mean)
outside = abs(axis_m) > exclusion_half_m;
if any(outside)
    value = max(profile(outside)) / max(core_mean, eps);
else
    value = NaN;
end
end

function value = near_edge_extreme(axis_m, profile, half_target_m, mode, core_mean)
band = abs(abs(axis_m) - half_target_m) <= 0.5 * half_target_m;
if ~any(band)
    value = NaN;
elseif mode > 0
    value = max(profile(band)) / max(core_mean, eps) - 1;
else
    inside = abs(axis_m) <= half_target_m;
    use = band & inside;
    if any(use)
        value = 1 - min(profile(use)) / max(core_mean, eps);
    else
        value = NaN;
    end
end
end

