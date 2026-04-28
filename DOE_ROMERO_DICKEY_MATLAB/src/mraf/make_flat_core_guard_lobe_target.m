function target = make_flat_core_guard_lobe_target(focal_x_m, focal_y_m, cfg)
% make_flat_core_guard_lobe_target Four-region MRAF target with shoulder guard.
% Regions: flat_core (constrained), shoulder_guard (soft cap),
% lobe_reservoir (free), outer_suppress (weak damping).

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
target_half_x_m = cfg.target_half_x_m;
target_half_y_m = cfg.target_half_y_m;

flat_half_x_m = cfg.mraf.flat_fraction_x * target_half_x_m;
flat_half_y_m = cfg.mraf.flat_fraction_y * target_half_y_m;
guard_half_x_m = cfg.mraf.guard_fraction_x * target_half_x_m;
guard_half_y_m = cfg.mraf.guard_fraction_y * target_half_y_m;
lobe_half_x_m = cfg.mraf.lobe_fraction_x * target_half_x_m;
lobe_half_y_m = cfg.mraf.lobe_fraction_y * target_half_y_m;

flat_core = abs(XF_m) <= flat_half_x_m & abs(YF_m) <= flat_half_y_m;
guard_box = abs(XF_m) <= guard_half_x_m & abs(YF_m) <= guard_half_y_m;
lobe_box = abs(XF_m) <= lobe_half_x_m & abs(YF_m) <= lobe_half_y_m;
shoulder_guard = guard_box & ~flat_core;
lobe_reservoir = lobe_box & ~guard_box;
outer_suppress = ~lobe_box;

A_target = nan(size(XF_m));
A_target(flat_core) = 1;
I_target = A_target.^2;

cap_outer = cfg.mraf.guard_cap_outer_intensity;
cap_intensity = nan(size(XF_m));
cap_intensity(shoulder_guard) = shoulder_guard_cap_map(XF_m(shoulder_guard), YF_m(shoulder_guard), flat_half_x_m, flat_half_y_m, guard_half_x_m, guard_half_y_m, cap_outer);
cap_amplitude = sqrt(max(cap_intensity, 0));

masks = struct();
masks.flat_core = flat_core;
masks.core = flat_core;
masks.signal = flat_core;
masks.finite_target = flat_core;
masks.target = flat_core;
masks.shoulder_guard = shoulder_guard;
masks.free = lobe_reservoir;
masks.free_edge = lobe_reservoir;
masks.free_noise = lobe_reservoir;
masks.lobe_reservoir = lobe_reservoir;
masks.side_lobe_reservoir = lobe_reservoir;
masks.outer_suppress = outer_suppress;
masks.noise = shoulder_guard | lobe_reservoir | outer_suppress;
masks.guard = lobe_box;
masks.edge = shoulder_guard | lobe_reservoir;
masks.core_rect = flat_core;
masks.wgs_update_scale = double(flat_core);
masks.guard_cap_amplitude = cap_amplitude;
masks.guard_cap_intensity = cap_intensity;
masks.wgs_y_edge_damping_band = false(size(flat_core));
masks.y_free_reservoir = false(size(flat_core));
masks.free_threshold = NaN;
masks.edge_low_threshold_effective = NaN;

target = struct();
target.amplitude = A_target;
target.intensity = I_target;
target.intensity_plot = I_target;
target.intensity_plot(~isfinite(target.intensity_plot)) = 0;
target.guard_cap_intensity = cap_intensity;
target.guard_cap_amplitude = cap_amplitude;
target.I_rd_norm = [];
target.I_target_geometry_norm = double(flat_core);
target.I_target_continuous = target.intensity_plot;
target.I_env = target.intensity_plot;
target.blend_map = double(flat_core);
target.free_mask = lobe_reservoir;
target.masks = masks;
target.mode = "flat_core_guard_lobe";
target.edge_mode = "none";
target.legacy_edge_mode = "none";
target.gamma = 1;
target.flat_low = NaN;
target.flat_high = NaN;
target.free_threshold = NaN;
target.pivot50_inner_power = 1;
target.pivot50_outer_power = 1;
target.target_x_scale = 1;
target.target_y_scale = 1;
target.remap_curve_x = linspace(0, 1, 501);
target.remap_curve_y = target.remap_curve_x;
target.remap_anchor_50 = 0.5;
target.core_mean_normalization = NaN;
target.flat_half_x_m = flat_half_x_m;
target.flat_half_y_m = flat_half_y_m;
target.guard_half_x_m = guard_half_x_m;
target.guard_half_y_m = guard_half_y_m;
target.lobe_half_x_m = lobe_half_x_m;
target.lobe_half_y_m = lobe_half_y_m;
target.flat_fraction_x = cfg.mraf.flat_fraction_x;
target.flat_fraction_y = cfg.mraf.flat_fraction_y;
target.guard_fraction_x = cfg.mraf.guard_fraction_x;
target.guard_fraction_y = cfg.mraf.guard_fraction_y;
target.lobe_fraction_x = cfg.mraf.lobe_fraction_x;
target.lobe_fraction_y = cfg.mraf.lobe_fraction_y;
target.guard_cap_outer_intensity = cfg.mraf.guard_cap_outer_intensity;
target.guard_blend = cfg.mraf.guard_blend;
end

function cap = shoulder_guard_cap_map(x_m, y_m, flat_half_x_m, flat_half_y_m, guard_half_x_m, guard_half_y_m, cap_outer)
rx = normalized_outside_fraction(abs(x_m), flat_half_x_m, guard_half_x_m);
ry = normalized_outside_fraction(abs(y_m), flat_half_y_m, guard_half_y_m);
t = max(rx, ry);
t = min(max(t, 0), 1);
smooth_t = t.^2 .* (3 - 2 .* t);
cap = 1 + (cap_outer - 1) .* smooth_t;
end

function r = normalized_outside_fraction(v, inner_half, outer_half)
r = (v - inner_half) ./ max(outer_half - inner_half, eps);
r(v <= inner_half) = 0;
r(v >= outer_half) = 1;
end
