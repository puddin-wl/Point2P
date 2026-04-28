% LEGACY:
% This script/function is kept for comparison only.
% The recommended main workflow is run_flat_core_free_edge_mraf.m.

function target = make_rd_derived_target(focal_x_m, focal_y_m, I_rd, cfg)
% make_rd_derived_target Create a smooth RD-derived amplitude target.

if strcmpi(string(cfg.mraf.target_mode), "hard_rectangle")
    target = make_hard_rectangle_target(focal_x_m, focal_y_m, I_rd, cfg);
    return;
end

I_rd = double(I_rd);
I_rd(~isfinite(I_rd)) = 0;
[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
core_rect = abs(XF_m) <= cfg.core_fraction * cfg.target_size_x_m & ...
    abs(YF_m) <= cfg.core_fraction * cfg.target_size_y_m;
core_mean = mean(I_rd(core_rect), 'omitnan');
if ~(isfinite(core_mean) && core_mean > 0)
    core_mean = max(I_rd(:));
end
I_norm = I_rd ./ max(core_mean, eps);
I_norm = min(max(I_norm, 0), 4);
target_x_scale = get_mraf_scale(cfg, 'target_x_scale');
target_y_scale = get_mraf_scale(cfg, 'target_y_scale');
I_geom = scaled_intensity_geometry(focal_x_m, focal_y_m, I_norm, target_x_scale, target_y_scale);
I_clip = min(max(I_geom, 0), 1);

masks = make_mraf_masks(focal_x_m, focal_y_m, I_geom, cfg);
flat_low = cfg.mraf.flat_low;
flat_high = cfg.mraf.flat_high;
blend_map = zeros(size(I_geom));
blend_map(I_geom >= flat_high) = 1;
blend_region = I_geom > flat_low & I_geom < flat_high;
t = (I_geom(blend_region) - flat_low) ./ max(flat_high - flat_low, eps);
blend_map(blend_region) = t.^2 .* (3 - 2 .* t);
blend_map(~masks.guard) = 0;

target_edge_mode = "identity";
if isfield(cfg.mraf, 'target_edge_mode')
    target_edge_mode = string(cfg.mraf.target_edge_mode);
end
if strcmpi(target_edge_mode, "pivot50")
    I_env = pivot50_remap(I_clip, cfg.mraf.pivot50_inner_power, cfg.mraf.pivot50_outer_power);
elseif strcmpi(target_edge_mode, "identity")
    I_env = I_clip;
elseif strcmpi(string(cfg.mraf.rd_target_edge_mode), "sigmoid")
    k = cfg.mraf.rd_target_sigmoid_k;
    t0 = cfg.mraf.rd_target_sigmoid_t;
    I_env = 1 ./ (1 + exp(-k .* (I_geom - t0)));
else
    I_env = I_clip .^ cfg.mraf.rd_target_gamma;
end
I_flat = ones(size(I_geom));
I_target_continuous = blend_map .* I_flat + (1 - blend_map) .* I_env;
I_target_continuous = min(max(I_target_continuous, 0), 1);

I_target = nan(size(I_geom));
I_target(masks.finite_target) = I_target_continuous(masks.finite_target);
I_target(masks.free) = NaN;
A_target = sqrt(max(I_target, 0));
A_target(~isfinite(I_target)) = NaN;

curve_x = linspace(0, 1, 501);
if strcmpi(target_edge_mode, "pivot50")
    curve_y = pivot50_remap(curve_x, cfg.mraf.pivot50_inner_power, cfg.mraf.pivot50_outer_power);
else
    curve_y = curve_x;
end

target = struct();
target.amplitude = A_target;
target.intensity = I_target;
target.intensity_plot = I_target;
target.intensity_plot(~isfinite(target.intensity_plot)) = 0;
target.I_rd_norm = I_norm;
target.I_target_geometry_norm = I_geom;
target.I_target_continuous = I_target_continuous;
target.I_env = I_env;
target.blend_map = blend_map;
target.free_mask = masks.free;
target.masks = masks;
target.mode = cfg.mraf.target_mode;
target.edge_mode = target_edge_mode;
target.legacy_edge_mode = cfg.mraf.rd_target_edge_mode;
target.gamma = cfg.mraf.rd_target_gamma;
target.flat_low = flat_low;
target.flat_high = flat_high;
target.free_threshold = cfg.mraf.free_threshold;
target.pivot50_inner_power = cfg.mraf.pivot50_inner_power;
target.pivot50_outer_power = cfg.mraf.pivot50_outer_power;
target.target_x_scale = target_x_scale;
target.target_y_scale = target_y_scale;
target.remap_curve_x = curve_x;
target.remap_curve_y = curve_y;
target.remap_anchor_50 = pivot50_remap(0.5, cfg.mraf.pivot50_inner_power, cfg.mraf.pivot50_outer_power);
target.core_mean_normalization = core_mean;
end

function scale = get_mraf_scale(cfg, field_name)
scale = 1;
if isfield(cfg.mraf, field_name)
    scale = double(cfg.mraf.(field_name));
end
if ~(isfinite(scale) && scale > 0)
    scale = 1;
end
end

function I_scaled = scaled_intensity_geometry(focal_x_m, focal_y_m, I_norm, x_scale, y_scale)
if abs(x_scale - 1) < 1e-12 && abs(y_scale - 1) < 1e-12
    I_scaled = I_norm;
    return;
end
[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
I_scaled = interp2(XF_m, YF_m, I_norm, XF_m ./ x_scale, YF_m ./ y_scale, 'linear', 0);
I_scaled(~isfinite(I_scaled)) = 0;
I_scaled = min(max(I_scaled, 0), 4);
end

function I_remap = pivot50_remap(I, inner_power, outer_power)
I = min(max(I, 0), 1);
I_remap = zeros(size(I));
outer = I < 0.5;
I_remap(outer) = 0.5 .* (2 .* I(outer)) .^ outer_power;
inner = ~outer;
I_remap(inner) = 0.5 + 0.5 .* (2 .* I(inner) - 1) .^ inner_power;
I_remap = min(max(I_remap, 0), 1);
end

function target = make_hard_rectangle_target(focal_x_m, focal_y_m, I_rd, cfg)
I_rd = double(I_rd);
I_rd(~isfinite(I_rd)) = 0;
[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
signal_mask = abs(XF_m) <= cfg.target_half_x_m & abs(YF_m) <= cfg.target_half_y_m;
free_mask = ~signal_mask;
I_target = nan(size(I_rd));
I_target(signal_mask) = 1;
A_target = sqrt(max(I_target, 0));
A_target(free_mask) = NaN;

masks = struct();
masks.signal = signal_mask;
masks.core = signal_mask;
masks.edge = false(size(signal_mask));
masks.free = free_mask;
masks.noise = free_mask;
masks.guard = true(size(signal_mask));
masks.finite_target = signal_mask;
masks.core_rect = signal_mask;
masks.free_threshold = NaN;
masks.edge_low_threshold_effective = NaN;

target = struct();
target.amplitude = A_target;
target.intensity = I_target;
target.intensity_plot = I_target;
target.intensity_plot(~isfinite(target.intensity_plot)) = 0;
target.I_rd_norm = I_rd ./ max(mean(I_rd(signal_mask), 'omitnan'), eps);
target.I_target_continuous = target.intensity_plot;
target.I_env = target.intensity_plot;
target.blend_map = double(signal_mask);
target.free_mask = free_mask;
target.masks = masks;
target.mode = "hard_rectangle";
target.edge_mode = "hard_rectangle";
target.legacy_edge_mode = "none";
target.gamma = 1;
target.flat_low = NaN;
target.flat_high = NaN;
target.free_threshold = NaN;
target.pivot50_inner_power = 1;
target.pivot50_outer_power = 1;
target.remap_curve_x = linspace(0, 1, 501);
target.remap_curve_y = target.remap_curve_x;
target.remap_anchor_50 = 0.5;
target.core_mean_normalization = mean(I_rd(signal_mask), 'omitnan');
target.hard_rectangle_half_x_m = cfg.target_half_x_m;
target.hard_rectangle_half_y_m = cfg.target_half_y_m;
target.hard_rectangle_size_x_um = cfg.target_size_x_m * 1e6;
target.hard_rectangle_size_y_um = cfg.target_size_y_m * 1e6;
end

