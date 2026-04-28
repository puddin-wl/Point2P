function target = make_flat_core_free_edge_target(focal_x_m, focal_y_m, cfg)
% make_flat_core_free_edge_target Build the recommended simple MRAF target.
% The only constrained far-field region is the central flat_core. The
% surrounding free_edge is an unconstrained side-lobe/transition reservoir,
% and the far outer region is weakly suppressed by project_farfield_mraf.

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);

target_half_x_m = cfg.target_half_x_m;
target_half_y_m = cfg.target_half_y_m;

flat_half_x_m = cfg.mraf.flat_fraction_x * target_half_x_m;
flat_half_y_m = cfg.mraf.flat_fraction_y * target_half_y_m;
free_half_x_m = cfg.mraf.free_fraction_x * target_half_x_m;
free_half_y_m = cfg.mraf.free_fraction_y * target_half_y_m;

flat_core = abs(XF_m) <= flat_half_x_m & abs(YF_m) <= flat_half_y_m;
free_box = abs(XF_m) <= free_half_x_m & abs(YF_m) <= free_half_y_m;
free_edge = free_box & ~flat_core;
outer_suppress = ~free_box;

A_target = nan(size(XF_m));
A_target(flat_core) = 1;
I_target = A_target.^2;

masks = struct();
masks.flat_core = flat_core;
masks.core = flat_core;
masks.signal = flat_core;
masks.finite_target = flat_core;
masks.target = flat_core;
masks.free = free_edge;
masks.free_edge = free_edge;
masks.free_noise = free_edge;
masks.side_lobe_reservoir = free_edge;
masks.outer_suppress = outer_suppress;
masks.noise = free_edge | outer_suppress;
masks.guard = free_box;
masks.edge = free_edge;
masks.core_rect = flat_core;
masks.wgs_update_scale = double(flat_core);
masks.wgs_y_edge_damping_band = false(size(flat_core));
masks.y_free_reservoir = false(size(flat_core));
masks.free_threshold = NaN;
masks.edge_low_threshold_effective = NaN;

target = struct();
target.amplitude = A_target;
target.intensity = I_target;
target.intensity_plot = I_target;
target.intensity_plot(~isfinite(target.intensity_plot)) = 0;
target.I_rd_norm = [];
target.I_target_geometry_norm = double(flat_core);
target.I_target_continuous = target.intensity_plot;
target.I_env = target.intensity_plot;
target.blend_map = double(flat_core);
target.free_mask = free_edge;
target.masks = masks;
target.mode = "flat_core_free_edge";
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
target.free_half_x_m = free_half_x_m;
target.free_half_y_m = free_half_y_m;
target.flat_fraction_x = cfg.mraf.flat_fraction_x;
target.flat_fraction_y = cfg.mraf.flat_fraction_y;
target.free_fraction_x = cfg.mraf.free_fraction_x;
target.free_fraction_y = cfg.mraf.free_fraction_y;
end
