function masks = make_mraf_masks(focal_x_m, focal_y_m, I_norm, cfg)
% make_mraf_masks Build RD-derived masks with an explicit free/noise threshold.
[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
core_rect = abs(XF_m) <= cfg.core_fraction * cfg.target_size_x_m & ...
    abs(YF_m) <= cfg.core_fraction * cfg.target_size_y_m;
guard_rect = abs(XF_m) <= cfg.guard_fraction * cfg.target_half_x_m & ...
    abs(YF_m) <= cfg.guard_fraction * cfg.target_half_y_m;

free_threshold = cfg.mraf.free_threshold;
edge_low_threshold = free_threshold;
core_mask = I_norm >= cfg.mraf.core_threshold & guard_rect;
signal_mask = I_norm >= cfg.mraf.signal_threshold & guard_rect;
target_mask = I_norm >= free_threshold & guard_rect;
free_noise_mask = I_norm < free_threshold & guard_rect;
outer_suppress_mask = ~guard_rect;
if isfield(cfg.mraf, 'y_free_reservoir_enable') && cfg.mraf.y_free_reservoir_enable
    reservoir_half_y_m = cfg.mraf.y_free_reservoir_half_y_um * 1e-6;
    reservoir_max_intensity = cfg.mraf.y_free_reservoir_max_intensity;
    y_reservoir = abs(YF_m) >= reservoir_half_y_m & I_norm < reservoir_max_intensity & guard_rect;
    free_noise_mask = free_noise_mask | y_reservoir;
    target_mask = target_mask & ~y_reservoir;
else
    y_reservoir = false(size(free_noise_mask));
end
edge_mask = I_norm >= edge_low_threshold & I_norm < cfg.mraf.core_threshold & guard_rect;
edge_mask = edge_mask & target_mask;
noise_mask = free_noise_mask;
wgs_update_scale = ones(size(I_norm));
if isfield(cfg.mraf, 'wgs_y_edge_damping_enable') && cfg.mraf.wgs_y_edge_damping_enable
    damping_band = I_norm >= cfg.mraf.wgs_y_edge_damping_low & ...
        I_norm <= cfg.mraf.wgs_y_edge_damping_high & guard_rect;
    damping_band = damping_band & abs(YF_m) >= cfg.target_half_y_m * 0.70;
    wgs_update_scale(damping_band) = cfg.mraf.wgs_y_edge_damping_factor;
else
    damping_band = false(size(I_norm));
end

if nnz(core_mask) < 16
    core_mask = core_rect;
end
if nnz(signal_mask) < 16
    signal_mask = I_norm >= cfg.mraf.signal_threshold & guard_rect;
end

masks = struct();
masks.core = core_mask;
masks.signal = signal_mask;
masks.edge = edge_mask;
masks.target = target_mask;
masks.free = free_noise_mask;
masks.free_noise = free_noise_mask;
masks.outer_suppress = outer_suppress_mask;
masks.noise = noise_mask;
masks.guard = guard_rect;
masks.finite_target = target_mask;
masks.core_rect = core_rect;
masks.free_threshold = free_threshold;
masks.edge_low_threshold_effective = edge_low_threshold;
masks.y_free_reservoir = y_reservoir;
masks.wgs_update_scale = wgs_update_scale;
masks.wgs_y_edge_damping_band = damping_band;
end
