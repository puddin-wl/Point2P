function target = make_slmsuite_like_mraf_target(focal_x_m, focal_y_m, cfg)
% make_slmsuite_like_mraf_target Build a slmsuite-style MRAF amplitude target.
%
% Semantics follow slmsuite/holodyne-style MRAF:
%   finite nonzero target amplitude: optimized signal region
%   NaN target amplitude: MRAF noise region, retains mraf_factor * current amplitude
%   zero target amplitude: null/background, projected to zero

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);

signal_size_x_um = cfg.target_size_x_m * 1e6;
signal_size_y_um = cfg.target_size_y_m * 1e6;
sample50_x_um = 36.11;
sample50_y_um = 22.74;
sample135_x_um = 40.51;
sample135_y_um = 27.88;
scale_x = signal_size_x_um / sample50_x_um;
scale_y = signal_size_y_um / sample50_y_um;
noise_size_x_um = sample135_x_um * scale_x;
noise_size_y_um = sample135_y_um * scale_y;

signal_half_x_m = 0.5 * signal_size_x_um * 1e-6;
signal_half_y_m = 0.5 * signal_size_y_um * 1e-6;
noise_half_x_m = 0.5 * noise_size_x_um * 1e-6;
noise_half_y_m = 0.5 * noise_size_y_um * 1e-6;

signal = abs(XF_m) <= signal_half_x_m & abs(YF_m) <= signal_half_y_m;
noise_box = abs(XF_m) <= noise_half_x_m & abs(YF_m) <= noise_half_y_m;
noise = noise_box & ~signal;
null = ~noise_box;

A_target = zeros(size(XF_m));
A_target(signal) = 1;
A_target(noise) = NaN;
A_target(null) = 0;

masks = struct();
masks.signal = signal;
masks.core = signal;
masks.flat_core = signal;
masks.noise = noise;
masks.free = noise;
masks.free_edge = noise;
masks.null = null;
masks.outer_suppress = null;
masks.finite_target = isfinite(A_target);
masks.wgs = signal;
masks.guard = noise_box;
masks.noise_box = noise_box;
masks.edge = noise;

target = struct();
target.amplitude = A_target;
target.intensity = A_target.^2;
target.intensity_plot = target.intensity;
target.blend_map = double(signal);
target.masks = masks;
target.mode = "slmsuite_like_mraf";
target.remap_anchor_50 = NaN;
target.spec = struct();
target.spec.signal_size_x_um = signal_size_x_um;
target.spec.signal_size_y_um = signal_size_y_um;
target.spec.noise_size_x_um = noise_size_x_um;
target.spec.noise_size_y_um = noise_size_y_um;
target.spec.signal_half_x_um = signal_half_x_m * 1e6;
target.spec.signal_half_y_um = signal_half_y_m * 1e6;
target.spec.noise_half_x_um = noise_half_x_m * 1e6;
target.spec.noise_half_y_um = noise_half_y_m * 1e6;
end
