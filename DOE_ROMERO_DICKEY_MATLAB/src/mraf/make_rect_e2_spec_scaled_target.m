function target = make_rect_e2_spec_scaled_target(focal_x_m, focal_y_m, cfg)
% make_rect_e2_spec_scaled_target Build scaled industrial e^-2 rectangular target.
%
% The target is defined directly in the focal plane by three rectangular
% iso-intensity contours scaled from the industrial template:
%   90%   : 285.04 x  97.94 um
%   50%   : 330.00 x 120.00 um
%   13.5% : 370.21 x 147.12 um
% Outside the 13.5% rectangle is don't-care (NaN target amplitude).

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
abs_x = abs(XF_m);
abs_y = abs(YF_m);

[size90_x_um, size90_y_um, size50_x_um, size50_y_um, size135_x_um, size135_y_um] = rect_e2_scaled_sizes_um(cfg);

half90_x_m = 0.5 * size90_x_um * 1e-6;
half90_y_m = 0.5 * size90_y_um * 1e-6;
half50_x_m = 0.5 * size50_x_um * 1e-6;
half50_y_m = 0.5 * size50_y_um * 1e-6;
half135_x_m = 0.5 * size135_x_um * 1e-6;
half135_y_m = 0.5 * size135_y_um * 1e-6;

inside90 = abs_x <= half90_x_m & abs_y <= half90_y_m;
inside50 = abs_x <= half50_x_m & abs_y <= half50_y_m;
inside135 = abs_x <= half135_x_m & abs_y <= half135_y_m;
transition90_50 = inside50 & ~inside90;
transition50_135 = inside135 & ~inside50;
transition = transition90_50 | transition50_135;
outside_e2 = ~inside135;

I_target = nan(size(XF_m));
I_target(inside90) = 1.0;

s90_50 = rectangular_band_coordinate(abs_x, abs_y, half90_x_m, half90_y_m, half50_x_m, half50_y_m);
u90_50 = smoothstep01(s90_50);
I_target(transition90_50) = 0.9 + (0.5 - 0.9) .* u90_50(transition90_50);

s50_135 = rectangular_band_coordinate(abs_x, abs_y, half50_x_m, half50_y_m, half135_x_m, half135_y_m);
u50_135 = smoothstep01(s50_135);
I_target(transition50_135) = 0.5 + (0.135 - 0.5) .* u50_135(transition50_135);

A_target = sqrt(max(I_target, 0));
A_target(outside_e2) = NaN;

masks = struct();
masks.plateau = inside90;
masks.flat_core = inside90;
masks.core = inside90;
masks.signal = inside90;
masks.finite_target = inside90;
masks.target = inside90;
masks.transition = transition;
masks.transition_90_50 = transition90_50;
masks.transition_50_135 = transition50_135;
masks.edge = transition;
masks.free = false(size(XF_m));
masks.free_edge = false(size(XF_m));
masks.free_noise = false(size(XF_m));
masks.side_lobe_reservoir = false(size(XF_m));
masks.outside_e2 = outside_e2;
masks.outer_suppress = outside_e2;
masks.noise = outside_e2;
masks.guard = inside135;
masks.target_13p5 = inside135;

target = struct();
target.amplitude = A_target;
target.intensity = I_target;
target.intensity_plot = I_target;
target.blend_map = double(inside90);
target.masks = masks;
target.mode = "rect_e2_spec_scaled";
target.remap_anchor_50 = NaN;
target.spec = struct();
target.spec.size90_x_um = size90_x_um;
target.spec.size90_y_um = size90_y_um;
target.spec.size50_x_um = size50_x_um;
target.spec.size50_y_um = size50_y_um;
target.spec.size135_x_um = size135_x_um;
target.spec.size135_y_um = size135_y_um;
target.spec.transition_width_135_90_x_um = 0.5 * (size135_x_um - size90_x_um);
target.spec.transition_width_135_90_y_um = 0.5 * (size135_y_um - size90_y_um);
target.spec.transition_cap_blend = get_mraf_field(cfg, 'transition_cap_blend', 0.35);
target.spec.outer_factor = get_mraf_field(cfg, 'outer_factor', 0.95);
end

function [size90_x_um, size90_y_um, size50_x_um, size50_y_um, size135_x_um, size135_y_um] = rect_e2_scaled_sizes_um(cfg)
sample50_x_um = 36.11;
sample50_y_um = 22.74;
sample135_x_um = 40.51;
sample135_y_um = 27.88;
sample_transition_um = 4.66;

size50_x_um = cfg.target_size_x_m * 1e6;
size50_y_um = cfg.target_size_y_m * 1e6;
scale_x = size50_x_um / sample50_x_um;
scale_y = size50_y_um / sample50_y_um;
size135_x_um = sample135_x_um * scale_x;
size135_y_um = sample135_y_um * scale_y;
size90_x_um = (sample135_x_um - 2 * sample_transition_um) * scale_x;
size90_y_um = (sample135_y_um - 2 * sample_transition_um) * scale_y;
end

function s = rectangular_band_coordinate(abs_x, abs_y, inner_half_x, inner_half_y, outer_half_x, outer_half_y)
sx = (abs_x - inner_half_x) ./ max(outer_half_x - inner_half_x, eps);
sy = (abs_y - inner_half_y) ./ max(outer_half_y - inner_half_y, eps);
s = max(sx, sy);
s = min(max(s, 0), 1);
end

function y = smoothstep01(x)
x = min(max(x, 0), 1);
y = x.^2 .* (3 - 2 .* x);
end

function value = get_mraf_field(cfg, name, default_value)
if isfield(cfg, 'mraf') && isfield(cfg.mraf, name)
    value = cfg.mraf.(name);
else
    value = default_value;
end
end
