function phase_data = generate_initial_phase(cfg, varargin)
% generate_initial_phase Generate the standalone Romero-Dickey initial phase.
%
% This module is intentionally separated from DOE_ROMERO_DICKEY_MATLAB/src so
% downstream code can request phase0 from one clear location instead of
% re-implementing the phase-generation chain inside MRAF scripts.
%
% Usage:
%   phase_data = generate_initial_phase(cfg)
%   phase_data = generate_initial_phase(cfg, 'do_forward', true)
%   phase_data = generate_initial_phase(cfg, 'output_dir', out_dir)
%
% Inputs:
%   cfg must contain the same phase/grid fields as default_config.m:
%   N, dx_doe_m, lambda_m, f_m, aperture_radius_m, input_1e2_radius_m,
%   input_1e_radius_m, Ro_x_m, Ro_y_m, phase_sign, phase_scale_x,
%   phase_scale_y, phase_method, focal_dx_m.
%
% Outputs:
%   phase_data.phase0_unwrapped_rad: unwrapped RD phase0 in radians.
%   phase_data.phase0_wrapped_rad: phase0 wrapped to [0, 2*pi), NaN outside aperture.
%   phase_data.input_amplitude: DOE-plane Gaussian amplitude.
%   phase_data.aperture_mask: clear aperture mask.
%   phase_data.focal_*: optional ideal FFT baseline if do_forward is true.

if nargin < 1 || isempty(cfg)
    cfg = default_initial_phase_config(fileparts(mfilename('fullpath')));
end
opts = parse_options(varargin{:});

[X_m, Y_m, x_m, y_m, grid] = ipg_make_grid(cfg);
[input_field, input_amplitude, input_intensity, aperture_mask] = ipg_gaussian_input_field(X_m, Y_m, cfg);
[phase0_unwrapped_rad, phase0_wrapped_rad, phase_info] = ipg_build_separable_phase_2d(X_m, Y_m, aperture_mask, cfg);

phase_data = struct();
phase_data.cfg = cfg;
phase_data.grid = grid;
phase_data.X_m = X_m;
phase_data.Y_m = Y_m;
phase_data.x_m = x_m;
phase_data.y_m = y_m;
phase_data.input_field = input_field;
phase_data.input_amplitude = input_amplitude;
phase_data.input_intensity = input_intensity;
phase_data.aperture_mask = aperture_mask;
phase_data.phase0_unwrapped_rad = phase0_unwrapped_rad;
phase_data.phase0_wrapped_rad = phase0_wrapped_rad;
phase_data.phase_info = phase_info;

if opts.do_forward
    doe_field = input_field .* exp(1i * phase0_unwrapped_rad);
    [focal_x_m, focal_y_m, intensity_norm, focal] = ipg_fresnel_to_focal_fft(doe_field, cfg.dx_doe_m, cfg.lambda_m, cfg.f_m);
    phase_data.focal_x_m = focal_x_m;
    phase_data.focal_y_m = focal_y_m;
    phase_data.initial_intensity_norm = intensity_norm;
    phase_data.focal = focal;
else
    focal_index = (-cfg.N/2):(cfg.N/2 - 1);
    phase_data.focal_x_m = focal_index * cfg.focal_dx_m;
    phase_data.focal_y_m = phase_data.focal_x_m;
end

if strlength(opts.output_dir) > 0
    ipg_save_outputs(char(opts.output_dir), phase_data, opts);
end
end

function opts = parse_options(varargin)
opts = struct();
opts.do_forward = false;
opts.output_dir = "";
opts.save_png = true;
opts.figure_dpi = 150;
if mod(numel(varargin), 2) ~= 0
    error('generate_initial_phase options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k + 1};
    switch name
        case "do_forward"
            opts.do_forward = logical(value);
        case "output_dir"
            opts.output_dir = string(value);
        case "save_png"
            opts.save_png = logical(value);
        case "figure_dpi"
            opts.figure_dpi = value;
        otherwise
            error('Unknown generate_initial_phase option: %s', name);
    end
end
end

function [X_m, Y_m, x_m, y_m, grid] = ipg_make_grid(cfg)
N = cfg.N;
dx_m = cfg.dx_doe_m;
index = (-N/2):(N/2 - 1);
x_m = index * dx_m;
y_m = index * dx_m;
[X_m, Y_m] = meshgrid(x_m, y_m);
grid = struct('N', N, 'dx_doe_m', dx_m, 'extent_m', N * dx_m, ...
    'x_min_m', min(x_m), 'x_max_m', max(x_m), 'y_min_m', min(y_m), 'y_max_m', max(y_m));
end

function [field, amplitude, intensity, aperture_mask] = ipg_gaussian_input_field(X_m, Y_m, cfg)
r2_m2 = X_m.^2 + Y_m.^2;
w_m = cfg.input_1e2_radius_m;
intensity = exp(-2 * r2_m2 / w_m^2);
amplitude = exp(-r2_m2 / w_m^2);
aperture_mask = r2_m2 <= cfg.aperture_radius_m^2;
field = amplitude .* aperture_mask;
intensity = intensity .* aperture_mask;
amplitude = amplitude .* aperture_mask;
end

function [phase_unwrapped_rad, phase_wrapped_rad, info] = ipg_build_separable_phase_2d(X_m, Y_m, aperture_mask, cfg)
[phase_x_rad, info_x] = ipg_romero_dickey_phase_1d(X_m, cfg.input_1e_radius_m, cfg.Ro_x_m, cfg.lambda_m, cfg.f_m);
[phase_y_rad, info_y] = ipg_romero_dickey_phase_1d(Y_m, cfg.input_1e_radius_m, cfg.Ro_y_m, cfg.lambda_m, cfg.f_m);
phase_unwrapped_rad = cfg.phase_sign .* (cfg.phase_scale_x .* phase_x_rad + cfg.phase_scale_y .* phase_y_rad);
phase_unwrapped_rad(~aperture_mask) = 0;
phase_wrapped_rad = mod(phase_unwrapped_rad, 2*pi);
phase_wrapped_rad(~aperture_mask) = NaN;
info = struct();
info.method = cfg.phase_method;
info.phase_sign = cfg.phase_sign;
info.phase_scale_x = cfg.phase_scale_x;
info.phase_scale_y = cfg.phase_scale_y;
info.x = info_x;
info.y = info_y;
info.unwrapped_min_rad = min(phase_unwrapped_rad(aperture_mask), [], 'all');
info.unwrapped_max_rad = max(phase_unwrapped_rad(aperture_mask), [], 'all');
end

function [phase_rad, info] = ipg_romero_dickey_phase_1d(x_m, ri_m, Ro_m, lambda_m, f_m)
xi = x_m ./ ri_m;
phi_dimless = xi .* sqrt(pi) ./ 2 .* erf(xi) + 0.5 .* exp(-xi.^2) - 0.5;
beta = 2 * pi * ri_m * Ro_m / (lambda_m * f_m);
phase_rad = beta .* phi_dimless;
info = struct();
info.ri_m = ri_m;
info.Ro_m = Ro_m;
info.beta = beta;
info.phi_dimless_min = min(phi_dimless(:));
info.phi_dimless_max = max(phi_dimless(:));
end

function [focal_x_m, focal_y_m, intensity_norm, focal] = ipg_fresnel_to_focal_fft(U_doe, dx_doe_m, lambda_m, f_m)
N = size(U_doe, 1);
df = 1 / (N * dx_doe_m);
freq = ((-N/2):(N/2 - 1)) * df;
focal_x_m = lambda_m * f_m * freq;
focal_y_m = focal_x_m;
U_focal = fftshift(fft2(ifftshift(U_doe))) * dx_doe_m^2 / (1i * lambda_m * f_m);
intensity = abs(U_focal).^2;
max_intensity = max(intensity(:));
if max_intensity > 0
    intensity_norm = intensity ./ max_intensity;
else
    intensity_norm = intensity;
end
focal = struct('U', U_focal, 'intensity_raw', intensity, 'intensity_norm', intensity_norm, ...
    'focal_dx_m', lambda_m * f_m / (N * dx_doe_m), 'normalization', "max=1");
end

function ipg_save_outputs(output_dir, phase_data, opts)
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
phase0_unwrapped_rad = phase_data.phase0_unwrapped_rad; %#ok<NASGU>
phase0_wrapped_rad = phase_data.phase0_wrapped_rad; %#ok<NASGU>
phase_info = phase_data.phase_info; %#ok<NASGU>
x_m = phase_data.x_m; y_m = phase_data.y_m; %#ok<NASGU>
focal_x_m = phase_data.focal_x_m; focal_y_m = phase_data.focal_y_m; %#ok<NASGU>
save(fullfile(output_dir, 'phase0.mat'), 'phase0_unwrapped_rad', 'phase0_wrapped_rad', 'phase_info', 'x_m', 'y_m', 'focal_x_m', 'focal_y_m', '-v7.3');
cfg = phase_data.cfg; grid = phase_data.grid; %#ok<NASGU>
save(fullfile(output_dir, 'config_snapshot.mat'), 'cfg', 'grid', '-v7.3');
if opts.save_png
    ipg_save_phase_png(fullfile(output_dir, 'phase0.png'), phase_data.phase0_wrapped_rad, opts.figure_dpi);
    if isfield(phase_data, 'initial_intensity_norm')
        ipg_save_intensity_png(fullfile(output_dir, 'initial_intensity.png'), phase_data.focal_x_m, phase_data.focal_y_m, phase_data.initial_intensity_norm, opts.figure_dpi);
        ipg_save_profile_png(fullfile(output_dir, 'initial_x_profile.png'), phase_data.focal_x_m, phase_data.initial_intensity_norm(round(numel(phase_data.focal_y_m)/2)+1, :), 'x', opts.figure_dpi);
        ipg_save_profile_png(fullfile(output_dir, 'initial_y_profile.png'), phase_data.focal_y_m, phase_data.initial_intensity_norm(:, round(numel(phase_data.focal_x_m)/2)+1).', 'y', opts.figure_dpi);
    end
end
end

function ipg_save_phase_png(path, phase_wrapped, dpi)
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(phase_wrapped); axis image off; colormap(gca, 'hsv'); colorbar;
title('Initial RD phase0 wrapped [0, 2\pi)', 'Interpreter', 'tex');
exportgraphics(fig, path, 'Resolution', dpi); close(fig);
end

function ipg_save_intensity_png(path, focal_x_m, focal_y_m, intensity_norm, dpi)
roi_x = abs(focal_x_m) <= 1.8 * 165e-6;
roi_y = abs(focal_y_m) <= 1.8 * 60e-6;
fig = figure('Visible', 'off', 'Color', 'w');
imagesc(focal_x_m(roi_x)*1e6, focal_y_m(roi_y)*1e6, intensity_norm(roi_y, roi_x));
axis image; colormap(gca, 'turbo'); colorbar; xlabel('x (um)'); ylabel('y (um)');
title('Initial RD focal intensity, ideal FFT');
exportgraphics(fig, path, 'Resolution', dpi); close(fig);
end

function ipg_save_profile_png(path, axis_m, profile, direction, dpi)
fig = figure('Visible', 'off', 'Color', 'w');
plot(axis_m*1e6, profile, 'b-', 'LineWidth', 1.2); grid on;
xlabel(sprintf('%s (um)', direction)); ylabel('normalized intensity');
xlim([-350 350]); ylim([0 1.2]); yline(0.5, 'k--'); yline(0.135, 'm--'); yline(0.9, 'g--');
title(sprintf('Initial RD %s center profile', upper(direction)));
exportgraphics(fig, path, 'Resolution', dpi); close(fig);
end
