function [phase_rad, info] = romero_dickey_phase_1d(x_m, ri_m, Ro_m, lambda_m, f_m)
% romero_dickey_phase_1d Compute the 1D Romero-Dickey Gaussian-to-flat-top phase.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Evaluate the analytical 1D RD stationary-phase formula used to construct
%   the separable 2D initial phase0.
% Inputs:
%   x_m:
%     DOE-plane coordinate in meters. This can be a vector or a 2D coordinate
%     array from meshgrid.
%   ri_m:
%     Gaussian 1/e amplitude radius used by the RD paper convention. In the
%     default config it is derived from the 5 mm 1/e^2 intensity diameter as
%     cfg.input_1e2_radius_m / sqrt(2).
%   Ro_m:
%     RD output scale for one target axis. In the default config Ro_x/Ro_y are
%     target_size_x/y / sqrt(pi), so the 330 x 120 um target enters here.
%   lambda_m, f_m:
%     Wavelength and Fourier-lens focal length.
% Outputs:
%   phase_rad:
%     Unwrapped phase contribution in radians for one axis.
% Physical meaning:
%   The dimensionless function phi_dimless implements the RD stationary-phase
%   Gaussian-to-flat-top mapping. beta = 2*pi*ri*Ro/(lambda*f) scales that
%   mapping by optical wavelength, focal length, input beam size, and target
%   output size.
% Used by:
%   build_separable_phase_2d.m.
% Notes:
%   This formula is analytical; it is not an iterative point-by-point mapping
%   in this MATLAB implementation.

% Normalized input coordinate. xi expresses DOE coordinate in units of the
% Gaussian 1/e amplitude radius, not in units of the clear aperture.
xi = x_m ./ ri_m;

% Romero-Dickey dimensionless phase primitive for Gaussian-to-flat-top beam
% shaping. The erf term describes cumulative Gaussian energy remapping; the
% exp(-xi^2) term completes the stationary-phase integral.
phi_dimless = xi .* sqrt(pi) ./ 2 .* erf(xi) + 0.5 .* exp(-xi.^2) - 0.5;

% beta contains all physical scaling: input beam radius, requested output
% scale, wavelength, and focal length. Do not multiply by 2*pi again later.
beta = 2 * pi * ri_m * Ro_m / (lambda_m * f_m);
phase_rad = beta .* phi_dimless;

info = struct();
info.ri_m = ri_m;
info.Ro_m = Ro_m;
info.beta = beta;
info.phi_dimless_min = min(phi_dimless(:));
info.phi_dimless_max = max(phi_dimless(:));
end
