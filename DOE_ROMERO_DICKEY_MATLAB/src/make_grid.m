function [X_m, Y_m, x_m, y_m, grid] = make_grid(cfg)
% make_grid Generate the centered DOE/SLM-plane physical coordinate grid.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Provide the physical x/y coordinate arrays used by the initial RD phase
%   formula and by the input Gaussian/aperture mask.
% Inputs:
%   cfg.N:
%     Computational grid size. Default is 2048 x 2048 samples.
%   cfg.dx_doe_m:
%     DOE-plane sample pitch in meters. It is derived in default_config.m from
%     the requested focal-plane sampling via dx = lambda*f/(N*focal_dx).
% Outputs:
%   x_m, y_m:
%     Centered 1D DOE-plane coordinates in meters.
%   X_m, Y_m:
%     2D meshgrid coordinates used by gaussian_input_field.m and
%     build_separable_phase_2d.m.
% Physical meaning:
%   This is the computational DOE plane. Its extent is not the same concept as
%   the 15 mm clear aperture or the 5 mm illuminated beam diameter; those are
%   applied later as masks/amplitudes on this grid.
% Used by:
%   run_one.m, run_mraf_one.m, export_phase_for_python.m, and
%   run_initial_phase_only.m.
% Notes:
%   Changing N or dx_doe_m changes the FFT/focal-plane coordinate relation and
%   therefore changes the sampled phase and propagation result.

N = cfg.N;
dx_m = cfg.dx_doe_m;
index = (-N/2):(N/2 - 1);
x_m = index * dx_m;
y_m = index * dx_m;
[X_m, Y_m] = meshgrid(x_m, y_m);

grid = struct();
grid.N = N;
grid.dx_doe_m = dx_m;
grid.extent_m = N * dx_m;
grid.x_min_m = min(x_m);
grid.x_max_m = max(x_m);
grid.y_min_m = min(y_m);
grid.y_max_m = max(y_m);
end
