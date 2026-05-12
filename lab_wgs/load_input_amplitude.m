function [amplitude, aperture_mask] = load_input_amplitude(cfg)
% Generate Gaussian input amplitude in the DOE plane.
% References: lab_test_f200mm/src/propagation.py:make_input_gaussian

N = cfg.N;
dx = cfg.dx_doe_m;

x = ((0:N-1) - N/2) * dx;
[X, Y] = meshgrid(x, x);

w = cfg.input_1e2_diameter_m / 2;  % 1/e^2 intensity radius
aperture_radius = cfg.aperture_diameter_m / 2;

r2 = X.^2 + Y.^2;
amplitude = exp(-r2 / w^2);
aperture_mask = r2 <= aperture_radius^2;
amplitude(~aperture_mask) = 0;

% Normalize to unit L2 norm (matches Python normalize_power with target_power=1)
nrm = sqrt(sum(abs(amplitude(:)).^2));
amplitude = amplitude / max(nrm, 1e-20);
