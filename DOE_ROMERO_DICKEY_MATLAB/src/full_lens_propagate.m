function [focal_x_m, focal_y_m, intensity_norm, focal] = full_lens_propagate(U_doe, X_m, Y_m, dx_doe_m, cfg)
% full_lens_propagate 计算 DOE->lens 自由传播、透镜相位、lens->焦平面传播。
%
% 第一步用角谱法传播 cfg.doe_to_lens_m，保持 DOE/lens 面采样不变。
% 第二步乘薄透镜二次相位 exp(-i*k*(x^2+y^2)/(2*f))。
% 第三步用单 FFT Fresnel 传播到 z=f 的焦平面，输出采样为 lambda*f/(N*dx)。

lambda_m = cfg.lambda_m;
f_m = cfg.f_m;
k_rad_m = 2 * pi / lambda_m;

U_lens_in = angular_spectrum_propagate(U_doe, dx_doe_m, lambda_m, cfg.doe_to_lens_m);
if cfg.apply_lens_aperture_in_full_model
    lens_aperture_mask = (X_m.^2 + Y_m.^2) <= cfg.lens_aperture_radius_m^2;
    U_lens_in = U_lens_in .* lens_aperture_mask;
else
    lens_aperture_mask = true(size(U_lens_in));
end
lens_phase = exp(-1i * k_rad_m .* (X_m.^2 + Y_m.^2) ./ (2 * f_m));
U_lens_out = U_lens_in .* lens_phase;

[focal_x_m, focal_y_m, U_focal] = fresnel_propagate_single_fft(U_lens_out, X_m, Y_m, dx_doe_m, lambda_m, f_m);
intensity = abs(U_focal).^2;
max_intensity = max(intensity(:));
if max_intensity > 0
    intensity_norm = intensity ./ max_intensity;
else
    intensity_norm = intensity;
end

focal = struct();
focal.U = U_focal;
focal.intensity_raw = intensity;
focal.intensity_norm = intensity_norm;
focal.focal_dx_m = lambda_m * f_m / (size(U_doe, 1) * dx_doe_m);
focal.normalization = "max=1";
focal.model = "full_doe_to_lens_lens_to_focus";
focal.doe_to_lens_m = cfg.doe_to_lens_m;
focal.apply_lens_aperture = cfg.apply_lens_aperture_in_full_model;
focal.lens_aperture_diameter_m = cfg.lens_aperture_diameter_m;
focal.lens_aperture_grid_fraction = nnz(lens_aperture_mask) / numel(lens_aperture_mask);
end

function [x_out_m, y_out_m, U_out] = fresnel_propagate_single_fft(U_in, X_m, Y_m, dx_m, lambda_m, z_m)
N = size(U_in, 1);
k_rad_m = 2 * pi / lambda_m;
df = 1 / (N * dx_m);
freq = ((-N/2):(N/2 - 1)) * df;
x_out_m = lambda_m * z_m * freq;
y_out_m = x_out_m;

input_chirp = exp(1i * k_rad_m .* (X_m.^2 + Y_m.^2) ./ (2 * z_m));
[X_out_m, Y_out_m] = meshgrid(x_out_m, y_out_m);
output_chirp = exp(1i * k_rad_m .* (X_out_m.^2 + Y_out_m.^2) ./ (2 * z_m));

U_out = exp(1i * k_rad_m * z_m) ./ (1i * lambda_m * z_m) .* output_chirp .* ...
    fftshift(fft2(ifftshift(U_in .* input_chirp))) * dx_m^2;
end
