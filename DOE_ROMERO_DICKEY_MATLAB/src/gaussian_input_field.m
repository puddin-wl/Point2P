function [field, amplitude, intensity, aperture_mask] = gaussian_input_field(X_m, Y_m, cfg)
% gaussian_input_field 生成 DOE 面 5 mm 1/e^2 强度直径高斯入射光场。
%
% 注意：15 mm 是 DOE clear aperture，不是入射光斑直径。
% 5 mm 是高斯光斑 1/e^2 强度直径；若 w 为 1/e^2 强度半径，
% 则 I = exp(-2*r^2/w^2)，amplitude = exp(-r^2/w^2)。

r2_m2 = X_m.^2 + Y_m.^2;
w_m = cfg.input_1e2_radius_m;
intensity = exp(-2 * r2_m2 / w_m^2);
amplitude = exp(-r2_m2 / w_m^2);
aperture_mask = r2_m2 <= cfg.aperture_radius_m^2;
field = amplitude .* aperture_mask;
intensity = intensity .* aperture_mask;
amplitude = amplitude .* aperture_mask;
end

