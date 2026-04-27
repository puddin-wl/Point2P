function [focal_x_m, focal_y_m, intensity_norm, focal] = fresnel_to_focal_fft(U_doe, dx_doe_m, lambda_m, f_m)
% fresnel_to_focal_fft 使用 Fourier lens 焦平面关系传播 DOE 后场。
%
% U_focal(xf,yf) proportional to FFT{U_doe(x,y)}。
% 坐标采样为 focal_dx = lambda*f/(N*dx_doe)。这里默认 DOE 已位于
% 有效 pupil/Fourier transform plane，不把 DOE 到 lens 的距离当焦距。

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

focal = struct();
focal.U = U_focal;
focal.intensity_raw = intensity;
focal.intensity_norm = intensity_norm;
focal.focal_dx_m = lambda_m * f_m / (N * dx_doe_m);
focal.normalization = "max=1";
end

