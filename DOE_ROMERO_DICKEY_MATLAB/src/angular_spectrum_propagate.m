function U_out = angular_spectrum_propagate(U_in, dx_m, lambda_m, z_m)
% angular_spectrum_propagate 使用角谱法传播标量复光场。
%
% 这是备用函数，可用于后续 DOE 到 lens 面自由传播、乘透镜二次相位、
% 再传播到焦平面的 full model。run_one 默认不调用该模型。

N = size(U_in, 1);
fx = ((-N/2):(N/2 - 1)) / (N * dx_m);
[FX, FY] = meshgrid(fx, fx);
k = 2 * pi / lambda_m;
argument = 1 - (lambda_m * FX).^2 - (lambda_m * FY).^2;
argument(argument < 0) = 0;
H = exp(1i * k * z_m .* sqrt(argument));
U_out = fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(U_in))) .* H)));
end

