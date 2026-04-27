function [phase_rad, info] = romero_dickey_phase_1d(x_m, ri_m, Ro_m, lambda_m, f_m)
% romero_dickey_phase_1d 计算 Romero-Dickey Gaussian-to-flattop 1D 解析相位。
%
% xi = x / ri，其中 ri 是论文定义的输入强度下降到 1/e peak 的半径。
% phase_rad = beta * phi_dimless，不要额外再乘 2*pi。

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

