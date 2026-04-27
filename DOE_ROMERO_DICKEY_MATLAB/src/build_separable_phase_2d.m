function [phase_unwrapped_rad, phase_wrapped_rad, info] = build_separable_phase_2d(X_m, Y_m, aperture_mask, cfg)
% build_separable_phase_2d 由 x/y 两个 Romero-Dickey 1D 相位构造 2D separable 相位。

[phase_x_rad, info_x] = romero_dickey_phase_1d(X_m, cfg.input_1e_radius_m, cfg.Ro_x_m, cfg.lambda_m, cfg.f_m);
[phase_y_rad, info_y] = romero_dickey_phase_1d(Y_m, cfg.input_1e_radius_m, cfg.Ro_y_m, cfg.lambda_m, cfg.f_m);

phase_unwrapped_rad = cfg.phase_sign .* (cfg.phase_scale_x .* phase_x_rad + cfg.phase_scale_y .* phase_y_rad);
phase_unwrapped_rad(~aperture_mask) = 0;
phase_wrapped_rad = wrap_phase_2pi(phase_unwrapped_rad);
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

