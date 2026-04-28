function [phase_unwrapped_rad, phase_wrapped_rad, info] = build_separable_phase_2d(X_m, Y_m, aperture_mask, cfg)
% build_separable_phase_2d Build the separable Romero-Dickey 2D phase.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Build the phase0 used by the baseline RD propagation and by MRAF as its
%   initial DOE/SLM phase.
% Inputs:
%   X_m, Y_m:
%     DOE/SLM-plane physical coordinate arrays in meters from make_grid.m.
%   aperture_mask:
%     Clear-aperture mask. In the default config this is 15 mm diameter. This
%     is the clear/computational aperture, not the illuminated Gaussian beam
%     diameter. The illuminated input beam diameter is cfg.input_1e2_diameter_m.
%   cfg:
%     Supplies lambda_m, f_m, input_1e_radius_m, Ro_x_m/Ro_y_m from the
%     requested 330 x 120 um target, plus phase sign/scale convention fields.
% Outputs:
%   phase_unwrapped_rad:
%     Unwrapped phase0 in radians. This is the phase used in exp(1i*phase0).
%   phase_wrapped_rad:
%     phase0 wrapped to [0, 2*pi) for DOE/SLM export and visualization; values
%     outside the aperture are NaN for plotting.
% Physical meaning:
%   This is an analytical separable Romero-Dickey stationary-phase initial
%   phase. It maps the input Gaussian beam toward a rectangular flat-top by
%   adding an x-only 1D RD mapping and a y-only 1D RD mapping. It is not an
%   iterative point-to-point solver in the current MATLAB code path.
% Used by:
%   run_one.m, run_mraf_one.m, export_phase_for_python.m, and
%   run_initial_phase_only.m.
% Notes:
%   The RD phase encodes target size through cfg.Ro_x_m/cfg.Ro_y_m. It does not
%   encode later MRAF target details such as noise/null regions or target masks.

% Evaluate the same 1D Romero-Dickey formula on the 2D coordinate arrays.
% phase_x_rad varies only with X_m and uses the x target scale Ro_x_m; phase_y
% varies only with Y_m and uses Ro_y_m. Both share the same Gaussian ri.
[phase_x_rad, info_x] = romero_dickey_phase_1d(X_m, cfg.input_1e_radius_m, cfg.Ro_x_m, cfg.lambda_m, cfg.f_m);
[phase_y_rad, info_y] = romero_dickey_phase_1d(Y_m, cfg.input_1e_radius_m, cfg.Ro_y_m, cfg.lambda_m, cfg.f_m);

% Separable 2D phase. phase_sign handles the Fourier sign convention; the two
% phase_scale fields are explicit tuning multipliers but default to 1.
phase_unwrapped_rad = cfg.phase_sign .* (cfg.phase_scale_x .* phase_x_rad + cfg.phase_scale_y .* phase_y_rad);

% Outside the clear aperture the input amplitude is zero. Keep the unwrapped
% phase finite there so U=input_amplitude.*exp(1i*phase) remains well defined.
phase_unwrapped_rad(~aperture_mask) = 0;

% Wrapped phase is for display/export. It is intentionally not used to remove
% the unwrapped phase from the numerical propagation path.
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
