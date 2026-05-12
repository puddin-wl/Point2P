function [phase_new, weights_new, focal_intensity] = wgs_iteration(phase, input_amp, weights, target, cfg, varargin)
% One iteration of hybrid-field WGS.
%
%   [phase_new, weights_new, I_focal] = wgs_iteration(phase, input_amp, weights, target, cfg)
%   runs pure simulation WGS (FFT both directions).
%
%   [phase_new, weights_new, I_focal] = wgs_iteration(..., I_meas)
%   runs hybrid-field WGS: measured intensity replaces FFT amplitude.
%
%   [phase_new, weights_new, I_focal] = wgs_iteration(..., I_meas, do_update)
%   do_update (logical, default true) controls whether WGS weights are updated
%   this iteration. Set false to skip the weight update (used for update_every).
%
% References: lab_test_f200mm/src/mraf_gs.py WGS path

p = inputParser;
p.addOptional('I_meas', [], @(x) isnumeric(x) && ismatrix(x));
p.addOptional('do_weight_update', true, @islogical);
p.parse(varargin{:});
I_meas = p.Results.I_meas;
do_update = p.Results.do_weight_update;

N = cfg.N;
n_total = N * N;

%% Forward propagation
field_doe = input_amp .* exp(1i * phase);
if isempty(I_meas)
    E_focal = fftshift(fft2(ifftshift(field_doe))) / sqrt(n_total);
    I_focal = abs(E_focal).^2;
else
    I_focal = I_meas;
    % Use FFT phase with measured amplitude
    E_focal_sim = fftshift(fft2(ifftshift(field_doe))) / sqrt(n_total);
    E_focal = sqrt(I_focal) .* exp(1i * angle(E_focal_sim));
end

A_focal = abs(E_focal);

%% Weight update
mask_flat = target.mask_flat;
mask_signal = target.mask_signal;
mask_bg = target.mask_bg;

if do_update
    amp_flat = A_focal(mask_flat);
    if ~isempty(amp_flat) && mean(amp_flat) > 0
        amp_mean = mean(amp_flat);
        ratio = amp_mean ./ max(amp_flat, 1e-12);
        updated = weights(mask_flat) .* ratio.^cfg.wgs_feedback_exponent;
        updated = max(cfg.wgs_weight_min, min(cfg.wgs_weight_max, updated));
        updated = updated / mean(updated);
        weights_new = weights;
        weights_new(mask_flat) = updated;
    else
        weights_new = weights;
    end
else
    weights_new = weights;
end

%% Build weighted target
A_target_eff = target.A_signal;
A_target_eff(mask_flat) = target.A_signal(mask_flat) .* weights_new(mask_flat);

%% Far-field projection (GS-like)
E_proj = E_focal;
% Signal region: replace amplitude with weighted target
E_proj(mask_signal) = A_target_eff(mask_signal) .* exp(1i * angle(E_focal(mask_signal)));
% Background: attenuate
E_proj(mask_bg) = E_proj(mask_bg) * cfg.bg_factor;

%% Inverse propagation
E_doe_new = fftshift(ifft2(ifftshift(E_proj))) * sqrt(n_total);
phase_new = mod(angle(E_doe_new), 2*pi);

%% Re-evaluate with new phase to get output intensity
field_new = input_amp .* exp(1i * phase_new);
E_focal_new = fftshift(fft2(ifftshift(field_new))) / sqrt(n_total);
focal_intensity = abs(E_focal_new).^2;
end
