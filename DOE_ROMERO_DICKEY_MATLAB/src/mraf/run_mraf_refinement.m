function result = run_mraf_refinement(input_amplitude, phase0, focal_x_m, focal_y_m, cfg)
% run_mraf_refinement Run MRAF refinement from an RD/Point2P initial phase.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   This function does not generate phase0. It receives phase0 from run_mraf_one.m
%   after build_separable_phase_2d.m has generated the Romero-Dickey initial
%   phase, then uses phase0 as the starting DOE/SLM phase for MRAF iterations.
% Inputs:
%   input_amplitude is the DOE-plane Gaussian amplitude from gaussian_input_field.m.
%   phase0 is the unwrapped RD initial phase in radians, same size as the DOE grid.
%   focal_x_m/focal_y_m define the ideal FFT focal-plane coordinates.
% Outputs:
%   result.rd.phase stores the original phase0; result.mraf.phase stores the
%   refined phase after MRAF. The RD baseline field is input_amplitude.*exp(1i*phase0).
% Physical meaning:
%   MRAF starts from the analytical RD phase and modifies it; phase0 is only an
%   initial condition, not a fixed final solution.
% Used by:
%   run_mraf_one.m and all MRAF probe scripts that call run_mraf_one.m.
% Notes:
%   Do not infer target-generation details from phase0. Current target masks are
%   built later in this function from cfg.mraf.target_mode.

input_power = sum(abs(input_amplitude).^2, 'all');
phase = double(phase0);

% RD baseline propagation: this is the first place phase0 becomes a complex
% DOE field for the ideal FFT model used by MRAF.
U_rd = mraf_forward_fft(input_amplitude .* exp(1i * phase));
I_rd_raw = abs(U_rd).^2;

if strcmpi(string(cfg.mraf.target_mode), "flat_core_free_edge")
    target = make_flat_core_free_edge_target(focal_x_m, focal_y_m, cfg);
elseif strcmpi(string(cfg.mraf.target_mode), "flat_core_guard_lobe")
    target = make_flat_core_guard_lobe_target(focal_x_m, focal_y_m, cfg);
elseif strcmpi(string(cfg.mraf.target_mode), "rect_e2_spec_scaled")
    target = make_rect_e2_spec_scaled_target(focal_x_m, focal_y_m, cfg);
elseif strcmpi(string(cfg.mraf.target_mode), "slmsuite_like_mraf")
    target = make_slmsuite_like_mraf_target(focal_x_m, focal_y_m, cfg);
else
    target = make_rd_derived_target(focal_x_m, focal_y_m, I_rd_raw, cfg);
end
rd_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, I_rd_raw, target.masks, input_power, cfg);

finite_mask = target.masks.finite_target & isfinite(target.amplitude) & target.amplitude > 0;
if isfield(target.masks, 'wgs')
    finite_mask = target.masks.wgs & isfinite(target.amplitude) & target.amplitude > 0;
end
weights = nan(size(target.amplitude));
weights(finite_mask) = 1;
iteration = zeros(cfg.mraf.n_iter, 14);

for iter = 1:cfg.mraf.n_iter
    U0 = input_amplitude .* exp(1i * phase);
    Uf = mraf_forward_fft(U0);
    Af = abs(Uf);
    if isfield(cfg.mraf, 'use_wgs') && cfg.mraf.use_wgs && strcmpi(string(cfg.mraf.method), "WGS-MRAF")
        weights = update_wgs_weights(weights, target.amplitude, Af, finite_mask, cfg);
    end
    Uf_proj = project_farfield_mraf(Uf, target.amplitude, target.masks, weights, cfg);
    U_back = mraf_backward_fft(Uf_proj);
    phase_candidate = angle(U_back);
    phase = angle((1 - cfg.mraf.phase_blend) .* exp(1i * phase) + cfg.mraf.phase_blend .* exp(1i * phase_candidate));

    Uf_iter = mraf_forward_fft(input_amplitude .* exp(1i * phase));
    iter_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, abs(Uf_iter).^2, target.masks, input_power, cfg);
    iteration(iter, :) = [iter, iter_metrics.flat_core_rms, iter_metrics.flat_core_pv, iter_metrics.flat_core_uniformity, ...
        iter_metrics.side_lobe_energy_ratio, iter_metrics.outer_energy_ratio, iter_metrics.size50_x_um, iter_metrics.size50_y_um, ...
        iter_metrics.transition_width_x_um, iter_metrics.transition_width_y_um, iter_metrics.core_rms, iter_metrics.shoulder_peak_x, ...
        iter_metrics.shoulder_peak_y, iter_metrics.outer_tail_peak_y_rel_to_core];
end

U_mraf = mraf_forward_fft(input_amplitude .* exp(1i * phase));
I_mraf_raw = abs(U_mraf).^2;
mraf_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, I_mraf_raw, target.masks, input_power, cfg);

result = struct();
result.target = target;
result.rd.U = U_rd;
result.rd.intensity_raw = I_rd_raw;
result.rd.intensity_norm = I_rd_raw ./ max(mean(I_rd_raw(target.masks.flat_core), 'omitnan'), eps);
result.rd.phase = phase0;
result.rd.phase_wrapped = wrap_mraf_phase(phase0);
result.rd.metrics = rd_metrics;
result.mraf.U = U_mraf;
result.mraf.intensity_raw = I_mraf_raw;
result.mraf.intensity_norm = I_mraf_raw ./ max(mean(I_mraf_raw(target.masks.flat_core), 'omitnan'), eps);
result.mraf.phase = phase;
result.mraf.phase_wrapped = wrap_mraf_phase(phase);
result.mraf.metrics = mraf_metrics;
result.iteration_metrics = array2table(iteration, 'VariableNames', {'iter','flat_core_rms','flat_core_pv','flat_core_uniformity','side_lobe_energy_ratio','outer_energy_ratio','size50_x_um','size50_y_um','transition_width_x_um','transition_width_y_um','core_rms','shoulder_peak_x','shoulder_peak_y','outer_tail_y'});
end
