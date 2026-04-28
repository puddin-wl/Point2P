function result = run_mraf_refinement(input_amplitude, phase0, focal_x_m, focal_y_m, cfg)
% run_mraf_refinement Lightweight RD-derived GS/WGS-MRAF refinement.
input_power = sum(abs(input_amplitude).^2, 'all');
phase = double(phase0);
U_rd = mraf_forward_fft(input_amplitude .* exp(1i * phase));
I_rd_raw = abs(U_rd).^2;
target = make_rd_derived_target(focal_x_m, focal_y_m, I_rd_raw, cfg);
rd_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, I_rd_raw, target.masks, input_power, cfg);

weights = ones(size(input_amplitude));
finite_mask = target.masks.finite_target & isfinite(target.amplitude);
iteration = zeros(cfg.mraf.n_iter, 12);

for iter = 1:cfg.mraf.n_iter
    U0 = input_amplitude .* exp(1i * phase);
    Uf = mraf_forward_fft(U0);
    Af = abs(Uf);
    if isfield(cfg.mraf, 'use_wgs') && cfg.mraf.use_wgs && strcmpi(string(cfg.mraf.method), "WGS-MRAF")
        weights = update_wgs_weights(weights, target.amplitude, Af, finite_mask, cfg, target.masks.wgs_update_scale);
    end
    Uf_proj = project_farfield_mraf(Uf, target.amplitude, target.masks, weights, cfg);
    U_back = mraf_backward_fft(Uf_proj);
    phase_new = angle(U_back);
    dphase = angle(exp(1i * (phase_new - phase)));
    phase = phase + cfg.mraf.phase_blend .* dphase;
    phase = angle(exp(1i * phase));

    Uf_iter = mraf_forward_fft(input_amplitude .* exp(1i * phase));
    iter_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, abs(Uf_iter).^2, target.masks, input_power, cfg);
    iteration(iter, :) = [iter, iter_metrics.size_50_x_um, iter_metrics.size_50_y_um, ...
        iter_metrics.transition_13p5_90_x_um, iter_metrics.transition_13p5_90_y_um, ...
        iter_metrics.size13p5_x_um, iter_metrics.size13p5_y_um, iter_metrics.core_rms, iter_metrics.shoulder_peak_x, iter_metrics.shoulder_peak_y, iter_metrics.side_lobe_peak_x_rel_to_core, iter_metrics.side_lobe_peak_y_rel_to_core];
end

U_mraf = mraf_forward_fft(input_amplitude .* exp(1i * phase));
I_mraf_raw = abs(U_mraf).^2;
mraf_metrics = compute_mraf_metrics(focal_x_m, focal_y_m, I_mraf_raw, target.masks, input_power, cfg);

result = struct();
result.target = target;
result.rd.U = U_rd;
result.rd.intensity_raw = I_rd_raw;
result.rd.intensity_norm = I_rd_raw ./ max(mean(I_rd_raw(target.masks.core), 'omitnan'), eps);
result.rd.phase = phase0;
result.rd.phase_wrapped = wrap_mraf_phase(phase0);
result.rd.metrics = rd_metrics;
result.mraf.U = U_mraf;
result.mraf.intensity_raw = I_mraf_raw;
result.mraf.intensity_norm = I_mraf_raw ./ max(mean(I_mraf_raw(target.masks.core), 'omitnan'), eps);
result.mraf.phase = phase;
result.mraf.phase_wrapped = wrap_mraf_phase(phase);
result.mraf.metrics = mraf_metrics;
result.iteration_metrics = array2table(iteration, 'VariableNames', {'iter','size50_x_um','size50_y_um','transition_13p5_90_x_um','transition_13p5_90_y_um','size13p5_x_um','size13p5_y_um','core_rms','shoulder_peak_x','shoulder_peak_y','side_lobe_x','side_lobe_y'});
end
