function weights = update_wgs_weights(weights, A_target, Af, finite_mask, cfg)
% update_wgs_weights Update WGS weights only on the finite flat-core target.
core_amp = Af(finite_mask);
mean_amp = mean(core_amp(:), 'omitnan') + eps;
ratio = mean_amp ./ (core_amp + eps);
core_weights = weights(finite_mask);
core_weights(~isfinite(core_weights)) = 1;
core_weights = core_weights .* ratio .^ cfg.mraf.wgs_exponent;
core_weights = min(max(core_weights, 0.5), 2.0);
mean_weight = mean(core_weights(:), 'omitnan');
if isfinite(mean_weight) && mean_weight > 0
    core_weights = core_weights ./ mean_weight;
end
weights(finite_mask) = core_weights;
weights(~finite_mask) = NaN;
end
