function weights = update_wgs_weights(weights, A_target, Af, finite_mask, cfg, update_scale)
% update_wgs_weights Conservative WGS amplitude weights on finite target pixels.
if nargin < 6 || isempty(update_scale)
    update_scale = ones(size(Af));
end
ratio = ones(size(Af));
ratio(finite_mask) = A_target(finite_mask) ./ max(Af(finite_mask), eps);
ratio(~isfinite(ratio)) = 1;
local_exponent = cfg.mraf.wgs_exponent .* update_scale;
weights(finite_mask) = weights(finite_mask) .* ratio(finite_mask) .^ local_exponent(finite_mask);
mean_weight = mean(weights(finite_mask), 'omitnan');
if isfinite(mean_weight) && mean_weight > 0
    weights(finite_mask) = weights(finite_mask) ./ mean_weight;
end
weights = min(max(weights, 0.1), 10);
end
