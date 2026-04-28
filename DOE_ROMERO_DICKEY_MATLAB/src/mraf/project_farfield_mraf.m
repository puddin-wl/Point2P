function Uf_proj = project_farfield_mraf(Uf, A_target, masks, weights, cfg)
% project_farfield_mraf Apply finite target amplitude and MRAF free-region projection.
Af_phase = exp(1i * angle(Uf));
finite_mask = masks.finite_target & isfinite(A_target);
if isfield(cfg.mraf, 'use_three_region_projection') && cfg.mraf.use_three_region_projection
    free_factor = cfg.mraf.free_factor;
    outer_factor = cfg.mraf.outer_factor;
    Uf_proj = Uf;
    if isfield(masks, 'free_noise')
        Uf_proj(masks.free_noise) = free_factor .* Uf(masks.free_noise);
    end
    if isfield(masks, 'outer_suppress')
        Uf_proj(masks.outer_suppress) = outer_factor .* Uf(masks.outer_suppress);
    end
else
    Uf_proj = cfg.mraf.mraf_factor .* Uf;
end
projected_amp = A_target;
if nargin >= 4 && ~isempty(weights)
    projected_amp(finite_mask) = projected_amp(finite_mask) .* weights(finite_mask);
end
Uf_proj(finite_mask) = projected_amp(finite_mask) .* Af_phase(finite_mask);
Uf_proj(~isfinite(Uf_proj)) = 0;
end
