function Uf_proj = project_farfield_mraf(Uf, A_target, masks, weights, cfg)
% project_farfield_mraf Apply flat-core/free-edge/outer-suppress MRAF projection.
current_phase = exp(1i * angle(Uf));
current_amp = abs(Uf);
finite_mask = masks.finite_target & isfinite(A_target);
if isfield(masks, 'wgs')
    finite_mask = masks.wgs & isfinite(A_target) & A_target > 0;
end

if isfield(cfg.mraf, 'use_three_region_projection') && cfg.mraf.use_three_region_projection
    projected_amp = current_amp;
    if isfield(masks, 'null') && isfield(masks, 'noise') && isfield(masks, 'signal')
        signal_mask = masks.signal & isfinite(A_target) & A_target > 0;
        noise_mask = masks.noise & isnan(A_target);
        null_mask = masks.null & isfinite(A_target) & A_target == 0;

        projected_amp(:) = 0;
        projected_amp(noise_mask) = cfg.mraf.mraf_factor .* current_amp(noise_mask);
        projected_amp(null_mask) = 0;

        target_amp = A_target(signal_mask);
        if nargin >= 4 && ~isempty(weights)
            signal_weights = weights(signal_mask);
            signal_weights(~isfinite(signal_weights)) = 1;
            target_amp = target_amp .* signal_weights;
        end
        projected_amp(signal_mask) = target_amp;
        Uf_proj = projected_amp .* current_phase;
        Uf_proj(~isfinite(Uf_proj)) = 0;
        return;
    end
    if isfield(masks, 'plateau') && isfield(masks, 'transition') && isfield(masks, 'outside_e2')
        plateau_mask = masks.plateau & finite_mask;
        transition_mask = masks.transition & isfinite(A_target);
        outside_e2_mask = masks.outside_e2;

        cap_blend = get_mraf_field(cfg, 'transition_cap_blend', 0.35);
        cap_amp = A_target;
        over_cap = transition_mask & current_amp > cap_amp;
        projected_amp(over_cap) = (1 - cap_blend) .* current_amp(over_cap) + cap_blend .* cap_amp(over_cap);
        projected_amp(outside_e2_mask) = current_amp(outside_e2_mask) .* cfg.mraf.outer_factor;

        target_amp = ones(nnz(plateau_mask), 1);
        if nargin >= 4 && ~isempty(weights)
            core_weights = weights(plateau_mask);
            core_weights(~isfinite(core_weights)) = 1;
            target_amp = target_amp .* core_weights(:);
        end
        projected_amp(plateau_mask) = target_amp;
        Uf_proj = projected_amp .* current_phase;
        Uf_proj(~isfinite(Uf_proj)) = 0;
        return;
    end
    if isfield(masks, 'shoulder_guard') && isfield(cfg.mraf, 'guard_cap_outer_intensity')
        if isfield(masks, 'guard_cap_amplitude')
            cap_amp = masks.guard_cap_amplitude;
        else
            cap_amp = shoulder_guard_cap_amplitude(size(current_amp), masks.shoulder_guard, cfg);
        end
        over_cap = masks.shoulder_guard & current_amp > cap_amp;
        capped_amp = (1 - cfg.mraf.guard_blend) .* current_amp(over_cap) + cfg.mraf.guard_blend .* cap_amp(over_cap);
        projected_amp(over_cap) = capped_amp;
    end
    if isfield(masks, 'free_edge')
        projected_amp(masks.free_edge) = current_amp(masks.free_edge) .* cfg.mraf.free_factor;
    elseif isfield(masks, 'free_noise')
        projected_amp(masks.free_noise) = current_amp(masks.free_noise) .* cfg.mraf.free_factor;
    end
    if isfield(masks, 'outer_suppress')
        projected_amp(masks.outer_suppress) = current_amp(masks.outer_suppress) .* cfg.mraf.outer_factor;
    end
else
    projected_amp = current_amp .* cfg.mraf.mraf_factor;
end

target_amp = A_target(finite_mask);
if nargin >= 4 && ~isempty(weights)
    core_weights = weights(finite_mask);
    core_weights(~isfinite(core_weights)) = 1;
    target_amp = target_amp .* core_weights;
end
projected_amp(finite_mask) = target_amp;
Uf_proj = projected_amp .* current_phase;
Uf_proj(~isfinite(Uf_proj)) = 0;
end

function value = get_mraf_field(cfg, name, default_value)
if isfield(cfg, 'mraf') && isfield(cfg.mraf, name)
    value = cfg.mraf.(name);
else
    value = default_value;
end
end

function cap_amp = shoulder_guard_cap_amplitude(array_size, shoulder_guard, cfg)
% Reconstruct a conservative cap amplitude map using mask distance surrogate.
% The exact target builder stores the physical cap map in target, but this
% projection API receives masks only, so use a smooth index-space distance to
% avoid adding more global state. The cap is only applied where current_amp
% exceeds it; lower amplitudes are never raised.
cap_amp = inf(array_size);
if ~any(shoulder_guard, 'all')
    return;
end
[rows, cols] = find(shoulder_guard);
row_min = min(rows); row_max = max(rows); col_min = min(cols); col_max = max(cols);
[CC, RR] = meshgrid(1:array_size(2), 1:array_size(1));
dx = min(abs(CC - col_min), abs(CC - col_max)) ./ max((col_max - col_min) / 2, eps);
dy = min(abs(RR - row_min), abs(RR - row_max)) ./ max((row_max - row_min) / 2, eps);
t = 1 - max(dx, dy);
t = min(max(t, 0), 1);
smooth_t = t.^2 .* (3 - 2 .* t);
cap_intensity = cfg.mraf.guard_cap_outer_intensity + (1 - cfg.mraf.guard_cap_outer_intensity) .* smooth_t;
cap_amp(shoulder_guard) = sqrt(max(cap_intensity(shoulder_guard), 0));
end
