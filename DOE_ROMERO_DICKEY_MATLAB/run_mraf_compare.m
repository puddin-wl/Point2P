function artifact_roots = run_mraf_compare()
% run_mraf_compare Run three small conservative smooth-target MRAF cases.
cases = cell(3, 1);

cases{1} = struct('artifact_group', "mraf_rd_derived_initial", 'variant_name', "mraf_smooth_gamma10_conservative", ...
    'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 0.7, ...
    'phase_blend', 0.25, 'method', "GS-MRAF");

cases{2} = struct('artifact_group', "mraf_rd_derived_initial", 'variant_name', "mraf_smooth_gamma11_conservative", ...
    'n_iter', 10, 'rd_target_gamma', 1.1, 'mraf_factor', 0.7, ...
    'phase_blend', 0.25, 'method', "GS-MRAF");

cases{3} = struct('artifact_group', "mraf_rd_derived_initial", 'variant_name', "mraf_smooth_gamma10_freekeep", ...
    'n_iter', 10, 'rd_target_gamma', 1.0, 'mraf_factor', 1.0, ...
    'phase_blend', 0.25, 'method', "GS-MRAF");

artifact_roots = strings(numel(cases), 1);
for k = 1:numel(cases)
    cfg = run_mraf_one(cases{k});
    artifact_roots(k) = string(cfg.save_root);
end

fprintf('\nCompleted MRAF comparison cases:\n');
for k = 1:numel(artifact_roots)
    fprintf('%d) %s\n', k, artifact_roots(k));
end
end