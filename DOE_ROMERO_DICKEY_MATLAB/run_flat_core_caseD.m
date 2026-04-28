function [cfg, result, fullprop] = run_flat_core_caseD()
% run_flat_core_caseD User-requested flat-core/free-edge CaseD.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

variant = struct('variant_name', "caseD_flat094_free108112_outer090", ...
    'preset_name', "flat_core_caseD", 'artifact_group', "flat_core_caseD", ...
    'target_mode', "flat_core_free_edge", 'method', "WGS-MRAF", 'use_wgs', true, ...
    'use_three_region_projection', true, ...
    'flat_fraction_x', 0.94, 'flat_fraction_y', 0.94, ...
    'free_fraction_x', 1.08, 'free_fraction_y', 1.12, ...
    'wgs_exponent', 0.10, 'n_iter', 20, 'phase_blend', 0.20, ...
    'free_factor', 1.00, 'outer_factor', 0.90);

[cfg, result, fullprop] = run_mraf_one(variant);
print_caseD_summary(cfg, result, fullprop);
end

function print_caseD_summary(cfg, result, fullprop)
m = result.mraf.metrics;
rd = result.rd.metrics;
fprintf('\nCaseD summary\n');
fprintf('Artifact: %s\n', cfg.save_root);
fprintf('flat_core_rms %.6g (RD %.6g, improvement %.6g)\n', m.flat_core_rms, rd.flat_core_rms, rd.flat_core_rms - m.flat_core_rms);
fprintf('flat_core_pv %.6g, uniformity %.6g\n', m.flat_core_pv, m.flat_core_uniformity);
fprintf('side_lobe_energy %.6g, outer_energy %.6g\n', m.side_lobe_energy_ratio, m.outer_energy_ratio);
fprintf('size50 %.3f x %.3f um, transition %.3f x %.3f um\n', m.size50_x_um, m.size50_y_um, m.transition_width_x_um, m.transition_width_y_um);
fprintf('shoulder %.6g x %.6g, outer_tail %.6g x %.6g\n', m.shoulder_peak_x, m.shoulder_peak_y, m.outer_tail_peak_x_rel_to_core, m.outer_tail_peak_y_rel_to_core);
fm = fullprop.mraf_metrics;
fprintf('FULL flat_core_rms %.6g, size50 %.3f x %.3f um\n', fm.flat_core_rms, fm.size50_x_um, fm.size50_y_um);
end
