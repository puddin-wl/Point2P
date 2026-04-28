function [cfg, result, fullprop] = run_flat_core_caseE_guard_lobe()
% run_flat_core_caseE_guard_lobe Four-region target with shoulder guard and lobe reservoir.
close all; clc;
project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));

variant = struct('variant_name', "CaseE_guard_lobe", ...
    'preset_name', "flat_core_guard_lobe", 'artifact_group', "flat_core_caseE_guard_lobe", ...
    'target_mode', "flat_core_guard_lobe", 'method', "WGS-MRAF", 'use_wgs', true, ...
    'use_three_region_projection', true, ...
    'flat_fraction_x', 0.90, 'flat_fraction_y', 0.90, ...
    'guard_fraction_x', 1.02, 'guard_fraction_y', 1.04, ...
    'lobe_fraction_x', 1.25, 'lobe_fraction_y', 1.35, ...
    'free_factor', 1.00, 'outer_factor', 0.80, ...
    'wgs_exponent', 0.08, 'n_iter', 20, 'phase_blend', 0.20, ...
    'guard_cap_outer_intensity', 0.12, 'guard_blend', 0.35);

[cfg, result, fullprop] = run_mraf_one(variant);
print_caseE_summary(cfg, result, fullprop);
end

function print_caseE_summary(cfg, result, fullprop)
m = result.mraf.metrics; rd = result.rd.metrics;
fprintf('\nCaseE_guard_lobe summary\n');
fprintf('Artifact: %s\n', cfg.save_root);
fprintf('flat_core_rms %.6g (RD %.6g, improvement %.6g), pv %.6g\n', m.flat_core_rms, rd.flat_core_rms, rd.flat_core_rms - m.flat_core_rms, m.flat_core_pv);
fprintf('shoulder_peak %.6g, shoulder_energy_ratio %.6g\n', m.shoulder_peak, m.shoulder_energy_ratio);
fprintf('lobe_peak %.6g, lobe_energy_ratio %.6g\n', m.lobe_peak, m.lobe_energy_ratio);
fprintf('side_lobe_energy %.6g, outer_energy %.6g\n', m.side_lobe_energy_ratio, m.outer_energy_ratio);
fprintf('size50 %.3f x %.3f um, transition %.3f x %.3f um\n', m.size50_x_um, m.size50_y_um, m.transition_width_x_um, m.transition_width_y_um);
fprintf('legacy shoulder %.6g x %.6g, outer_tail %.6g x %.6g\n', m.shoulder_peak_x, m.shoulder_peak_y, m.outer_tail_peak_x_rel_to_core, m.outer_tail_peak_y_rel_to_core);
fm = fullprop.mraf_metrics;
fprintf('FULL flat_core_rms %.6g, shoulder_peak %.6g, lobe_peak %.6g\n', fm.flat_core_rms, fm.shoulder_peak, fm.lobe_peak);
end
