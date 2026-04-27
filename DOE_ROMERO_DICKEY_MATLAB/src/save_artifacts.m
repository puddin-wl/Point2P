function save_artifacts(cfg, grid, phase_info, focal, target_info, metrics, ...
    x_m, y_m, focal_x_m, focal_y_m, input_amplitude, input_intensity, ...
    aperture_mask, phase_unwrapped_rad, phase_wrapped_rad, focal_intensity, target_intensity, prefix)
% save_artifacts 保存配置、相位、指标、剖面 CSV 和文本报告。

if nargin < 18
    prefix = '';
end
if ~exist(cfg.save_root, 'dir')
    mkdir(cfg.save_root);
end

focal_summary = rmfield(focal, {'U', 'intensity_raw', 'intensity_norm'});
save(fullfile(cfg.save_root, [prefix 'config_used.mat']), 'cfg', 'grid', 'phase_info', 'focal_summary', 'target_info');
save(fullfile(cfg.save_root, [prefix 'phase_unwrapped.mat']), 'phase_unwrapped_rad', '-v7.3');
save(fullfile(cfg.save_root, [prefix 'phase_wrapped.mat']), 'phase_wrapped_rad', '-v7.3');

write_struct_text(fullfile(cfg.save_root, [prefix 'config_used.txt']), cfg);
write_beta_report(fullfile(cfg.save_root, [prefix 'beta_report.txt']), cfg);
write_metrics_report(fullfile(cfg.save_root, [prefix 'metrics.txt']), metrics);

x_table = table(focal_x_m(:), metrics.x_profile(:), metrics.target_x_profile(:), ...
    'VariableNames', {'x_m', 'intensity_norm', 'target'});
y_table = table(focal_y_m(:), metrics.y_profile(:), metrics.target_y_profile(:), ...
    'VariableNames', {'y_m', 'intensity_norm', 'target'});
writetable(x_table, fullfile(cfg.save_root, [prefix 'x_profile.csv']));
writetable(y_table, fullfile(cfg.save_root, [prefix 'y_profile.csv']));

save(fullfile(cfg.save_root, [prefix 'workspace_summary.mat']), 'x_m', 'y_m', 'focal_x_m', 'focal_y_m', ...
    'input_amplitude', 'input_intensity', 'aperture_mask', 'focal_intensity', 'target_intensity', 'metrics', '-v7.3');
end

function write_struct_text(path, s)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Config used for Romero-Dickey DOE run\n');
fprintf(fid, 'Generated: %s\n\n', char(datetime('now')));
fields = fieldnames(s);
for k = 1:numel(fields)
    name = fields{k};
    value = s.(name);
    fprintf(fid, '%s = %s\n', name, value_to_string(value));
end
end

function write_beta_report(path, cfg)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Romero-Dickey beta report\n');
fprintf(fid, '==========================\n\n');
fprintf(fid, 'lambda = %.9g m\n', cfg.lambda_m);
fprintf(fid, 'f = %.9g m\n', cfg.f_m);
fprintf(fid, 'input 1/e^2 intensity diameter = %.9g m\n', cfg.input_1e2_diameter_m);
fprintf(fid, 'input 1/e^2 intensity radius w = %.9g m\n', cfg.input_1e2_radius_m);
fprintf(fid, 'paper ri = input_1e2_radius_m / sqrt(2) = %.9g m\n', cfg.input_1e_radius_m);
fprintf(fid, 'target size x = %.9g m, y = %.9g m\n', cfg.target_size_x_m, cfg.target_size_y_m);
fprintf(fid, 'Ro definition = %s\n', cfg.Ro_definition);
fprintf(fid, 'Ro_x = target_size_x_m / sqrt(pi) = %.9g m\n', cfg.Ro_x_m);
fprintf(fid, 'Ro_y = target_size_y_m / sqrt(pi) = %.9g m\n\n', cfg.Ro_y_m);
fprintf(fid, 'beta = 2*pi*ri*Ro/(lambda*f)\n');
fprintf(fid, 'beta_x = %.9f\n', cfg.beta_x);
fprintf(fid, 'beta_y = %.9f\n\n', cfg.beta_y);
fprintf(fid, 'Interpretation:\n');
fprintf(fid, '- beta is the key dimensionless difficulty parameter in the Romero-Dickey stationary-phase approximation.\n');
fprintf(fid, '- Larger beta usually gives a more geometrical-optics-like mapping and sharper flat-top behavior.\n');
fprintf(fid, '- beta_y is smaller because the requested y size is only 120 um, so y-edge sharpness and uniformity are expected to be harder.\n');
fprintf(fid, '- The 15 mm aperture is a clear aperture mask, not the 5 mm Gaussian beam diameter.\n');
end

function write_metrics_report(path, metrics)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Focal-plane metrics\n');
fprintf(fid, '===================\n\n');
fprintf(fid, 'Normalization: %s\n\n', metrics.normalization);
fprintf(fid, 'FWHM_x = %.6f um\n', metrics.fwhm_x_m * 1e6);
fprintf(fid, 'FWHM_y = %.6f um\n', metrics.fwhm_y_m * 1e6);
fprintf(fid, 'Transition 13-90 x = %.6f um\n', metrics.transition_width_13_90_x_m * 1e6);
fprintf(fid, 'Transition 13-90 y = %.6f um\n\n', metrics.transition_width_13_90_y_m * 1e6);
fprintf(fid, 'Core mean = %.9g\n', metrics.core_mean);
fprintf(fid, 'Core RMS relative = %.9g\n', metrics.core_rms);
fprintf(fid, 'Peak-to-valley relative = %.9g\n', metrics.peak_to_valley);
fprintf(fid, 'Center profile std x relative = %.9g\n', metrics.center_profile_std_x);
fprintf(fid, 'Center profile std y relative = %.9g\n\n', metrics.center_profile_std_y);
fprintf(fid, 'Efficiency target ROI = %.9g\n', metrics.efficiency_target);
fprintf(fid, 'Efficiency guard ROI = %.9g\n\n', metrics.efficiency_guard);
fprintf(fid, 'Side-lobe peak x relative to core = %.9g\n', metrics.side_lobe_peak_x_rel_to_core);
fprintf(fid, 'Side-lobe peak y relative to core = %.9g\n', metrics.side_lobe_peak_y_rel_to_core);
fprintf(fid, 'Overshoot near edge x = %.9g\n', metrics.overshoot_near_edge_x);
fprintf(fid, 'Overshoot near edge y = %.9g\n', metrics.overshoot_near_edge_y);
fprintf(fid, 'Undershoot inside x = %.9g\n', metrics.undershoot_inside_x);
fprintf(fid, 'Undershoot inside y = %.9g\n', metrics.undershoot_inside_y);
end

function text = value_to_string(value)
if isstring(value) || ischar(value)
    text = char(value);
elseif isnumeric(value) || islogical(value)
    if isscalar(value)
        text = sprintf('%.12g', value);
    else
        text = mat2str(value);
    end
else
    text = sprintf('<%s>', class(value));
end
end
