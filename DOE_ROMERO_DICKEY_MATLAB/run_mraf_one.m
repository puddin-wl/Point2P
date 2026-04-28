function [cfg, result, fullprop] = run_mraf_one(variant)
% run_mraf_one Run RD-derived MRAF refinement after Romero-Dickey baseline.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   This is the main MRAF entry that generates phase0 before any MRAF target or
%   projection is applied.
% Inputs:
%   default_config.m supplies the DOE grid, wavelength, focal length, 15 mm clear
%   aperture, 5 mm incident Gaussian beam, and 330 x 120 um target size.
% Outputs:
%   phase_rd_unwrapped is passed as phase0 to run_mraf_refinement.m. The wrapped
%   version is saved as common/phase_rd_wrapped.mat and visualized as PNG.
% Physical meaning:
%   phase0 is the analytical separable Romero-Dickey initial DOE phase. MRAF
%   starts from it and may change the final phase.
% Used by:
%   All run_mraf_* probe scripts that call run_mraf_one.m.
% Notes:
%   The target_mode selected for MRAF is not used to compute phase0; phase0 is
%   always generated before target construction.
if nargin < 1
    variant = struct();
end
close all; clc;

project_root = fileparts(mfilename('fullpath'));
addpath(fullfile(project_root, 'config'));
addpath(fullfile(project_root, 'src'));
addpath(fullfile(project_root, 'src', 'mraf'));
addpath(fullfile(fileparts(project_root), 'initial_phase_generation'));

cfg = default_config(project_root);
cfg = apply_variant(cfg, variant);
stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
group_root = fullfile(project_root, 'artifacts', char(string(cfg.mraf.artifact_group)));
if ~exist(group_root, 'dir'), mkdir(group_root); end
cfg.save_root = fullfile(group_root, char(stamp + "_" + string(cfg.mraf.variant_name)));
cfg = prepare_artifact_dirs(cfg);

fprintf('\n=== Point2P RD-derived MRAF refinement ===\n');
fprintf('Variant: %s\n', string(cfg.mraf.variant_name));
fprintf('N = %d, focal dx = %.3f um/pixel\n', cfg.N, cfg.focal_dx_m * 1e6);
fprintf('MRAF: iter=%d, factor=%.3f, phase_blend=%.3f, gamma=%.3f, method=%s\n', ...
    cfg.mraf.n_iter, cfg.mraf.mraf_factor, cfg.mraf.phase_blend, cfg.mraf.rd_target_gamma, string(cfg.mraf.method));
fprintf('Save root: %s\n\n', cfg.save_root);

[X_m, Y_m, x_m, y_m, grid] = make_grid(cfg); %#ok<ASGLU>
phase_data = generate_initial_phase(cfg);
input_amplitude = phase_data.input_amplitude;
phase_rd_unwrapped = phase_data.phase0_unwrapped_rad;
phase_rd_wrapped = phase_data.phase0_wrapped_rad;
phase_info = phase_data.phase_info; %#ok<NASGU>

focal_index = (-cfg.N/2):(cfg.N/2 - 1);
focal_x_m = focal_index * cfg.focal_dx_m;
focal_y_m = focal_x_m;

fprintf('Step 1/4: ideal FFT RD baseline and smooth RD-derived target...\n');
result = run_mraf_refinement(input_amplitude, phase_rd_unwrapped, focal_x_m, focal_y_m, cfg);

fprintf('Step 2/4: full propagation verification...\n');
U_rd_doe = input_amplitude .* exp(1i * phase_rd_unwrapped);
U_mraf_doe = input_amplitude .* exp(1i * result.mraf.phase);
[full_x_m, full_y_m, ~, full_rd_focal] = full_lens_propagate(U_rd_doe, X_m, Y_m, cfg.dx_doe_m, cfg);
[~, ~, ~, full_mraf_focal] = full_lens_propagate(U_mraf_doe, X_m, Y_m, cfg.dx_doe_m, cfg);
full_masks = result.target.masks;
input_power = sum(abs(input_amplitude).^2, 'all');
full_rd_metrics = compute_mraf_metrics(full_x_m, full_y_m, full_rd_focal.intensity_raw, full_masks, input_power, cfg);
full_mraf_metrics = compute_mraf_metrics(full_x_m, full_y_m, full_mraf_focal.intensity_raw, full_masks, input_power, cfg);
fullprop = struct();
fullprop.cfg = cfg;
fullprop.focal_x_m = full_x_m;
fullprop.focal_y_m = full_y_m;
fullprop.rd_metrics = full_rd_metrics;
fullprop.mraf_metrics = full_mraf_metrics;
fullprop.rd_intensity_norm = full_rd_focal.intensity_raw ./ max(mean(full_rd_focal.intensity_raw(full_masks.core), 'omitnan'), eps);
fullprop.mraf_intensity_norm = full_mraf_focal.intensity_raw ./ max(mean(full_mraf_focal.intensity_raw(full_masks.core), 'omitnan'), eps);
fullprop.model = "DOE -> angular spectrum 200 mm -> lens pupil/quadratic phase -> single-FFT focal propagation 429 mm";

fprintf('Step 3/4: saving artifacts and reports...\n');
save(fullfile(cfg.artifact_dirs.common, 'phase_rd_wrapped.mat'), 'phase_rd_wrapped', '-v7.3');
phase_mraf_wrapped = result.mraf.phase_wrapped;
save(fullfile(cfg.artifact_dirs.common, 'phase_mraf_wrapped.mat'), 'phase_mraf_wrapped', '-v7.3');
writetable(result.iteration_metrics, fullfile(cfg.artifact_dirs.common, 'iteration_metrics.csv'));
write_config_text(fullfile(cfg.artifact_dirs.common, 'config_used.txt'), cfg);
write_config_json(fullfile(cfg.artifact_dirs.common, 'config_used.json'), cfg);
write_metrics_text(fullfile(cfg.artifact_dirs.ideal_fft, 'metrics_rd.txt'), result.rd.metrics, 'RD baseline ideal FFT metrics');
write_metrics_text(fullfile(cfg.artifact_dirs.ideal_fft, 'metrics_mraf.txt'), result.mraf.metrics, 'MRAF ideal FFT metrics');
write_compare_text(fullfile(cfg.artifact_dirs.ideal_fft, 'metrics_compare.txt'), result.rd.metrics, result.mraf.metrics, 'Ideal FFT RD vs MRAF metrics');
write_metrics_text(fullfile(cfg.artifact_dirs.full_propagation, 'fullprop_metrics_rd.txt'), full_rd_metrics, 'RD full propagation metrics');
write_metrics_text(fullfile(cfg.artifact_dirs.full_propagation, 'fullprop_metrics_mraf.txt'), full_mraf_metrics, 'MRAF full propagation metrics');
write_compare_text(fullfile(cfg.artifact_dirs.full_propagation, 'fullprop_metrics_compare.txt'), full_rd_metrics, full_mraf_metrics, 'Full propagation RD vs MRAF metrics');
write_mraf_report(fullfile(cfg.artifact_dirs.common, 'mraf_report.txt'), cfg, result.rd.metrics, result.mraf.metrics, full_rd_metrics, full_mraf_metrics, fullprop.model, result.target, focal_x_m, focal_y_m);
write_mask_report(fullfile(cfg.artifact_dirs.target_masks, 'mask_report.txt'), cfg, result.target, focal_x_m, focal_y_m);

center_x = round(numel(focal_x_m)/2)+1; center_y = round(numel(focal_y_m)/2)+1;
x_profile_table = table(focal_x_m(:), result.rd.metrics.x_profile(:), result.target.intensity(center_y,:).', result.mraf.metrics.x_profile(:), ...
    'VariableNames', {'x_m','rd_intensity_core_norm','target_intensity','mraf_intensity_core_norm'});
y_profile_table = table(focal_y_m(:), result.rd.metrics.y_profile(:), result.target.intensity(:,center_x), result.mraf.metrics.y_profile(:), ...
    'VariableNames', {'y_m','rd_intensity_core_norm','target_intensity','mraf_intensity_core_norm'});
writetable(x_profile_table, fullfile(cfg.artifact_dirs.ideal_fft, 'x_profile_rd_target_mraf.csv'));
writetable(y_profile_table, fullfile(cfg.artifact_dirs.ideal_fft, 'y_profile_rd_target_mraf.csv'));

fprintf('Step 4/4: plotting diagnostics...\n');
plot_mraf_diagnostics(cfg, focal_x_m, focal_y_m, result.rd, result.target, result.mraf, fullprop);

fprintf('\nDone. MRAF artifacts written to:\n%s\n', cfg.save_root);
fprintf('Ideal shoulder x/y RD: %.6g / %.6g, MRAF: %.6g / %.6g\n', result.rd.metrics.shoulder_peak_x, result.rd.metrics.shoulder_peak_y, result.mraf.metrics.shoulder_peak_x, result.mraf.metrics.shoulder_peak_y);
fprintf('Ideal core RMS RD: %.6g, MRAF: %.6g\n', result.rd.metrics.core_rms, result.mraf.metrics.core_rms);
end

function cfg = apply_variant(cfg, variant)
if ischar(variant) || isstring(variant)
    cfg.mraf.variant_name = string(variant);
    cfg.mraf.preset_name = string(variant);
    return;
end
if ~isstruct(variant), return; end
fields = fieldnames(variant);
for k = 1:numel(fields)
    name = fields{k};
    if strcmp(name, 'variant_name')
        cfg.mraf.variant_name = string(variant.(name));
    elseif isfield(cfg.mraf, name)
        cfg.mraf.(name) = variant.(name);
    else
        cfg.(name) = variant.(name);
    end
end
end

function cfg = prepare_artifact_dirs(cfg)
if ~exist(cfg.save_root, 'dir'), mkdir(cfg.save_root); end
cfg.artifact_dirs = struct();
cfg.artifact_dirs.common = fullfile(cfg.save_root, 'common');
cfg.artifact_dirs.target_masks = fullfile(cfg.save_root, 'target_masks');
cfg.artifact_dirs.ideal_fft = fullfile(cfg.save_root, 'ideal_fft');
cfg.artifact_dirs.full_propagation = fullfile(cfg.save_root, 'full_propagation');
names = fieldnames(cfg.artifact_dirs);
for k = 1:numel(names)
    if ~exist(cfg.artifact_dirs.(names{k}), 'dir'), mkdir(cfg.artifact_dirs.(names{k})); end
end
end

function write_config_json(path, cfg)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
try
    fprintf(fid, '%s', jsonencode(cfg, 'PrettyPrint', true));
catch
    fprintf(fid, '%s', jsonencode(cfg));
end
end

function write_config_text(path, cfg)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Point2P MRAF config\nGenerated: %s\n\n', char(datetime('now')));
write_struct(fid, cfg, '');
end

function write_struct(fid, s, prefix)
fields = fieldnames(s);
for k = 1:numel(fields)
    name = fields{k}; value = s.(name); key = [prefix name];
    if isstruct(value)
        write_struct(fid, value, [key '.']);
    else
        fprintf(fid, '%s = %s\n', key, value_to_string(value));
    end
end
end

function write_metrics_text(path, metrics, title_text)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n%s\n\n', title_text, repmat('=', 1, strlength(title_text)));
write_metric_fields(fid, metrics);
end

function write_compare_text(path, rd, mraf, title_text)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n%s\n\n', title_text, repmat('=', 1, strlength(title_text)));
fprintf(fid, 'Region definitions are recorded in common/mraf_report.txt and target_masks/mask_report.txt.\n\n');
names = metric_names();
fprintf(fid, '%-38s %16s %16s %16s\n', 'metric', 'RD', 'MRAF', 'delta');
for k = 1:numel(names)
    name = names{k}; a = rd.(name); b = mraf.(name);
    fprintf(fid, '%-38s %16.8g %16.8g %16.8g\n', name, a, b, b-a);
end
end

function write_mraf_report(path, cfg, rd, mraf, full_rd, full_mraf, full_model, target, focal_x_m, focal_y_m)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Point2P RD-derived MRAF report\n===============================\n\n');
fprintf(fid, 'Preset: %s\n', string(cfg.mraf.preset_name));
fprintf(fid, 'Variant: %s\n', string(cfg.mraf.variant_name));
fprintf(fid, 'Initial phase: Romero-Dickey / Point2P separable phase.\n');
fprintf(fid, 'Target mode: %s; finite target is amplitude, free/noise region is NaN.\n', cfg.mraf.target_mode);
if strcmpi(string(cfg.mraf.target_mode), "flat_core_free_edge")
    fprintf(fid, 'Flat-core/free-edge target: only flat_core has target amplitude=1; free_edge is unconstrained; outer_suppress is weakly damped.\n');
    fprintf(fid, 'flat_fraction x/y %.3f / %.3f, free_fraction x/y %.3f / %.3f.\n', cfg.mraf.flat_fraction_x, cfg.mraf.flat_fraction_y, cfg.mraf.free_fraction_x, cfg.mraf.free_fraction_y);
elseif strcmpi(string(cfg.mraf.target_mode), "hard_rectangle")
    fprintf(fid, 'Hard rectangle target: signal |x|<=%.3f um, |y|<=%.3f um; outside is MRAF free region (NaN amplitude).\n', cfg.target_half_x_m*1e6, cfg.target_half_y_m*1e6);
else
    fprintf(fid, 'Smooth target: I_env = min(I_rd_norm,1)^gamma, I_flat = 1, I_target = blend*I_flat + (1-blend)*I_env.\n');
end
fprintf(fid, 'Smoothstep blend: flat_low=%.3f, flat_high=%.3f; t=(I-flat_low)/(flat_high-flat_low); blend=t^2*(3-2*t).\n', cfg.mraf.flat_low, cfg.mraf.flat_high);
fprintf(fid, 'Current baseline preset is caseC_stable: gamma=1.0, mraf_factor=1.0, phase_blend=0.25, n_iter=10, GS-MRAF, WGS off.\n');
fprintf(fid, 'Primary annotated metrics: ???????50%%?, ???????13.5%%?, ??????13.5%%-90%%?.\n');
fprintf(fid, 'Default goal is shoulder suppression first; transition narrowing is secondary.\n');
fprintf(fid, 'FFT convention: forward fftshift(fft2(fftshift(U)))/N, backward ifftshift(ifft2(ifftshift(Uf)))*N.\n');
fprintf(fid, 'Focal sampling dx_f = lambda*f/(N*dx_doe) = %.6f um/pixel.\n', cfg.focal_dx_m * 1e6);
fprintf(fid, 'DOE window = %.6f mm; 15 mm clear aperture is a circular mask inside this window.\n', cfg.doe_grid_extent_m * 1e3);
fprintf(fid, 'Full propagation model actually run: %s.\n\n', full_model);
fprintf(fid, 'Knobs: n_iter=%d, mraf_factor=%.3f, phase_blend=%.3f, gamma=%.3f, method=%s.\n', cfg.mraf.n_iter, cfg.mraf.mraf_factor, cfg.mraf.phase_blend, cfg.mraf.rd_target_gamma, string(cfg.mraf.method));
fprintf(fid, 'Projection mode: use_three_region=%d, free_factor=%.3f, outer_factor=%.3f. In three-region mode mraf_factor is legacy only.\n', cfg.mraf.use_three_region_projection, cfg.mraf.free_factor, cfg.mraf.outer_factor);
fprintf(fid, 'Target edge mode: %s; pivot50 inner_power=%.4f, outer_power=%.4f, remap(0.5)=%.6f.\n', string(cfg.mraf.target_edge_mode), cfg.mraf.pivot50_inner_power, cfg.mraf.pivot50_outer_power, target.remap_anchor_50);
fprintf(fid, 'Target geometry scale: x_scale=%.5f, y_scale=%.5f; target/masks use scaled RD geometry while RD baseline remains unscaled.\n', cfg.mraf.target_x_scale, cfg.mraf.target_y_scale);
fprintf(fid, 'Side-lobe trade knobs: y_free_reservoir=%d half_y=%.3f um max_I=%.3f; wgs_y_edge_damping=%d factor=%.3f band_I=[%.3f, %.3f].\n', ...
    cfg.mraf.y_free_reservoir_enable, cfg.mraf.y_free_reservoir_half_y_um, cfg.mraf.y_free_reservoir_max_intensity, ...
    cfg.mraf.wgs_y_edge_damping_enable, cfg.mraf.wgs_y_edge_damping_factor, cfg.mraf.wgs_y_edge_damping_low, cfg.mraf.wgs_y_edge_damping_high);
fprintf(fid, 'Region thresholds on scaled target geometry: core>=%.3f, signal>=%.3f, smooth blend %.3f..%.3f, free/noise<%.3f inside guard.\n', cfg.mraf.core_threshold, cfg.mraf.signal_threshold, cfg.mraf.flat_low, cfg.mraf.flat_high, cfg.mraf.free_threshold);
fprintf(fid, 'Three-region masks: target=guard & I>=%.3f; free_noise=guard & I<%.3f; outer_suppress=outside guard.\n', cfg.mraf.free_threshold, cfg.mraf.free_threshold);
fprintf(fid, 'Three-region projection: target uses target amplitude, free_noise keeps free_factor*Uf, outer_suppress keeps outer_factor*Uf.\n');
write_mask_summary(fid, target, focal_x_m, focal_y_m, cfg);
fprintf(fid, '\n');
fprintf(fid, 'Quick ideal summary:\n');
fprintf(fid, 'RD shoulder x/y %.8g / %.8g -> MRAF %.8g / %.8g\n', rd.shoulder_peak_x, rd.shoulder_peak_y, mraf.shoulder_peak_x, mraf.shoulder_peak_y);
fprintf(fid, 'RD core_rms %.8g -> MRAF %.8g\n', rd.core_rms, mraf.core_rms);
fprintf(fid, 'RD size50 x/y %.4f / %.4f um -> MRAF %.4f / %.4f um\n', rd.size50_x_um, rd.size50_y_um, mraf.size50_x_um, mraf.size50_y_um);
fprintf(fid, 'RD size13.5 x/y %.4f / %.4f um -> MRAF %.4f / %.4f um\n', rd.size13p5_x_um, rd.size13p5_y_um, mraf.size13p5_x_um, mraf.size13p5_y_um);
fprintf(fid, 'RD transition 13.5-90 x/y %.4f / %.4f um -> MRAF %.4f / %.4f um\n', rd.transition_13p5_90_x_um, rd.transition_13p5_90_y_um, mraf.transition_13p5_90_x_um, mraf.transition_13p5_90_y_um);
if rd.has_crossing_warning || mraf.has_crossing_warning, fprintf(fid, 'WARNING: one or more profile crossings were not found; NaN metrics are expected for those entries.\n'); end
fprintf(fid, '\nQuick full-prop summary:\n');
fprintf(fid, 'RD shoulder x/y %.8g / %.8g -> MRAF %.8g / %.8g\n', full_rd.shoulder_peak_x, full_rd.shoulder_peak_y, full_mraf.shoulder_peak_x, full_mraf.shoulder_peak_y);
fprintf(fid, 'RD core_rms %.8g -> MRAF %.8g\n', full_rd.core_rms, full_mraf.core_rms);
end


function write_mask_report(path, cfg, target, focal_x_m, focal_y_m)
fid = fopen(path, 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'MRAF mask report\n================\n\n');
if strcmpi(string(cfg.mraf.target_mode), "hard_rectangle")
    fprintf(fid, 'target_mode = hard_rectangle\n');
    fprintf(fid, 'signal rectangle size = %.3f x %.3f um\n', cfg.target_size_x_m*1e6, cfg.target_size_y_m*1e6);
    fprintf(fid, 'signal region bounds = |x|<=%.3f um, |y|<=%.3f um\n', cfg.target_half_x_m*1e6, cfg.target_half_y_m*1e6);
    fprintf(fid, 'outside signal rectangle is MRAF free region; target amplitude is NaN there\n');
    fprintf(fid, 'n_iter = %d, mraf_factor = %.3f, phase_blend = %.3f\n', cfg.mraf.n_iter, cfg.mraf.mraf_factor, cfg.mraf.phase_blend);
elseif strcmpi(string(cfg.mraf.target_mode), "rect_e2_spec_scaled")
    fprintf(fid, 'target_mode = rect_e2_spec_scaled\n');
    fprintf(fid, 'scaled industrial e^-2 rectangular target; outside 13.5%% is don''t-care with weak outer suppression\n');
    fprintf(fid, 'target 90%% size = %.3f x %.3f um\n', target.spec.size90_x_um, target.spec.size90_y_um);
    fprintf(fid, 'target 50%% size = %.3f x %.3f um\n', target.spec.size50_x_um, target.spec.size50_y_um);
    fprintf(fid, 'target 13.5%% size = %.3f x %.3f um\n', target.spec.size135_x_um, target.spec.size135_y_um);
    fprintf(fid, 'target transition 13.5%%-90%% = %.3f x %.3f um\n', target.spec.transition_width_135_90_x_um, target.spec.transition_width_135_90_y_um);
    fprintf(fid, 'projection: plateau strong WGS target, transition soft cap blend=%.3f, outside_e2 outer_factor=%.3f\n', cfg.mraf.transition_cap_blend, cfg.mraf.outer_factor);
elseif strcmpi(string(cfg.mraf.target_mode), "slmsuite_like_mraf")
    fprintf(fid, 'target_mode = slmsuite_like_mraf\n');
    fprintf(fid, 'target amplitude semantics: signal finite=1, noise=NaN, null/background=0\n');
    fprintf(fid, 'signal rectangle size = %.3f x %.3f um\n', target.spec.signal_size_x_um, target.spec.signal_size_y_um);
    fprintf(fid, 'noise rectangle size = %.3f x %.3f um\n', target.spec.noise_size_x_um, target.spec.noise_size_y_um);
    fprintf(fid, 'signal bounds = |x|<=%.3f um, |y|<=%.3f um\n', target.spec.signal_half_x_um, target.spec.signal_half_y_um);
    fprintf(fid, 'noise bounds = |x|<=%.3f um, |y|<=%.3f um\n', target.spec.noise_half_x_um, target.spec.noise_half_y_um);
    fprintf(fid, 'projection: signal uses WGS target, noise keeps mraf_factor * current amplitude, null is zero; mraf_factor=%.3f\n', cfg.mraf.mraf_factor);
else
    fprintf(fid, 'core_mask: I_rd_norm >= %.3f inside guard\n', cfg.mraf.core_threshold);
end
fprintf(fid, 'signal_mask: I_rd_norm >= %.3f inside guard\n', cfg.mraf.signal_threshold);
fprintf(fid, 'smooth blend: flat_low=%.3f, flat_high=%.3f, smoothstep t^2*(3-2*t)\n', cfg.mraf.flat_low, cfg.mraf.flat_high);
fprintf(fid, 'target_edge_mode: %s\n', string(cfg.mraf.target_edge_mode));
fprintf(fid, 'pivot50 inner_power = %.4f, outer_power = %.4f, remap(0.5)=%.6f\n', cfg.mraf.pivot50_inner_power, cfg.mraf.pivot50_outer_power, target.remap_anchor_50);
fprintf(fid, 'target geometry scale: x_scale=%.5f, y_scale=%.5f\n', cfg.mraf.target_x_scale, cfg.mraf.target_y_scale);
fprintf(fid, 'projection: use_three_region=%d, free_factor=%.3f, outer_factor=%.3f, legacy mraf_factor=%.3f\n', cfg.mraf.use_three_region_projection, cfg.mraf.free_factor, cfg.mraf.outer_factor, cfg.mraf.mraf_factor);
fprintf(fid, 'y_free_reservoir: enable=%d, half_y=%.3f um, max_I=%.3f, pixels=%d\n', cfg.mraf.y_free_reservoir_enable, cfg.mraf.y_free_reservoir_half_y_um, cfg.mraf.y_free_reservoir_max_intensity, mask_count(target, 'y_free_reservoir'));
fprintf(fid, 'wgs_y_edge_damping: enable=%d, factor=%.3f, I band=[%.3f, %.3f], pixels=%d\n', cfg.mraf.wgs_y_edge_damping_enable, cfg.mraf.wgs_y_edge_damping_factor, cfg.mraf.wgs_y_edge_damping_low, cfg.mraf.wgs_y_edge_damping_high, mask_count(target, 'wgs_y_edge_damping_band'));
fprintf(fid, 'edge_mask: %.3f <= scaled target geometry < %.3f inside guard and not free/noise\n', cfg.mraf.free_threshold, cfg.mraf.core_threshold);
fprintf(fid, 'target/free/outer: target=guard & I>=%.3f; free_noise=guard & I<%.3f; outer_suppress=outside guard\n', cfg.mraf.free_threshold, cfg.mraf.free_threshold);
fprintf(fid, 'outer_tail search x: abs(x) > %.3f um; y: abs(y) > %.3f um on center profiles; legacy side_lobe fields equal outer_tail.\n', cfg.side_lobe_exclusion_fraction * cfg.target_half_x_m * 1e6, cfg.side_lobe_exclusion_fraction * cfg.target_half_y_m * 1e6);
fprintf(fid, 'true_side_lobe search: derivative local maxima after 13.5%% crossing + 2.5 um margin, with 5-point smoothing, min height 0.02, min prominence 0.005.\n');
fprintf(fid, 'shoulder search uses 0.70*half50 to 1.05*half50 on both sides of each center profile.\n\n');
write_mask_summary(fid, target, focal_x_m, focal_y_m, cfg);
end

function write_mask_summary(fid, target, focal_x_m, focal_y_m, cfg)
fields = {'plateau','flat_core','core','signal','edge','transition','transition_90_50','transition_50_135','outside_e2','target','free_edge','free_noise','free','side_lobe_reservoir','outer_suppress','noise','guard','target_13p5','finite_target'};
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(target.masks, name), continue; end
    mask = target.masks.(name);
    [width_x_um, width_y_um] = mask_extent_um(mask, focal_x_m, focal_y_m);
    fprintf(fid, '%-14s pixels = %10d, fraction = %.6f, bbox = %.3f x %.3f um\n', name, nnz(mask), nnz(mask)/numel(mask), width_x_um, width_y_um);
end
fprintf(fid, 'outer_tail exclusion half widths: x = %.3f um, y = %.3f um; legacy side_lobe fields equal outer_tail\n', cfg.side_lobe_exclusion_fraction * cfg.target_half_x_m * 1e6, cfg.side_lobe_exclusion_fraction * cfg.target_half_y_m * 1e6);
fprintf(fid, 'true_side_lobe is a derivative local maximum beyond the 13.5%% crossing, not a monotonic tail sample.\n');
fprintf(fid, 'shoulder band is computed from measured half50: [0.70, 1.05] * half50 on left/right and bottom/top.\n');
end

function n = mask_count(target, name)
if isfield(target.masks, name)
    n = nnz(target.masks.(name));
else
    n = 0;
end
end

function [width_x_um, width_y_um] = mask_extent_um(mask, focal_x_m, focal_y_m)
cols = any(mask, 1); rows = any(mask, 2);
if any(cols), width_x_um = (max(focal_x_m(cols)) - min(focal_x_m(cols))) * 1e6; else, width_x_um = NaN; end
if any(rows), width_y_um = (max(focal_y_m(rows)) - min(focal_y_m(rows))) * 1e6; else, width_y_um = NaN; end
end

function write_metric_fields(fid, metrics)
names = metric_names();
for k = 1:numel(names)
    fprintf(fid, '%s = %.10g\n', names{k}, metrics.(names{k}));
end
end

function names = metric_names()
names = {'measured_size90_x_um','measured_size90_y_um','measured_size50_x_um','measured_size50_y_um','measured_size135_x_um','measured_size135_y_um', ...
    'size50_x_um','size50_y_um','size_50_x_um','size_50_y_um','size13p5_x_um','size13p5_y_um', ...
    'transition_13p5_90_x_um','transition_13p5_90_y_um','transition_13p5_90_x_left_um','transition_13p5_90_x_right_um', ...
    'transition_13p5_90_y_left_um','transition_13p5_90_y_right_um','transition_width_13_90_x_um','transition_width_13_90_y_um', ...
    'transition_width_135_90_x_um','transition_width_135_90_y_um', ...
    'transition_width_10_90_x_um','transition_width_10_90_y_um','transition_width_x_um','transition_width_y_um', ...
    'flat_core_mean','flat_core_rms','flat_core_pv','flat_core_uniformity','signal_rms','signal_pv','side_lobe_energy_ratio','outer_energy_ratio','noise_energy_ratio','null_energy_ratio','e2_efficiency', ...
    'shoulder_peak','shoulder_energy_ratio','lobe_peak','lobe_energy_ratio', ...
    'core_rms','rms_90','rms_uniformity_90_region','e2_diffraction_efficiency','shoulder_peak_90_50','shoulder_peak_50_135','peak_to_valley', ...
    'shoulder_peak_x','shoulder_peak_y','shoulder_peak_left','shoulder_peak_right','shoulder_peak_bottom','shoulder_peak_top', ...
    'overshoot_x','overshoot_y','undershoot_x','undershoot_y','outer_tail_peak_x_rel_to_core','outer_tail_peak_y_rel_to_core', ...
    'side_lobe_peak_x_rel_to_core','side_lobe_peak_y_rel_to_core','true_side_lobe_peak_x_rel_to_core','true_side_lobe_peak_y_rel_to_core', ...
    'true_side_lobe_pos_x_um','true_side_lobe_pos_y_um','has_true_side_lobe_x','has_true_side_lobe_y', ...
    'efficiency_inside_signal','efficiency_inside_guard','total_power_conservation_check','has_crossing_warning'};
end

function text = value_to_string(value)
if isstring(value) || ischar(value)
    text = char(value);
elseif isnumeric(value) || islogical(value)
    if isscalar(value), text = sprintf('%.12g', value); else, text = mat2str(value); end
else
    text = sprintf('<%s>', class(value));
end
end
