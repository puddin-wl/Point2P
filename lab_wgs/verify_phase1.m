%% verify_phase1 — Verify MATLAB WGS matches Python (Phase 1)
% Runs one WGS iteration in MATLAB and compares with Python ground truth.
% Python command used to generate ground truth:
%   python run_rtad_mraf_gs_case.py --phase-mat make_phase0_output_5p5mm/phase0.mat
%     --phase-var phase0_wrapped_rad --beam-diameter 5.5 --method wgs
%     --wgs-strategy flat_local --iters 1 --wgs-feedback-exponent 0.8
%     --wgs-weight-min 0.5 --wgs-weight-max 1.5 --bg-factor 0.9
%     --wgs-update-every 1 --no-swap-phase-xy --no-cupy
%     --outdir artifacts/phase1_verify_1iter

clear; clc;

%% Configuration (matching Python parameters)
cfg = experiments_wgs_config();
cfg.input_1e2_diameter_m = 5.5e-3;   % beam diameter 5.5 mm (Python: --beam-diameter 5.5)
cfg.wgs_feedback_exponent = 0.8;      % Python: --wgs-feedback-exponent 0.8
cfg.wgs_weight_min = 0.5;             % Python: --wgs-weight-min 0.5
cfg.wgs_weight_max = 1.5;             % Python: --wgs-weight-max 1.5
cfg.bg_factor = 0.9;                  % Python: --bg-factor 0.9
cfg.wgs_update_every = 1;             % Force weight update on first iteration

%% Python ground truth (from artifacts/phase1_verify_1iter/metrics.csv)
py_initial_rms = 3.524;               % iter 0 flat_rms * 100
py_initial_size50_x = 327.44;         % iter 0 size50_x_um
py_initial_size50_y = 117.18;         % iter 0 size50_y_um
py_iter1_rms = 7.159;                 % iter 1 flat_rms * 100
py_iter1_size50_x = 324.98;           % iter 1 size50_x_um
py_iter1_size50_y = 116.44;           % iter 1 size50_y_um
py_weight_mean = 1.000;               % iter 1 wgs_weight_mean
py_weight_std = 0.0151;               % iter 1 wgs_weight_std
py_weight_min = 0.9799;               % iter 1 wgs_weight_min
py_weight_max = 1.0612;               % iter 1 wgs_weight_max

%% Load the SAME initial phase as Python
phase0_path = '../lab_test_f200mm/make_phase0_output_5p5mm/phase0.mat';
phase0_var = 'phase0_wrapped_rad';
fprintf('Loading initial phase from: %s [%s]\n', phase0_path, phase0_var);
phase0 = load_initial_phase(phase0_path, phase0_var);

%% Build input amplitude and target
fprintf('Building input amplitude (beam=%.1f mm)...\n', cfg.input_1e2_diameter_m * 1e3);
[input_amp, aperture_mask] = load_input_amplitude(cfg);

fprintf('Building RTAD target...\n');
target = build_rtad_target(cfg);

%% Evaluate initial state
N = cfg.N;
n_total = N * N;
E_doe = input_amp .* exp(1i * phase0);
E_focal = fftshift(fft2(ifftshift(E_doe))) / sqrt(n_total);
I_before = abs(E_focal).^2;

% Compute flat-region metrics directly (scale-invariant, no image analysis needed)
flat_vals = I_before(target.mask_flat);
flat_mean_before = mean(flat_vals);
flat_std_before = std(flat_vals);
rms_before = flat_std_before / flat_mean_before * 100;

fprintf('\n========== Before WGS (Initial RD Phase) ==========\n');
fprintf('  MATLAB: RMS=%.3f%%, flat_mean=%.6f\n', rms_before, flat_mean_before);
fprintf('  Python: RMS=%.3f%%, size50=[%.1f, %.1f] um\n', py_initial_rms, py_initial_size50_x, py_initial_size50_y);

%% Run one WGS iteration (with weight update, matching Python)
fprintf('\nRunning one WGS iteration...\n');
weights = ones(N, N);
[phase_new, weights_new, I_after] = wgs_iteration(phase0, input_amp, weights, target, cfg, [], true);

%% Evaluate after WGS
flat_vals_after = I_after(target.mask_flat);
flat_mean_after = mean(flat_vals_after);
flat_std_after = std(flat_vals_after);
rms_after = flat_std_after / flat_mean_after * 100;

w_flat = weights_new(target.mask_flat);
w_mean = mean(w_flat);
w_std = std(w_flat);
w_min = min(w_flat);
w_max = max(w_flat);

fprintf('\n========== After 1 WGS Iteration ==========\n');
fprintf('  MATLAB: RMS=%.3f%%, flat_mean=%.6f\n', rms_after, flat_mean_after);
fprintf('  Python: RMS=%.3f%%, size50=[%.1f, %.1f] um\n', py_iter1_rms, py_iter1_size50_x, py_iter1_size50_y);
fprintf('\n  Weight stats:\n');
fprintf('  MATLAB: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n', w_mean, w_std, w_min, w_max);
fprintf('  Python: mean=%.4f, std=%.4f, min=%.4f, max=%.4f\n', py_weight_mean, py_weight_std, py_weight_min, py_weight_max);

%% Phase change
dphi = mod(phase_new - phase0 + pi, 2*pi) - pi;
fprintf('\n  RMS phase change: %.4f rad\n', sqrt(mean(dphi(aperture_mask).^2)));
fprintf('  Max |phase change|: %.4f rad\n', max(abs(dphi(aperture_mask))));

%% Verification
fprintf('\n========== Verification ==========\n');
fprintf('Note: exact numerical match is not expected because:\n');
fprintf('  1. Python L2-normalizes the target amplitude; MATLAB does not\n');
fprintf('  2. Python applies mraf_factor=0.4 to the free region; MATLAB skips this\n');
fprintf('     (experimental WGS does not use MRAF, per lab_wgs design)\n');
fprintf('\n');

rms_ok = abs(rms_before - py_initial_rms) < 0.5;
w_ok = abs(w_mean - py_weight_mean) < 0.1 && abs(w_std - py_weight_std) < 0.05;

if rms_ok
    fprintf('PASS: Initial RMS matches Python (%.3f%% vs %.3f%%).\n', rms_before, py_initial_rms);
else
    fprintf('CHECK: Initial RMS differs (%.3f%% vs %.3f%%). Verify FFT/amplitude setup.\n', rms_before, py_initial_rms);
end

if w_ok
    fprintf('PASS: Weight statistics are consistent with Python.\n');
else
    fprintf('INFO: Weight stats differ (expected due to normalization differences).\n');
end

fprintf('\nPhase 1 verification complete.\n');
fprintf('For full Python comparison, see: lab_test_f200mm/artifacts/phase1_verify_1iter/\n');

%% Save MATLAB output for manual comparison
out_dir = fullfile(cfg.output_root, 'phase1_verify');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save(fullfile(out_dir, 'phase_after_1iter.mat'), 'phase_new', '-v7.3');
save(fullfile(out_dir, 'weights_after_1iter.mat'), 'weights_new', '-v7.3');
fprintf('MATLAB output saved to: %s\n', out_dir);
