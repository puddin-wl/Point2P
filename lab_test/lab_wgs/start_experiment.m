%% start_experiment — Run experimental WGS starting from simulation-refined phase
% This is the intended workflow:
%   1. Load the WGS-refined phase from Stage 2 (simulation, RMS ~0.18%)
%   2. Run experimental WGS loop to correct real-world errors
%   3. Output the final lab-optimized phase

clear; clc;

%% Configuration
cfg = experiments_wgs_config();

%% Load simulation-refined phase (Stage 2 output)
% Default: f200mm, beam=7mm, RMS=0.18%, 200 WGS iterations
phase_path = cfg.phase0_path;
if isempty(phase_path) || ~exist(phase_path, 'file')
    error('Phase file not found: %s\nSet cfg.phase0_path in experiments_wgs_config.m', phase_path);
end

fprintf('Loading simulation-refined phase...\n');
fprintf('  Path: %s\n', phase_path);
phase0 = load_initial_phase(phase_path);
fprintf('  Shape: %d x %d, range: [%.4f, %.4f] rad\n', size(phase0, 1), size(phase0, 2), min(phase0(:)), max(phase0(:)));

%% Build input amplitude and target
fprintf('Building input amplitude (beam=%.1f mm)...\n', cfg.input_1e2_diameter_m * 1e3);
input_amp = load_input_amplitude(cfg);

fprintf('Building RTAD target...\n');
target = build_rtad_target(cfg);

%% Preview initial focal plane
N = cfg.N;
E_doe = input_amp .* exp(1i * phase0);
E_focal = fftshift(fft2(ifftshift(E_doe))) / sqrt(N*N);
I_initial = abs(E_focal).^2;

flat_vals = I_initial(target.mask_flat);
rms_initial = std(flat_vals) / mean(flat_vals) * 100;
fprintf('\nInitial state (simulation-refined phase):\n');
fprintf('  RMS in flat region: %.2f%%\n', rms_initial);
fprintf('  Expected: ~0.18%% (from simulation WGS)\n');

%% Run experimental WGS
sim_mode = true;   % true=simulation (FFT), false=hardware (camera+SLM)

fprintf('\n========== Starting experimental WGS ==========\n');
if sim_mode, mode_str = 'SIMULATION'; else, mode_str = 'HARDWARE'; end
fprintf('Mode: %s\n', mode_str);
fprintf('Max iterations: %d, target RMS: %.1f%%\n', cfg.max_iters, cfg.rms_target * 100);
fprintf('Feedback exponent: %.1f, weight range: [%.1f, %.1f]\n', ...
    cfg.wgs_feedback_exponent, cfg.wgs_weight_min, cfg.wgs_weight_max);
fprintf('Update every: %d iterations\n\n', cfg.wgs_update_every);

[phase_final, history] = run_experimental_wgs(cfg, phase0, input_amp, target, sim_mode);

fprintf('\n========== Done ==========\n');
fprintf('Final RMS: %.2f%% (initial: %.2f%%)\n', history.rms(end) * 100, rms_initial);
fprintf('Final phase saved to: %s/phase_final.mat\n', cfg.output_root);
