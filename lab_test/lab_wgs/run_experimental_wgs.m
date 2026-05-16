function [phase_final, history] = run_experimental_wgs(cfg, initial_phase, input_amp, target, sim_mode)
% Run the experimental WGS feedback loop.
%
%   [phase_final, history] = run_experimental_wgs(cfg, initial_phase, input_amp, target)
%   runs in simulation mode (FFT both directions).
%
%   [phase_final, history] = run_experimental_wgs(cfg, initial_phase, input_amp, target, false)
%   runs in hardware mode: SLM -> camera -> analyze -> WGS -> SLM ...
%
% References: lab_wgs/readme.md algorithm section, cam_in_loop.m hardware interface

if nargin < 5, sim_mode = true; end

N = cfg.N;
max_iters = cfg.max_iters;
update_every = cfg.wgs_update_every;
mask_flat = target.mask_flat;

% Initialize
phase = initial_phase;
weights = ones(N, N);

% Compute initial focal plane
E_doe = input_amp .* exp(1i * phase);
E_focal = fftshift(fft2(ifftshift(E_doe))) / sqrt(N*N);
I_focal = abs(E_focal).^2;

% Direct RMS from mask_flat
flat_vals = I_focal(mask_flat);
rms0 = std(flat_vals) / mean(flat_vals);

history.iter = 0;
history.rms = rms0;
history.phase_snapshots = {};
history.exposure_us = [];

fprintf('=== Experimental WGS Loop (%s mode) ===\n', ...
    conditional(sim_mode, 'SIMULATION', 'HARDWARE'));
fprintf('Iter %2d: RMS=%.3f%%\n', 0, rms0 * 100);

%% Hardware init
if ~sim_mode
    % Send initial phase to SLM
    display_slm_phase(phase, cfg);
    pause(1.0);

    % Init camera
    vid = capture_focal_image(cfg, []);

    % Auto-exposure: first mask out zero-order region from initial snap
    img_init = double(getsnapshot(vid));
    [~, max_idx] = max(img_init(:));
    [cy, cx] = ind2sub(size(img_init), max_idx);
    [cam_H, cam_W] = size(img_init);
    [CamX, CamY] = meshgrid(1:cam_W, 1:cam_H);
    mask_no_zero = ((CamX - cx).^2 + (CamY - cy).^2) > 150^2;

    auto_exposure(vid, 200, mask_no_zero);
    history.exposure_us(1) = getselectedsource(vid).ExposureTime;
end

%% Main loop
for it = 1:max_iters
    do_update = (it > 0) && (mod(it, update_every) == 0);

    if sim_mode
        [phase, weights, I_focal] = wgs_iteration(phase, input_amp, weights, target, cfg, [], do_update);
    else
        % 1. Capture camera image
        img = capture_focal_image(cfg, vid);

        % 2. Analyze: locate flat-top region
        result = analyze_captured_image(img, cfg.camera_pixel_um);
        if isfield(result, 'error')
            warning('analyze_captured_image failed: %s. Skipping weight update.', result.error);
            continue;
        end

        % 3. Map camera ROI to simulation grid
        bbox = [result.flat_bbox(1), result.flat_bbox(2), ...
                result.flat_bbox(3), result.flat_bbox(4)];
        [I_meas, ~, ~] = map_camera_to_focal(img, result.flat_center_px, bbox, cfg);

        % 4. WGS iteration with measured intensity
        [phase, weights, I_focal] = wgs_iteration(phase, input_amp, weights, target, cfg, I_meas, do_update);

        % 5. Display new phase on SLM
        display_slm_phase(phase, cfg);
        pause(0.3);

        % 6. Check saturation
        if result.saturated_percent > 1
            fprintf('  [WARN] %.1f%% saturated! Reduce exposure.\n', result.saturated_percent);
        end
    end

    % RMS from flat region
    flat_vals = I_focal(mask_flat);
    rms_now = std(flat_vals) / mean(flat_vals);

    history.iter(end+1) = it;
    history.rms(end+1) = rms_now;

    fprintf('Iter %2d: RMS=%.3f%%', it, rms_now * 100);
    if do_update
        wf = weights(mask_flat);
        fprintf('  [W: mean=%.3f std=%.3f]', mean(wf), std(wf));
    end
    fprintf('\n');

    % Convergence check
    if rms_now < cfg.rms_target
        fprintf('Target RMS %.2f%% reached at iteration %d.\n', cfg.rms_target * 100, it);
        break;
    end
    if it > cfg.convergence_window
        recent = history.rms(end-cfg.convergence_window+1:end);
        if max(recent) - min(recent) < 1e-4 && rms_now < 0.05
            fprintf('Converged.\n');
            break;
        end
    end

    if mod(it, 10) == 0
        history.phase_snapshots{end+1} = phase;
    end
end

%% Cleanup
if ~sim_mode && exist('vid', 'var') && isvalid(vid)
    stop(vid);
    delete(vid);
    clear vid;
    imaqreset;
    fprintf('Camera released.\n');
end

% Save
phase_final = phase;
if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
save(fullfile(cfg.output_root, 'phase_final.mat'), 'phase_final', '-v7.3');
save(fullfile(cfg.output_root, 'wgs_history.mat'), 'history', '-v7.3');
fprintf('Saved to %s/\n', cfg.output_root);
end


function s = conditional(cond, t, f)
if cond, s = t; else, s = f; end
end
