% run_all_beam_diameters Run initial phase generation for 5, 6, and 7 mm beams.
%
% Run this script from this folder or from MATLAB with this folder on path.

function run_all_beam_diameters()
    module_root = fileparts(mfilename('fullpath'));
    addpath(module_root);

    beam_diameters_mm = [5, 6, 7];

    for i = 1:length(beam_diameters_mm)
        d_mm = beam_diameters_mm(i);

        cfg = default_initial_phase_config(module_root);
        cfg.input_1e2_diameter_m = d_mm * 1e-3;
        cfg.input_1e2_radius_m = cfg.input_1e2_diameter_m / 2;
        cfg.input_1e_radius_m = cfg.input_1e2_radius_m / sqrt(2);
        cfg.beta_x = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_x_m / (cfg.lambda_m * cfg.f_m);
        cfg.beta_y = 2 * pi * cfg.input_1e_radius_m * cfg.Ro_y_m / (cfg.lambda_m * cfg.f_m);

        stamp = string(datetime('now', 'Format', 'yyyyMMdd-HHmmss'));
        out_dir = fullfile(cfg.output_root, char(stamp));

        fprintf('\n========================================\n');
        fprintf('Generating initial phase for D = %d mm beam\n', d_mm);
        fprintf('beta_x = %.6f, beta_y = %.6f\n', cfg.beta_x, cfg.beta_y);
        fprintf('f = %.1f mm, dx_doe = %.4f um\n', cfg.f_m * 1e3, cfg.dx_doe_m * 1e6);
        fprintf('========================================\n');

        phase_data = generate_initial_phase(cfg, ...
            'do_forward', true, ...
            'output_dir', out_dir, ...
            'figure_dpi', cfg.figure_dpi);

        fprintf('Output folder: %s\n', out_dir);
        fprintf('phase0 size: %d x %d\n\n', ...
            size(phase_data.phase0_unwrapped_rad, 1), ...
            size(phase_data.phase0_unwrapped_rad, 2));
    end

    fprintf('All three beam diameters completed.\n');
end
