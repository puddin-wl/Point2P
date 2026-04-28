function write_diagnostics_report(output_path, diagnostics)
% write_diagnostics_report Write key result diagnostics to a text file.

fid = fopen(output_path, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Focal result diagnostics\n');
fprintf(fid, '========================\n\n');
fprintf(fid, 'normalization = %s\n', string(diagnostics.normalization));
fprintf(fid, 'core_mean_raw = %.10g\n', diagnostics.core_mean_raw);
fprintf(fid, 'core_rms = %.10g\n', diagnostics.core_rms);
fprintf(fid, 'core_pv = %.10g\n', diagnostics.core_pv);
fprintf(fid, 'core_uniformity = %.10g\n\n', diagnostics.core_uniformity);

fprintf(fid, 'size90_x_um = %.10g\n', diagnostics.size90_x_um);
fprintf(fid, 'size90_y_um = %.10g\n', diagnostics.size90_y_um);
fprintf(fid, 'size50_x_um = %.10g\n', diagnostics.size50_x_um);
fprintf(fid, 'size50_y_um = %.10g\n', diagnostics.size50_y_um);
fprintf(fid, 'size13p5_x_um = %.10g\n', diagnostics.size13p5_x_um);
fprintf(fid, 'size13p5_y_um = %.10g\n\n', diagnostics.size13p5_y_um);

fprintf(fid, 'transition_13p5_90_x_um = %.10g\n', diagnostics.transition_13p5_90_x_um);
fprintf(fid, 'transition_13p5_90_y_um = %.10g\n', diagnostics.transition_13p5_90_y_um);
fprintf(fid, 'transition_13p5_90_x_left_um = %.10g\n', diagnostics.transition_13p5_90_x_left_um);
fprintf(fid, 'transition_13p5_90_x_right_um = %.10g\n', diagnostics.transition_13p5_90_x_right_um);
fprintf(fid, 'transition_13p5_90_y_bottom_um = %.10g\n', diagnostics.transition_13p5_90_y_bottom_um);
fprintf(fid, 'transition_13p5_90_y_top_um = %.10g\n\n', diagnostics.transition_13p5_90_y_top_um);

fprintf(fid, 'shoulder_peak_x = %.10g\n', diagnostics.shoulder_peak_x);
fprintf(fid, 'shoulder_peak_y = %.10g\n', diagnostics.shoulder_peak_y);
fprintf(fid, 'side_lobe_peak_x_rel_to_core = %.10g\n', diagnostics.side_lobe_peak_x_rel_to_core);
fprintf(fid, 'side_lobe_peak_y_rel_to_core = %.10g\n', diagnostics.side_lobe_peak_y_rel_to_core);
fprintf(fid, 'true_side_lobe_peak_x_rel_to_core = %.10g\n', diagnostics.true_side_lobe_peak_x_rel_to_core);
fprintf(fid, 'true_side_lobe_pos_x_um = %.10g\n', diagnostics.true_side_lobe_pos_x_um);
fprintf(fid, 'has_true_side_lobe_x = %.0f\n', diagnostics.has_true_side_lobe_x);
fprintf(fid, 'true_side_lobe_peak_y_rel_to_core = %.10g\n', diagnostics.true_side_lobe_peak_y_rel_to_core);
fprintf(fid, 'true_side_lobe_pos_y_um = %.10g\n', diagnostics.true_side_lobe_pos_y_um);
fprintf(fid, 'has_true_side_lobe_y = %.0f\n\n', diagnostics.has_true_side_lobe_y);

fprintf(fid, 'energy_inside_signal = %.10g\n', diagnostics.energy_inside_signal);
fprintf(fid, 'energy_inside_guard = %.10g\n', diagnostics.energy_inside_guard);
fprintf(fid, 'noise_energy_ratio = %.10g\n', diagnostics.noise_energy_ratio);
fprintf(fid, 'null_energy_ratio = %.10g\n', diagnostics.null_energy_ratio);
fprintf(fid, 'total_power = %.10g\n', diagnostics.total_power);
end
