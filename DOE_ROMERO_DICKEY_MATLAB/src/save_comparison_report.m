function save_comparison_report(cfg, ideal_metrics, full_metrics)
% save_comparison_report 保存 ideal FFT baseline 与 full propagation 的关键指标对比。

path = fullfile(cfg.save_root, 'comparison_ideal_vs_full.txt');
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Ideal FFT baseline vs full propagation comparison\n');
fprintf(fid, '=================================================\n\n');
fprintf(fid, 'Ideal baseline: DOE plane is assumed to be the Fourier transform pupil plane.\n');
fprintf(fid, 'Full propagation: DOE -> %.6g m free space -> thin lens -> focal plane.\n\n', cfg.doe_to_lens_m);
fprintf(fid, 'Full model lens aperture enabled = %d\n', cfg.apply_lens_aperture_in_full_model);
fprintf(fid, 'Full model lens aperture diameter = %.9g m\n\n', cfg.lens_aperture_diameter_m);
write_metric_row(fid, 'FWHM_x_um', ideal_metrics.fwhm_x_m * 1e6, full_metrics.fwhm_x_m * 1e6);
write_metric_row(fid, 'FWHM_y_um', ideal_metrics.fwhm_y_m * 1e6, full_metrics.fwhm_y_m * 1e6);
write_metric_row(fid, 'transition_13_90_x_um', ideal_metrics.transition_width_13_90_x_m * 1e6, full_metrics.transition_width_13_90_x_m * 1e6);
write_metric_row(fid, 'transition_13_90_y_um', ideal_metrics.transition_width_13_90_y_m * 1e6, full_metrics.transition_width_13_90_y_m * 1e6);
write_metric_row(fid, 'core_rms', ideal_metrics.core_rms, full_metrics.core_rms);
write_metric_row(fid, 'peak_to_valley', ideal_metrics.peak_to_valley, full_metrics.peak_to_valley);
write_metric_row(fid, 'efficiency_target', ideal_metrics.efficiency_target, full_metrics.efficiency_target);
write_metric_row(fid, 'efficiency_guard', ideal_metrics.efficiency_guard, full_metrics.efficiency_guard);
write_metric_row(fid, 'side_lobe_x', ideal_metrics.side_lobe_peak_x_rel_to_core, full_metrics.side_lobe_peak_x_rel_to_core);
write_metric_row(fid, 'side_lobe_y', ideal_metrics.side_lobe_peak_y_rel_to_core, full_metrics.side_lobe_peak_y_rel_to_core);
write_metric_row(fid, 'overshoot_x', ideal_metrics.overshoot_near_edge_x, full_metrics.overshoot_near_edge_x);
write_metric_row(fid, 'overshoot_y', ideal_metrics.overshoot_near_edge_y, full_metrics.overshoot_near_edge_y);
fprintf(fid, '\nInterpretation:\n');
fprintf(fid, '- Use ideal_* outputs first to judge whether the Romero-Dickey separable phase is behaving as intended.\n');
fprintf(fid, '- Use full_* outputs to evaluate degradation caused by the 200 mm DOE-to-lens propagation and real lens placement.\n');
fprintf(fid, '- With infinite lens aperture and ideal scalar propagation, free-space propagation mainly adds phase and the focal intensity can be nearly identical to the ideal FFT baseline.\n');
fprintf(fid, '- The default full model therefore applies a configurable lens aperture to expose practical clipping/deviation effects.\n');
end

function write_metric_row(fid, name, ideal_value, full_value)
fprintf(fid, '%-28s ideal = %.9g, full = %.9g, full-ideal = %.9g\n', name, ideal_value, full_value, full_value - ideal_value);
end
