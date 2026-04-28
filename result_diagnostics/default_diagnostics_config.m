function cfg = default_diagnostics_config()
% default_diagnostics_config Minimal config for focal-plane result diagnostics.
%
% This config is intentionally independent from MRAF target/projection code.
% It only describes the physical target size and default profile levels used
% when processing an already-computed focal-plane intensity result.

cfg = struct();
cfg.target_size_x_m = 330e-6;
cfg.target_size_y_m = 120e-6;
cfg.target_half_x_m = cfg.target_size_x_m / 2;
cfg.target_half_y_m = cfg.target_size_y_m / 2;

cfg.core_fraction = 0.45;
cfg.guard_fraction = 1.50;
cfg.side_lobe_exclusion_fraction = 1.15;
cfg.profile_levels = [0.90, 0.50, 0.135];
cfg.transition_low = 0.135;
cfg.transition_high = 0.90;
cfg.shoulder_inner_fraction = 0.70;
cfg.shoulder_outer_fraction = 1.05;
cfg.true_side_lobe_min_height = 0.02;
cfg.true_side_lobe_min_prominence = 0.005;
cfg.true_side_lobe_smooth_window = 5;
cfg.figure_dpi = 150;
cfg.center_profile_xlim_factor = 2.2;
end
