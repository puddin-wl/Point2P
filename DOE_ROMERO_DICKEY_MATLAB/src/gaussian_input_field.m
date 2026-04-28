function [field, amplitude, intensity, aperture_mask] = gaussian_input_field(X_m, Y_m, cfg)
% gaussian_input_field Build the DOE-plane incident Gaussian field and aperture.
%
% INITIAL PHASE GENERATION:
% Purpose:
%   Generate the input amplitude that accompanies phase0 and the aperture mask
%   that limits where phase0 is physically active.
% Inputs:
%   X_m, Y_m:
%     DOE-plane coordinates from make_grid.m.
%   cfg.input_1e2_diameter_m:
%     Illuminated Gaussian beam diameter at 1/e^2 intensity. Default: 5 mm.
%   cfg.aperture_diameter_m:
%     Clear aperture / computational pupil diameter. Default: 15 mm. This is
%     not the same as the illuminated beam diameter.
% Outputs:
%   field:
%     Real input field amplitude after applying the clear aperture.
%   amplitude, intensity:
%     Gaussian amplitude and intensity on the DOE plane after aperture masking.
%   aperture_mask:
%     Circular clear aperture mask used by build_separable_phase_2d.m.
% Physical meaning:
%   The initial phase is designed for a Gaussian input beam. This function
%   supplies that beam envelope; the phase itself is added later by
%   field .* exp(1i*phase0).
% Used by:
%   run_one.m, run_mraf_one.m, export_phase_for_python.m, and
%   run_initial_phase_only.m.
% Notes:
%   The 5 mm default is the 1/e^2 intensity diameter, so w below is the 1/e^2
%   intensity radius and amplitude = exp(-r^2/w^2).

r2_m2 = X_m.^2 + Y_m.^2;
w_m = cfg.input_1e2_radius_m;
intensity = exp(-2 * r2_m2 / w_m^2);
amplitude = exp(-r2_m2 / w_m^2);
aperture_mask = r2_m2 <= cfg.aperture_radius_m^2;
field = amplitude .* aperture_mask;
intensity = intensity .* aperture_mask;
amplitude = amplitude .* aperture_mask;
end
