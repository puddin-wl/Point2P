function phase_wrapped = wrap_mraf_phase(phase_rad)
phase_wrapped = mod(phase_rad, 2*pi);
end