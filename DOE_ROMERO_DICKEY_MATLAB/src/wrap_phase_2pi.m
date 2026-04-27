function phase_wrapped_rad = wrap_phase_2pi(phase_rad)
% wrap_phase_2pi 将弧度相位包裹到 [0, 2*pi)。

phase_wrapped_rad = mod(phase_rad, 2*pi);
end

