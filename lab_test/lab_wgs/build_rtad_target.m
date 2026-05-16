function target = build_rtad_target(cfg)
% Build RTAD flat-top rectangular target with all masks.
% References: lab_test_f200mm/src/rtad_target.py:make_rtad_rect_target

N = cfg.N;
dx = cfg.focal_dx_um;
dy = cfg.focal_dy_um;

W50_um = cfg.W50_um;
H50_um = cfg.H50_um;
delta_x_um = cfg.delta_x_um;
delta_y_um = cfg.delta_y_um;
guard_x_um = cfg.guard_x_um;
guard_y_um = cfg.guard_y_um;
release_level = cfg.release_level;

% Coordinate axes
x_axis = ((0:N-1) - N/2) * dx;
y_axis = ((0:N-1) - N/2) * dy;
[X, Y] = meshgrid(x_axis, y_axis);
absX = abs(X);
absY = abs(Y);

% Geometry
a50 = W50_um / 2;
b50 = H50_um / 2;
a0 = a50 - delta_x_um;
a1 = a50 + delta_x_um;
b0 = b50 - delta_y_um;
b1 = b50 + delta_y_um;
a2 = a1 + guard_x_um;
b2 = b1 + guard_y_um;

% Raised cosine edges
Ix = raised_cosine_edge_1d(absX, a0, a1);
Iy = raised_cosine_edge_1d(absY, b0, b1);
I_full = Ix .* Iy;
A_full = sqrt(I_full);

% Masks
mask_flat = (absX <= a0) & (absY <= b0);
mask_template_support = I_full > 0;
mask_signal = I_full >= release_level;
mask_signal = mask_signal | mask_flat;
mask_edge_lock = mask_signal & ~mask_flat;
mask_guard_window = (absX <= a2) & (absY <= b2);
mask_free = mask_guard_window & ~mask_signal;
mask_bg_far = ~(mask_signal | mask_free);

% Compatibility aliases
mask_edge = mask_edge_lock;
mask_support = mask_signal;
mask_bg = mask_bg_far;

% Signal amplitude
A_signal = A_full;
A_signal(~mask_signal) = 0;

% Center profiles
ix0 = N/2 + 1;
iy0 = N/2 + 1;
profile_x_I = I_full(iy0, :);
profile_y_I = I_full(:, ix0)';

% Assemble output struct
target.I_full = I_full;
target.A_full = A_full;
target.A_signal = A_signal;
target.mask_flat = mask_flat;
target.mask_edge_lock = mask_edge_lock;
target.mask_signal = mask_signal;
target.mask_template_support = mask_template_support;
target.mask_guard_window = mask_guard_window;
target.mask_bg_far = mask_bg_far;
target.mask_edge = mask_edge;
target.mask_support = mask_support;
target.mask_bg = mask_bg;
target.mask_free = mask_free;
target.x_um = x_axis;
target.y_um = y_axis;
target.params.W50_um = W50_um;
target.params.H50_um = H50_um;
target.params.a0_um = a0;
target.params.a1_um = a1;
target.params.b0_um = b0;
target.params.b1_um = b1;
target.params.a2_um = a2;
target.params.b2_um = b2;
target.params.release_level = release_level;
end


function C = raised_cosine_edge_1d(u, u0, u1)
% Evaluate a 1D raised-cosine falling edge.
% C=1 for u<=u0, C=0 for u>=u1, cosine transition in between.

C = zeros(size(u), 'single');
C(u <= u0) = 1;
idx = (u > u0) & (u < u1);
t = (u(idx) - u0) / (u1 - u0);
C(idx) = 0.5 * (1 + cos(pi * t));
end
