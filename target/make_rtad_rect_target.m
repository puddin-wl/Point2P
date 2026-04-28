function T = make_rtad_rect_target(varargin)
%MAKE_RTAD_RECT_TARGET Generate a rectangular RTAD target for flat-top DOE refinement.
%
% This function generates only the target amplitude/intensity and masks.
% It does not run GS/MRAF/WGS iterations.
%
% Default target:
%   - size50 = 330 um x 120 um, interpreted as INTENSITY 50% width/height.
%   - descending edge = raised-cosine transition.
%   - output amplitude target = sqrt(output intensity target).
%
% Coordinate convention:
%   - Matrix rows correspond to y/eta.
%   - Matrix columns correspond to x/xi.
%   - T.I_target(row, col) = I(eta, xi).
%
% Example:
%   T = make_rtad_rect_target( ...
%       'N', [2048 2048], ...
%       'dx_um', 0.5, ...
%       'dy_um', 0.5, ...
%       'W50_um', 330, ...
%       'H50_um', 120, ...
%       'delta_x_um', 15, ...
%       'delta_y_um', 8, ...
%       'mode', 'separable', ...
%       'make_plots', true, ...
%       'save_dir', 'artifacts/rtad_target_test');
%
% Important:
%   MRAF normally constrains amplitude, while many beam-size metrics use
%   intensity. This function fixes the 50% boundary in INTENSITY, then
%   converts intensity to amplitude by A = sqrt(I).

p = inputParser;
p.FunctionName = 'make_rtad_rect_target';

addParameter(p, 'N', [2048 2048]);             % [Ny Nx] or scalar
addParameter(p, 'dx_um', 1.0);                 % output-plane x pixel size, um
addParameter(p, 'dy_um', 1.0);                 % output-plane y pixel size, um
addParameter(p, 'x_um', []);                   % optional x axis vector, um
addParameter(p, 'y_um', []);                   % optional y axis vector, um

addParameter(p, 'W50_um', 330.0);              % intensity 50% full width, um
addParameter(p, 'H50_um', 120.0);              % intensity 50% full height, um
addParameter(p, 'delta_x_um', 15.0);           % half transition width around x 50% boundary, um
addParameter(p, 'delta_y_um', 8.0);            % half transition width around y 50% boundary, um
addParameter(p, 'center_x_um', 0.0);           % target center x, um
addParameter(p, 'center_y_um', 0.0);           % target center y, um

% mode:
%   'separable' : I = Cx(x) * Cy(y). Smoothest corners; first recommended test.
%   'box'       : I = C(max-normalized-distance). Strict rectangular 50% contour, including corners.
addParameter(p, 'mode', 'separable');

% Free/guard band is not part of target intensity. It is only output as a mask
% for later MRAF usage.
addParameter(p, 'guard_x_um', 20.0);
addParameter(p, 'guard_y_um', 12.0);

addParameter(p, 'make_plots', false);
addParameter(p, 'save_dir', '');
addParameter(p, 'prefix', 'rtad_rect');

parse(p, varargin{:});
S = p.Results;

% -------------------------
% Build coordinate grid
% -------------------------
if isempty(S.x_um) || isempty(S.y_um)
    N = S.N;
    if isscalar(N)
        Ny = N;
        Nx = N;
    else
        Ny = N(1);
        Nx = N(2);
    end

    x_um = ((1:Nx) - (Nx + 1) / 2) * S.dx_um;
    y_um = ((1:Ny) - (Ny + 1) / 2) * S.dy_um;
else
    x_um = S.x_um(:).';
    y_um = S.y_um(:);
    Ny = numel(y_um);
    Nx = numel(x_um);
end

[X, Y] = meshgrid(x_um - S.center_x_um, y_um - S.center_y_um);
absX = abs(X);
absY = abs(Y);

% -------------------------
% Target geometry
% -------------------------
a50 = S.W50_um / 2;
b50 = S.H50_um / 2;

dx_edge = S.delta_x_um;
dy_edge = S.delta_y_um;

if dx_edge <= 0 || dy_edge <= 0
    error('delta_x_um and delta_y_um must be positive.');
end

if dx_edge >= a50 || dy_edge >= b50
    error('delta_x_um/delta_y_um are too large: core half-size would become non-positive.');
end

a0 = a50 - dx_edge;
a1 = a50 + dx_edge;
b0 = b50 - dy_edge;
b1 = b50 + dy_edge;

mode = lower(string(S.mode));

% -------------------------
% Generate intensity target
% -------------------------
switch mode
    case "separable"
        Ix = raised_cosine_edge(absX, a0, a1);
        Iy = raised_cosine_edge(absY, b0, b1);
        I_target = Ix .* Iy;

        mask_flat = absX <= a0 & absY <= b0;
        mask_support = I_target > 0;
        mask_edge = mask_support & ~mask_flat;

    case "box"
        % q = 0 at the flat-core boundary, q = 1 at the size50 boundary,
        % q = 2 at the outer zero boundary.
        qx = (absX - a0) / dx_edge;
        qy = (absY - b0) / dy_edge;
        q = max(qx, qy);

        I_target = zeros(Ny, Nx);
        I_target(q <= 0) = 1;

        trans = q > 0 & q < 2;
        I_target(trans) = 0.5 * (1 + cos(pi * q(trans) / 2));
        I_target(q >= 2) = 0;

        mask_flat = q <= 0;
        mask_support = q < 2;
        mask_edge = mask_support & ~mask_flat;

    otherwise
        error('Unknown mode "%s". Use "separable" or "box".', S.mode);
end

% Numeric safety.
I_target = max(0, min(1, I_target));
A_target = sqrt(I_target);

mask_bg = ~mask_support;

% Optional guard/free ring for later MRAF use.
a2 = a1 + S.guard_x_um;
b2 = b1 + S.guard_y_um;
mask_guard_support = absX <= a2 & absY <= b2;
mask_free = mask_guard_support & ~mask_support;

% -------------------------
% Profiles and summary
% -------------------------
[~, ix0] = min(abs(x_um - S.center_x_um));
[~, iy0] = min(abs(y_um - S.center_y_um));

profile_x_I = I_target(iy0, :);
profile_y_I = I_target(:, ix0);
profile_x_A = A_target(iy0, :);
profile_y_A = A_target(:, ix0);

measured.W50_x_um = estimate_width_at_level(x_um - S.center_x_um, profile_x_I, 0.5);
measured.H50_y_um = estimate_width_at_level(y_um - S.center_y_um, profile_y_I, 0.5);
measured.W90_x_um = estimate_width_at_level(x_um - S.center_x_um, profile_x_I, 0.9);
measured.H90_y_um = estimate_width_at_level(y_um - S.center_y_um, profile_y_I, 0.9);
measured.W10_x_um = estimate_width_at_level(x_um - S.center_x_um, profile_x_I, 0.1);
measured.H10_y_um = estimate_width_at_level(y_um - S.center_y_um, profile_y_I, 0.1);

params = struct();
params.W50_um = S.W50_um;
params.H50_um = S.H50_um;
params.a50_um = a50;
params.b50_um = b50;
params.delta_x_um = dx_edge;
params.delta_y_um = dy_edge;
params.a0_um = a0;
params.a1_um = a1;
params.b0_um = b0;
params.b1_um = b1;
params.guard_x_um = S.guard_x_um;
params.guard_y_um = S.guard_y_um;
params.a2_um = a2;
params.b2_um = b2;
params.mode = char(mode);
params.Ny = Ny;
params.Nx = Nx;
params.dx_um = median(diff(x_um));
params.dy_um = median(diff(y_um));
params.center_x_um = S.center_x_um;
params.center_y_um = S.center_y_um;

T = struct();
T.I_target = I_target;
T.A_target = A_target;
T.mask_flat = mask_flat;
T.mask_edge = mask_edge;
T.mask_support = mask_support;
T.mask_bg = mask_bg;
T.mask_free = mask_free;
T.mask_guard_support = mask_guard_support;
T.x_um = x_um;
T.y_um = y_um;
T.params = params;
T.measured = measured;
T.profiles.x_um = x_um - S.center_x_um;
T.profiles.y_um = y_um - S.center_y_um;
T.profiles.x_I = profile_x_I;
T.profiles.y_I = profile_y_I;
T.profiles.x_A = profile_x_A;
T.profiles.y_A = profile_y_A;

% -------------------------
% Save outputs / plots
% -------------------------
if ~isempty(S.save_dir)
    if ~exist(S.save_dir, 'dir')
        mkdir(S.save_dir);
    end
    save(fullfile(S.save_dir, [S.prefix '_target.mat']), 'T', '-v7.3');
end

if S.make_plots
    make_target_plots(T, S.save_dir, S.prefix);
end

print_summary(T);

end

% =========================================================================
% Local helpers
% =========================================================================
function C = raised_cosine_edge(u, u0, u1)
%RAISED_COSINE_EDGE 1 inside u<=u0, 0 outside u>=u1, cosine transition.
C = zeros(size(u));
C(u <= u0) = 1;
idx = u > u0 & u < u1;
t = (u(idx) - u0) / (u1 - u0);
C(idx) = 0.5 * (1 + cos(pi * t));
C(u >= u1) = 0;
end

function width = estimate_width_at_level(axis_um, profile, level)
%ESTIMATE_WIDTH_AT_LEVEL Estimate full width by linear interpolation.
% Assumes a centered, roughly single-peaked profile.
axis_um = axis_um(:).';
profile = profile(:).';

[~, i0] = min(abs(axis_um));

% Left crossing.
left_axis = axis_um(1:i0);
left_prof = profile(1:i0);
idxL = find(left_prof < level, 1, 'last');
if isempty(idxL) || idxL == numel(left_prof)
    xL = NaN;
else
    x1 = left_axis(idxL);     y1 = left_prof(idxL);
    x2 = left_axis(idxL + 1); y2 = left_prof(idxL + 1);
    xL = interp_crossing(x1, y1, x2, y2, level);
end

% Right crossing.
right_axis = axis_um(i0:end);
right_prof = profile(i0:end);
idxR = find(right_prof < level, 1, 'first');
if isempty(idxR) || idxR == 1
    xR = NaN;
else
    x1 = right_axis(idxR - 1); y1 = right_prof(idxR - 1);
    x2 = right_axis(idxR);     y2 = right_prof(idxR);
    xR = interp_crossing(x1, y1, x2, y2, level);
end

if isnan(xL) || isnan(xR)
    width = NaN;
else
    width = xR - xL;
end
end

function x = interp_crossing(x1, y1, x2, y2, level)
if abs(y2 - y1) < eps
    x = 0.5 * (x1 + x2);
else
    x = x1 + (level - y1) * (x2 - x1) / (y2 - y1);
end
end

function make_target_plots(T, save_dir, prefix)
x = T.x_um;
y = T.y_um;
P = T.params;

fig1 = figure('Color', 'w', 'Name', 'RTAD target intensity');
imagesc(x, y, T.I_target);
axis image;
set(gca, 'YDir', 'normal');
colorbar;
xlabel('x / um');
ylabel('y / um');
title(sprintf('RTAD target intensity, mode=%s', P.mode));
hold on;
contour(x, y, T.I_target, [0.5 0.5], 'LineWidth', 1.2);
plot_rect(P.a0_um, P.b0_um, '--');
plot_rect(P.a50_um, P.b50_um, '-');
plot_rect(P.a1_um, P.b1_um, ':');
legend({'I=0.5 contour', 'flat core', 'size50', 'outer edge'}, 'Location', 'bestoutside');

fig2 = figure('Color', 'w', 'Name', 'RTAD center profiles');
plot(T.profiles.x_um, T.profiles.x_I, 'LineWidth', 1.4);
hold on;
plot(T.profiles.y_um, T.profiles.y_I, 'LineWidth', 1.4);
yline(0.5, ':');
yline(0.9, ':');
yline(0.1, ':');
xline(-P.a0_um, '--'); xline(P.a0_um, '--');
xline(-P.a50_um, '-'); xline(P.a50_um, '-');
xline(-P.a1_um, ':'); xline(P.a1_um, ':');
xline(-P.b0_um, '--'); xline(P.b0_um, '--');
xline(-P.b50_um, '-'); xline(P.b50_um, '-');
xline(-P.b1_um, ':'); xline(P.b1_um, ':');
grid on;
xlabel('coordinate / um');
ylabel('normalized intensity');
title('RTAD center profiles');
legend({'x center profile', 'y center profile'}, 'Location', 'best');

fig3 = figure('Color', 'w', 'Name', 'RTAD masks');
mask_show = zeros(size(T.I_target));
mask_show(T.mask_flat) = 1;
mask_show(T.mask_edge) = 2;
mask_show(T.mask_free) = 3;
imagesc(x, y, mask_show);
axis image;
set(gca, 'YDir', 'normal');
colorbar;
xlabel('x / um');
ylabel('y / um');
title('Masks: 1=flat, 2=edge, 3=free guard band');

if ~isempty(save_dir)
    exportgraphics(fig1, fullfile(save_dir, [prefix '_intensity.png']), 'Resolution', 180);
    exportgraphics(fig2, fullfile(save_dir, [prefix '_profiles.png']), 'Resolution', 180);
    exportgraphics(fig3, fullfile(save_dir, [prefix '_masks.png']), 'Resolution', 180);
end
end

function plot_rect(a, b, style)
xx = [-a, a, a, -a, -a];
yy = [-b, -b, b, b, -b];
plot(xx, yy, style, 'LineWidth', 1.2);
end

function print_summary(T)
P = T.params;
M = T.measured;
fprintf('\nRTAD rectangular target generated.\n');
fprintf('  mode              : %s\n', P.mode);
fprintf('  grid              : Ny=%d, Nx=%d, dx=%.6g um, dy=%.6g um\n', P.Ny, P.Nx, P.dx_um, P.dy_um);
fprintf('  target size50      : W50=%.3f um, H50=%.3f um\n', P.W50_um, P.H50_um);
fprintf('  measured size50    : W50=%.3f um, H50=%.3f um\n', M.W50_x_um, M.H50_y_um);
fprintf('  flat core halfsize : a0=%.3f um, b0=%.3f um\n', P.a0_um, P.b0_um);
fprintf('  size50 halfsize    : a50=%.3f um, b50=%.3f um\n', P.a50_um, P.b50_um);
fprintf('  outer halfsize     : a1=%.3f um, b1=%.3f um\n', P.a1_um, P.b1_um);
fprintf('  guard halfsize     : a2=%.3f um, b2=%.3f um\n', P.a2_um, P.b2_um);
fprintf('  mask pixels        : flat=%d, edge=%d, free=%d, bg=%d\n\n', ...
    nnz(T.mask_flat), nnz(T.mask_edge), nnz(T.mask_free), nnz(T.mask_bg));
end
