function result = analyze_captured_image(img, pixel_um)
% Analyze a captured flat-top spot image using the gradient-edge method.
%
%   result = analyze_captured_image(img, pixel_um)
%
%   img: 2D image array (double)
%   pixel_um: camera pixel size in um (default 3.45)
%
%   result fields: background, flat_center_px, flat_mean, flat_std, rms_pct,
%                  size50_px, size50_um, size90_px, size13p5_px, aspect_ratio,
%                  saturation_pct, hotspot_ratio_x/y, gradient_edges_x/y
%
% References: fig_analysis/analyze_captured.py

if nargin < 2
    pixel_um = 3.45;
end

img = double(img);
[ny, nx] = size(img);

%% Background estimation (four corners)
corners = [reshape(img(1:30, 1:30), [], 1);
           reshape(img(1:30, nx-29:nx), [], 1);
           reshape(img(ny-29:ny, 1:30), [], 1);
           reshape(img(ny-29:ny, nx-29:nx), [], 1)];
bg_val = median(corners);
bg_std = std(corners);

%% Find bright region
threshold = bg_val + 10 * bg_std;
bright = img > threshold;

% Morphological closing (binary_closing)
bright = bwmorph(bright, 'close');
bright = imfill(bright, 'holes');

% Label connected components
[labeled, nf] = bwlabel(bright);

if nf == 0
    warning('No bright region found in image.');
    result = struct('background', bg_val, 'error', 'No bright region found');
    return;
end

% Select largest component
comp_sizes = zeros(1, nf);
for i = 1:nf
    comp_sizes(i) = sum(labeled(:) == i);
end
[~, main_id] = max(comp_sizes);
mask = labeled == main_id;

%% Center of mass
[y_idx, x_idx] = find(mask);
vals = img(mask);
total_mass = sum(vals);
if total_mass > 0
    cx = sum(x_idx .* vals) / total_mass;
    cy = sum(y_idx .* vals) / total_mass;
else
    cx = mean(x_idx);
    cy = mean(y_idx);
end

ymin = min(y_idx); ymax = max(y_idx);
xmin = min(x_idx); xmax = max(x_idx);
bbox_w = xmax - xmin + 1;
bbox_h = ymax - ymin + 1;

%% Flat region statistics
flat_vals = img(mask);
flat_mean = mean(flat_vals);
flat_std = std(flat_vals);
flat_min = min(flat_vals);
flat_max = max(flat_vals);
n_sat = sum(flat_vals >= 255);

rms_pct = flat_std / flat_mean * 100;
pv_pct = (flat_max - flat_min) / flat_mean * 100;
saturation_pct = 100 * n_sat / numel(flat_vals);

%% Gradient-edge width analysis
cxi = round(cx);
cyi = round(cy);

prof_x = img(cyi, :);
prof_y = img(:, cxi)';

% Smooth profiles (7-pixel uniform filter)
prof_x_s = movmean(prof_x, 7);
prof_y_s = movmean(prof_y, 7);

% Background on profiles
x_bg = median([prof_x(1:min(50,nx)), prof_x(max(1,nx-49):nx)]);
y_bg = median([prof_y(1:min(50,ny)), prof_y(max(1,ny-49):ny)]);

% Gradient magnitude
gx = abs(gradient(prof_x_s));
gy = abs(gradient(prof_y_s));

% Find gradient edges bracketing the flat-top
[lx, rx] = find_edge_pair(gx, cxi);
[ly, ry] = find_edge_pair(gy, cyi);

% Flat level = median between gradient edges (exclude edge transition zones)
margin = 8;
flat_x_idx = max(1, lx+margin) : min(nx, rx-margin);
flat_y_idx = max(1, ly+margin) : min(ny, ry-margin);

if ~isempty(flat_x_idx)
    flat_level_x = median(prof_x_s(flat_x_idx));
    hotspot_x = max(prof_x_s(flat_x_idx)) / flat_level_x;
else
    flat_level_x = median(flat_vals);
    hotspot_x = 1.0;
end

if ~isempty(flat_y_idx)
    flat_level_y = median(prof_y_s(flat_y_idx));
    hotspot_y = max(prof_y_s(flat_y_idx)) / flat_level_y;
else
    flat_level_y = median(flat_vals);
    hotspot_y = 1.0;
end

% Find width crossings
x_cross = find_crossings(prof_x_s, x_bg, flat_level_x, cxi);
y_cross = find_crossings(prof_y_s, y_bg, flat_level_y, cyi);

w50_x = x_cross.W50;
w50_y = y_cross.W50;
w90_x = x_cross.W90;
w90_y = y_cross.W90;
w13_x = x_cross.W135;
w13_y = y_cross.W135;

aspect = w50_x / w50_y;

%% Assemble result
result = struct();
result.background = bg_val;
result.background_std = bg_std;
result.flat_center_px = [cx, cy];
result.flat_bbox = [xmin, ymin, bbox_w, bbox_h];
result.flat_mean = flat_mean;
result.flat_std = flat_std;
result.rms_percent = rms_pct;
result.pv_percent = pv_pct;
result.saturated_percent = saturation_pct;
result.gradient_edges_x_px = [lx, rx];
result.gradient_edges_y_px = [ly, ry];
result.flat_level_x = flat_level_x;
result.flat_level_y = flat_level_y;
result.hotspot_ratio_x = hotspot_x;
result.hotspot_ratio_y = hotspot_y;
result.size50_px = [w50_x, w50_y];
result.size90_px = [w90_x, w90_y];
result.size13p5_px = [w13_x, w13_y];
result.aspect_ratio = aspect;
result.pixel_um = pixel_um;
result.size50_um = [w50_x * pixel_um, w50_y * pixel_um];
result.size13p5_um = [w13_x * pixel_um, w13_y * pixel_um];
result.overexposed = saturation_pct > 1;
result.n_saturated = n_sat;
end


function [left_edge, right_edge] = find_edge_pair(grad, center)
% Find the two dominant gradient peaks bracketing the flat-top.
N = length(grad);
search_radius = 300;

left_start = max(1, center - search_radius);
left_grad = grad(left_start:min(N, center));
right_grad = grad(max(1, center):min(N, center + search_radius));

thr = 0.15 * max(grad);

left_peaks = find_peaks(left_grad, thr, left_start - 1);
right_peaks = find_peaks(right_grad, thr, max(1, center) - 1);

if isempty(left_peaks)
    left_edge = max(1, center - 50);
else
    left_peaks = sortrows(left_peaks, 2, 'descend');
    left_edge = min(left_peaks(1:min(3,end), 1));
end

if isempty(right_peaks)
    right_edge = min(N, center + 50);
else
    right_peaks = sortrows(right_peaks, 2, 'descend');
    right_edge = max(right_peaks(1:min(3,end), 1));
end
end


function peaks = find_peaks(arr, threshold, offset)
% Find local maxima above threshold.
peaks = zeros(0, 2);
n = length(arr);
for i = 3:n-2
    if arr(i) > threshold && ...
       arr(i) >= arr(i-1) && arr(i) >= arr(i-2) && ...
       arr(i) > arr(i+1) && arr(i) > arr(i+2)
        peaks(end+1, :) = [i + offset, arr(i)]; %#ok<AGROW>
    end
end
end


function crossings = find_crossings(prof, bg, flat_level, center)
% Find 90%, 50%, 13.5% threshold crossings relative to flat level.
crossings = struct();

for level_info = [struct('label', 'W90',  'frac', 0.9);
                   struct('label', 'W50',  'frac', 0.5);
                   struct('label', 'W135', 'frac', 0.135)]'
    frac = level_info.frac;
    threshold = bg + frac * (flat_level - bg);
    n = length(prof);

    % Left crossing
    left_pos = center;
    for i = center:-1:2
        if prof(i) < threshold
            if abs(prof(i) - prof(i+1)) > 0.1
                left_pos = i + (threshold - prof(i)) / (prof(i+1) - prof(i));
            else
                left_pos = i;
            end
            break;
        end
    end

    % Right crossing
    right_pos = center;
    for i = center:n-1
        if prof(i+1) < threshold
            if abs(prof(i+1) - prof(i)) > 0.1
                right_pos = i + 1 - (prof(i+1) - threshold) / (prof(i+1) - prof(i));
            else
                right_pos = i + 1;
            end
            break;
        end
    end

    crossings.(level_info.label) = right_pos - left_pos;
end
end
