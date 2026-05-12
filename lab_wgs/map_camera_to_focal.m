function [I_mapped, x_sim_idx, y_sim_idx] = map_camera_to_focal(I_cam, cam_center_px, cam_bbox, cfg)
% Map camera image to simulation focal-plane grid.
%
%   [I_mapped, x_idx, y_idx] = map_camera_to_focal(I_cam, cam_center_px, cam_bbox, cfg)
%
%   I_cam: camera image (2D array)
%   cam_center_px: [cx, cy] center of flat-top in camera pixels
%   cam_bbox: [xmin, ymin, width, height] bounding box in camera pixels
%   cfg: configuration struct from experiments_wgs_config
%
%   I_mapped: intensity mapped onto the simulation 2048×2048 grid (NaN outside ROI)
%   x_sim_idx, y_sim_idx: simulation grid index ranges corresponding to camera ROI
%
% The mapping uses the center-of-mass and bounding box from analyze_captured_image,
% then maps camera pixels onto the focal-plane physical coordinates (um), which
% directly correspond to the simulation grid.

N = cfg.N;
focal_dx_um = cfg.focal_dx_um;
pixel_um = cfg.camera_pixel_um;

cx_cam = cam_center_px(1);
cy_cam = cam_center_px(2);
xmin = cam_bbox(1);
ymin = cam_bbox(2);
bbox_w = cam_bbox(3);
bbox_h = cam_bbox(4);

% Bounding box center in camera pixels
bbox_cx = xmin + bbox_w / 2;
bbox_cy = ymin + bbox_h / 2;

% Offset of flat-top center relative to bbox center (in camera pixels)
offset_x_px = cx_cam - bbox_cx;
offset_y_px = cy_cam - bbox_cy;

% Offset in um on the focal plane
offset_x_um = offset_x_px * pixel_um;
offset_y_um = offset_y_px * pixel_um;

% Simulation grid center index
sim_center = N / 2 + 1;

% Mapping: camera pixel position (relative to center) → focal um → simulation index
% Camera pixel i → physical position = (i - cx_cam) * pixel_um
% Simulation index j → physical position = (j - sim_center) * focal_dx_um
%
% So: (i - cx_cam) * pixel_um = (j - sim_center) * focal_dx_um
% → j = sim_center + (i - cx_cam) * pixel_um / focal_dx_um

scale = pixel_um / focal_dx_um;

% Determine simulation grid bounds for the camera ROI
cam_cols = xmin:xmin+bbox_w-1;
cam_rows = ymin:ymin+bbox_h-1;

sim_x_start = round(sim_center + (cam_cols(1) - cx_cam) * scale);
sim_x_end = round(sim_center + (cam_cols(end) - cx_cam) * scale);
sim_y_start = round(sim_center + (cam_rows(1) - cy_cam) * scale);
sim_y_end = round(sim_center + (cam_rows(end) - cy_cam) * scale);

% Clip to valid simulation grid
sim_x_start = max(1, min(N, sim_x_start));
sim_x_end = max(1, min(N, sim_x_end));
sim_y_start = max(1, min(N, sim_y_start));
sim_y_end = max(1, min(N, sim_y_end));

x_sim_idx = sim_x_start:sim_x_end;
y_sim_idx = sim_y_start:sim_y_end;

% Initialize mapped intensity with NaN
I_mapped = nan(N, N);

% Map: for each camera pixel in ROI, find nearest simulation pixel
% Use sub2ind-friendly 2D interpolation
[XX_cam, YY_cam] = meshgrid(cam_cols, cam_rows);
XX_sim = sim_center + (XX_cam - cx_cam) * scale;
YY_sim = sim_center + (YY_cam - cy_cam) * scale;

% Round to nearest simulation pixel
XX_sim = round(XX_sim);
YY_sim = round(YY_sim);

% Bin camera pixels into simulation grid (average multiple hits)
I_accum = zeros(N, N);
counts = zeros(N, N);

roi = I_cam(ymin:ymin+bbox_h-1, xmin:xmin+bbox_w-1);

for iy = 1:bbox_h
    for ix = 1:bbox_w
        sx = XX_sim(iy, ix);
        sy = YY_sim(iy, ix);
        if sx >= 1 && sx <= N && sy >= 1 && sy <= N
            I_accum(sy, sx) = I_accum(sy, sx) + double(roi(iy, ix));
            counts(sy, sx) = counts(sy, sx) + 1;
        end
    end
end

mask = counts > 0;
I_mapped(mask) = I_accum(mask) ./ counts(mask);
end
