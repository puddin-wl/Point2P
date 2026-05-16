function gray_img = display_slm_phase(phase, cfg)
% Convert computational DOE phase to SLM-loadable 8-bit grayscale BMP,
% write to disk, and send to SLM via SecondDll.
%
%   gray_img = display_slm_phase(phase, cfg)
%
%   phase: 2048x2048 phase array in [0, 2pi)
%   cfg: configuration struct from experiments_wgs_config
%   gray_img: 1080x1920 uint8 image sent to SLM
%
% Steps:
%   1. Crop central region matching SLM physical size (12.288 x 6.912 mm)
%   2. Complex-field cubic interpolation to 1920x1080 (avoids 2pi-wrap ringing)
%   3. Wrap to [0, 2pi), map to [0, 255]
%   4. Write BMP, call SecondDll to display on SLM
%
% References: lab_test_f200mm/save_slm_phase.py, cam_in_loop.m SLM interface

N_comp = cfg.N;
dx_doe_um = cfg.dx_doe_m * 1e6;
slm_w = cfg.slm_width;
slm_h = cfg.slm_height;
slm_pitch = cfg.slm_pitch_um;

% Physical size of SLM active area in mm
slm_phys_x_mm = slm_w * slm_pitch / 1000;
slm_phys_y_mm = slm_h * slm_pitch / 1000;

% Number of computational pixels that fit in SLM physical area
crop_w = floor(slm_phys_x_mm / (dx_doe_um / 1000));
crop_h = floor(slm_phys_y_mm / (dx_doe_um / 1000));

% Crop central region
cx = N_comp / 2;
cy = N_comp / 2;
x0 = cx - floor(crop_w / 2);
x1 = x0 + crop_w - 1;
y0 = cy - floor(crop_h / 2);
y1 = y0 + crop_h - 1;

phase_cropped = phase(y0:y1, x0:x1);

% Complex-field interpolation (continuous across 2pi phase wraps)
[y_crop, x_crop] = size(phase_cropped);
x_out = linspace(1, x_crop, slm_w);
y_out = linspace(1, y_crop, slm_h);
[X_out, Y_out] = meshgrid(x_out, y_out);

complex_field = exp(1i * phase_cropped);
real_part = interp2(double(real(complex_field)), X_out, Y_out, 'cubic');
imag_part = interp2(double(imag(complex_field)), X_out, Y_out, 'cubic');
phase_slm = atan2(imag_part, real_part);

% Wrap and convert to 8-bit
phase_slm = mod(phase_slm, 2 * pi);
gray_img = uint8(phase_slm / (2 * pi) * 255);

fprintf('[SLM] %dx%d -> %dx%d (crop %dx%d, zoom %.3fx %.3fx) | gray [%d,%d]\n', ...
    N_comp, N_comp, slm_w, slm_h, crop_w, crop_h, ...
    slm_w/crop_w, slm_h/crop_h, min(gray_img(:)), max(gray_img(:)));

% Write BMP and send to SLM
slm_dir = fullfile(cfg.output_root, 'slm_phases');
if ~exist(slm_dir, 'dir'), mkdir(slm_dir); end
bmp_path = fullfile(slm_dir, 'current_phase.bmp');
imwrite(gray_img, bmp_path, 'bmp');

% Load SLM DLL if needed
if ~libisloaded('SecondDll')
    loadlibrary('SecondDll.dll', 'SecondDll.h');
end

calllib('SecondDll', 'saShowImageFromFilePath', ...
    bmp_path, 0, slm_w, 0, slm_w, slm_h, 1);
end
