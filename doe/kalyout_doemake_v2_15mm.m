clc;
clear;
close all;

% V2 WGS phase -> 15 mm x 15 mm square, four-mask binary DOE layout.
% Requires the MATLAB GDSII toolbox functions:
%   gds_structure, gds_element, gds_library, write_gds_library
%
% Important: GDS dbunit is a coordinate precision, not an etch depth.
% The physical etch depths are documented per GDS layer below.

script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, 'v2_wgs_15mm_doe_20260730');
data_file = fullfile(data_dir, 'doe_phase_v2_wgs_15mm_1024.mat');
output_file = fullfile(data_dir, 'doe_v2_wgs_15mm_square_4mask.gds');

required_functions = {'gds_structure', 'gds_element', ...
    'gds_library', 'write_gds_library'};
for k = 1:numel(required_functions)
    assert(~isempty(which(required_functions{k})), ...
        'Missing MATLAB GDSII toolbox function: %s', required_functions{k});
end

data = load(data_file, '-mat');
assert(isfield(data, 'doe_phase'), 'Input MAT must contain doe_phase.');
pupil_start = double(data.doe_phase);

Npupil = 1024;
DoeSizeUm = 15000.0;
PixelPitchUm = DoeSizeUm / Npupil; % 14.6484375 um/pixel
assert(isequal(size(pupil_start), [Npupil, Npupil]), ...
    'doe_phase must be 1024 x 1024.');

% Match the supplied kalyout_doemake.m material and quantization settings.
wavelength_nm = 532.0;
n_substrate = 1.458;
n_ambient = 1.00029;
pupil_height_nm = 77.5;

% Full square phase: intentionally no circular aperture mask.
pupil_rezero = pupil_start - min(pupil_start(:));
pupilrad = mod(pupil_rezero, 2*pi);
pupildata_nm = wavelength_nm / (n_substrate - n_ambient) ...
    * pupilrad / (2*pi);
pupil_step = floor(pupildata_nm / pupil_height_nm);
pupil_step = uint8(min(max(pupil_step, 0), 15));

% Four binary masks. Layer depth increments are:
%   layer 1: 8 * 77.5 = 620.0 nm
%   layer 2: 4 * 77.5 = 310.0 nm
%   layer 3: 2 * 77.5 = 155.0 nm
%   layer 4: 1 * 77.5 =  77.5 nm
mask_bits = [8, 4, 2, 1];
gs = gds_structure('DOE_V2_WGS_15MM');

for layer_idx = 1:numel(mask_bits)
    bit_value = mask_bits(layer_idx);
    layer_mask = bitand(pupil_step, uint8(bit_value)) ~= 0;
    polygons = cell(0, 1);

    % Merge adjacent selected pixels into horizontal rectangles. Matrix rows
    % are y and columns are x, matching the simulation's axis convention.
    for row = 1:Npupil
        transitions = diff([false, layer_mask(row, :), false]);
        starts = find(transitions == 1);
        stops = find(transitions == -1) - 1;
        for run_idx = 1:numel(starts)
            x0 = (starts(run_idx) - 1) * PixelPitchUm;
            x1 = stops(run_idx) * PixelPitchUm;
            y0 = (row - 1) * PixelPitchUm;
            y1 = row * PixelPitchUm;
            polygons{end + 1, 1} = [ ... %#ok<SAGROW>
                x0, y0; x1, y0; x1, y1; x0, y1; x0, y0];
        end
    end
    gs(end + 1) = gds_element('boundary', 'xy', polygons, ...
        'layer', layer_idx);
    fprintf('Layer %d (bit %d, %.1f nm): %d rectangles\n', ...
        layer_idx, bit_value, bit_value * pupil_height_nm, numel(polygons));
end

% Reference outline only; layer 100 is not an etch mask.
outline_width_um = 1.0;
outline = {
    [0, 0; DoeSizeUm, 0; DoeSizeUm, outline_width_um; ...
        0, outline_width_um; 0, 0];
    [0, DoeSizeUm-outline_width_um; DoeSizeUm, DoeSizeUm-outline_width_um; ...
        DoeSizeUm, DoeSizeUm; 0, DoeSizeUm; 0, DoeSizeUm-outline_width_um];
    [0, outline_width_um; outline_width_um, outline_width_um; ...
        outline_width_um, DoeSizeUm-outline_width_um; ...
        0, DoeSizeUm-outline_width_um; 0, outline_width_um];
    [DoeSizeUm-outline_width_um, outline_width_um; ...
        DoeSizeUm, outline_width_um; DoeSizeUm, DoeSizeUm-outline_width_um; ...
        DoeSizeUm-outline_width_um, DoeSizeUm-outline_width_um; ...
        DoeSizeUm-outline_width_um, outline_width_um]
};
gs(end + 1) = gds_element('boundary', 'xy', outline, 'layer', 100);

% Coordinates are in micrometres; database precision is 1 nm.
glib = gds_library('DOE_V2_WGS_15MM_LIB', ...
    'uunit', 1e-6, 'dbunit', 1e-9, gs);
write_gds_library(glib, output_file);

fprintf('Wrote %s\n', output_file);
fprintf('DOE size: %.6f mm x %.6f mm\n', DoeSizeUm/1000, DoeSizeUm/1000);
fprintf('Pixel pitch: %.7f um\n', PixelPitchUm);
