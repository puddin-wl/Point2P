function phase = load_initial_phase(phase_path, var_name)
% Load simulation WGS refined phase from .mat or .npy file.
%
%   phase = load_initial_phase(phase_path) loads from a .mat file (auto-detect
%   variable) or .npy file.
%
%   phase = load_initial_phase(phase_path, var_name) uses the specified variable
%   name for .mat files.
%
% References: lab_test_f200mm/artifacts/<ts>/phase_refined.npy or .mat

if nargin < 2
    var_name = '';
end

[~, ~, ext] = fileparts(phase_path);

if strcmpi(ext, '.npy')
    phase = read_npy(phase_path);
else
    data = load(phase_path);
    if isempty(var_name)
        % Auto-detect: find the 2048x2048 variable
        fns = fieldnames(data);
        found = false;
        for i = 1:numel(fns)
            v = data.(fns{i});
            if isnumeric(v) && ismatrix(v) && all(size(v) == 2048)
                phase = double(v);
                found = true;
                break;
            end
        end
        if ~found
            % Fallback: use first numeric 2D array
            for i = 1:numel(fns)
                v = data.(fns{i});
                if isnumeric(v) && ismatrix(v) && size(v,1) > 100
                    phase = double(v);
                    found = true;
                    break;
                end
            end
        end
        if ~found
            error('No 2048x2048 array found in %s. Specify var_name.', phase_path);
        end
    else
        phase = double(data.(var_name));
    end
end

% Ensure [0, 2π) wrapping and replace NaN (aperture exterior) with 0
phase(isnan(phase)) = 0;
phase = mod(phase, 2*pi);
end


function data = read_npy(filepath)
% Minimal .npy reader for float32/float64 arrays.
% Only handles C-order, little-endian files written by numpy.save.

fid = fopen(filepath, 'rb');
if fid < 0
    error('Cannot open file: %s', filepath);
end
cleanup = onCleanup(@() fclose(fid));

% Magic string
magic = fread(fid, 6, 'uint8=>char')';
if ~strcmp(magic, sprintf('\x93NUMPY'))
    error('Not a valid .npy file: %s', filepath);
end

% Header
major = fread(fid, 1, 'uint8');
minor = fread(fid, 1, 'uint8');
header_len = fread(fid, 1, 'uint16');
header_str = fread(fid, header_len, 'uint8=>char')';

% Parse dtype and shape from header dict (Python literal)
% Header looks like: {'descr': '<f4', 'fortran_order': False, 'shape': (2048, 2048)}
dtype_match = regexp(header_str, '''descr'':\s*''(\S+)''', 'tokens');
shape_match = regexp(header_str, '''shape'':\s*\(([^)]*)\)', 'tokens');

if isempty(dtype_match) || isempty(shape_match)
    error('Cannot parse .npy header: %s', header_str);
end

dtype_str = dtype_match{1}{1};
shape_str = shape_match{1}{1};

shape_parts = strsplit(strtrim(shape_str), ',');
dims = zeros(1, numel(shape_parts));
for i = 1:numel(shape_parts)
    part = strtrim(shape_parts{i});
    if ~isempty(part)
        dims(i) = str2double(part);
    end
end
dims(dims == 0) = [];

% Determine MATLAB type
switch dtype_str
    case {'<f4', '<f'}
        matlab_type = 'float32';
        byte_size = 4;
    case {'<f8', '<f8'}
        matlab_type = 'float64';
        byte_size = 8;
    case {'<i4'}
        matlab_type = 'int32';
        byte_size = 4;
    case {'<i8'}
        matlab_type = 'int64';
        byte_size = 8;
    otherwise
        error('Unsupported dtype: %s', dtype_str);
end

data = fread(fid, prod(dims), [matlab_type '=>' matlab_type]);
data = reshape(data, dims);
end
