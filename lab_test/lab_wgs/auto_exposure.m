function [exposure_us, peak_val] = auto_exposure(vid, target_peak, mask)
% Auto-adjust camera exposure to reach target peak intensity.
%
%   [exposure_us, peak_val] = auto_exposure(vid, target_peak, mask)
%
%   vid: videoinput object (must be started)
%   target_peak: desired peak pixel value (default 200)
%   mask: optional logical mask (same size as image) to exclude regions
%
% Pattern: iterative binary search, adapted from cam_in_loop.m.

if nargin < 2 || isempty(target_peak), target_peak = 200; end
if nargin < 3, mask = []; end

src = getselectedsource(vid);
try
    exp_info = propinfo(src, 'ExposureTime');
    exp_max = min(exp_info.ConstraintValue(2), 1000000);
    exp_min = max(exp_info.ConstraintValue(1), 10);
catch
    exp_max = 1000000;
    exp_min = 10;
end

max_attempts = 20;
for attempt = 1:max_attempts
    img_test = double(getsnapshot(vid));
    if ~isempty(mask)
        img_test = img_test .* mask;
    end
    current_peak = max(img_test(:));

    if current_peak >= 250
        src.ExposureTime = max(round(src.ExposureTime * 0.6), exp_min);
    elseif current_peak < target_peak - 20
        gain = min((target_peak + 20) / max(current_peak, 1), 2.5);
        src.ExposureTime = min(round(src.ExposureTime * gain), exp_max);
    else
        break;
    end
    pause(0.05);
end

exposure_us = src.ExposureTime;
peak_val = max(double(getsnapshot(vid)), [], 'all');

fprintf('[exposure] %.0f us, peak=%.0f (target=%.0f)\n', exposure_us, peak_val, target_peak);
end
