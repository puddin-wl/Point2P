function [img, vid] = capture_focal_image(cfg, vid)
% Capture focal-plane image from camera (gentl adapter, Mono8).
%
%   img = capture_focal_image(cfg)        % init camera, capture, close
%   [img, vid] = capture_focal_image(cfg, vid)  % reuse existing camera session
%
% Pattern adapted from cam_in_loop.m.

persistent vid_persistent

if nargin < 2 || isempty(vid)
    if isempty(vid_persistent) || ~isvalid(vid_persistent)
        vid_persistent = init_camera();
    end
    vid = vid_persistent;
end

img = double(getsnapshot(vid));
end


function vid = init_camera()
% Initialize gentl camera with auto-exposure to find a good starting point.
fprintf('[camera] Initializing gentl Mono8 camera...\n');
vid = videoinput("gentl", 1, "Mono8");
src = getselectedsource(vid);
triggerconfig(vid, 'manual');

% Get exposure limits
try
    exp_info = propinfo(src, 'ExposureTime');
    exp_max = min(exp_info.ConstraintValue(2), 1000000);
    exp_min = max(exp_info.ConstraintValue(1), 10);
catch
    exp_max = 1000000;
    exp_min = 10;
end

% Start with moderate exposure
src.ExposureTime = 2000;
start(vid);
pause(0.5);

fprintf('[camera] Ready. Exposure range: [%.0f, %.0f] us\n', exp_min, exp_max);
end
