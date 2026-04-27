function [target_intensity, info] = rectangular_target(focal_x_m, focal_y_m, cfg)
% rectangular_target 在焦平面坐标上生成矩形平顶目标强度。
%
% hard 模式为硬矩形；logistic 模式使用可调边缘宽度生成软边。

[XF_m, YF_m] = meshgrid(focal_x_m, focal_y_m);
half_x_m = cfg.target_half_x_m;
half_y_m = cfg.target_half_y_m;

switch lower(char(cfg.target_edge))
    case 'hard'
        target_intensity = double(abs(XF_m) <= half_x_m & abs(YF_m) <= half_y_m);
    case 'logistic'
        sx_m = max(cfg.transition_width_x_m, eps);
        sy_m = max(cfg.transition_width_y_m, eps);
        edge_x = 1 ./ (1 + exp((abs(XF_m) - half_x_m) ./ sx_m));
        edge_y = 1 ./ (1 + exp((abs(YF_m) - half_y_m) ./ sy_m));
        target_intensity = edge_x .* edge_y;
        target_intensity = target_intensity ./ max(target_intensity(:));
    otherwise
        error('Unknown target_edge: %s', cfg.target_edge);
end

info = struct();
info.edge = cfg.target_edge;
info.half_x_m = half_x_m;
info.half_y_m = half_y_m;
info.size_x_m = cfg.target_size_x_m;
info.size_y_m = cfg.target_size_y_m;
end

