function [X_m, Y_m, x_m, y_m, grid] = make_grid(cfg)
% make_grid 生成以零为中心的 DOE 面二维物理坐标网格。

N = cfg.N;
dx_m = cfg.dx_doe_m;
index = (-N/2):(N/2 - 1);
x_m = index * dx_m;
y_m = index * dx_m;
[X_m, Y_m] = meshgrid(x_m, y_m);

grid = struct();
grid.N = N;
grid.dx_doe_m = dx_m;
grid.extent_m = N * dx_m;
grid.x_min_m = min(x_m);
grid.x_max_m = max(x_m);
grid.y_min_m = min(y_m);
grid.y_max_m = max(y_m);
end

