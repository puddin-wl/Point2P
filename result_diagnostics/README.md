# MATLAB 焦面结果诊断

本目录只分析已经计算好的焦面强度，不生成相位、不运行 MRAF/WGS，也不修改目标。

## 文件

- `compute_focal_diagnostics.m`：计算中心剖面、阈值宽度、过渡宽度、核心 RMS/PV、肩峰、旁瓣和能量比例；
- `plot_center_profiles_diagnostics.m`：绘制中心剖面；
- `write_diagnostics_report.m`：写出文本报告；
- `default_diagnostics_config.m`：默认阈值和期望尺寸；
- `run_diagnostics_example.m`：示例入口。

## 基本用法

```matlab
cfg = default_diagnostics_config();
diagnostics = compute_focal_diagnostics( ...
    focal_x_m, focal_y_m, intensity, 'cfg', cfg);
plot_center_profiles_diagnostics(diagnostics, 'output.png');
write_diagnostics_report(diagnostics, 'report.txt');
```

## 指标约定

- `size50`、`size90` 和 `size13p5` 按归一化强度阈值求中心剖面交点；
- 核心 RMS 只在固定平顶核心区域计算；
- e⁻² 效率和核心区功率不是同一个指标；
- 旁瓣必须在主过渡区之外寻找，不能把边缘本身当作旁瓣。

输出默认写入 `artifacts/`，不进入 Git。

