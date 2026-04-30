# Real-World Simulation Plan

本文件记录第二阶段“真实入射条件 / 真实光路误差模拟”的计划。当前只记录方案，
暂不实施代码。

## Summary

第二阶段目录：

```text
E:\program\Point2P\real_world_simulation
```

这一阶段不重新优化 DOE，而是固定第一阶段最终相位：

```text
E:\program\Point2P\rtad_mraf_gs_python\artifacts\fixed_baseline_bg0p9_initial_compare_20260429-175148\phase_refined.npy
```

然后改变真实光路或入射条件，评估当前 DOE 对误差的敏感性。第一版只做光路和
入射误差，不加入制造误差。

## Core Idea

第一阶段已经完成 DOE 设计流程：

1. 先用 Point2P 生成点对点初始相位 `phase0`。
2. 再以 `phase0` 为初值，用 WGS 优化得到最终 `phase_refined`。

第二阶段要回答的问题是：

```text
如果这片已经设计好的 DOE 放到真实光路里，
入射光、焦平面、孔径和装调条件发生偏差时，
最终匀化光斑会变差多少？
```

所以第二阶段默认是容差评估，而不是重新设计 DOE。

## Planned Project Structure

后续实施时，计划在本文件夹下建立：

```text
real_world_simulation/
  README.md
  config_default.py
  run_real_world_sweep.py
  src/
    field_models.py
    propagation.py
    metrics.py
    plotting.py
  artifacts/
```

各模块职责：

- `config_default.py`：保存 baseline 路径、物理参数和默认 sweep 范围。
- `run_real_world_sweep.py`：运行单个 sweep 或全部第一批 sweep。
- `src/field_models.py`：生成真实入射场，包括光斑尺寸、位置偏移、发散/会聚、
  入射倾角和孔径裁切。
- `src/propagation.py`：复用 nominal Fourier-lens 传播，并新增离焦传播。
- `src/metrics.py`：复用第一阶段诊断指标，并额外记录入射通过率。
- `src/plotting.py`：生成强度图、中心剖面 overlay、指标趋势图和规格书风格 PDF。

## Baseline Inputs

默认输入参数：

```text
wavelength = 532 nm
focal_length = 429 mm
input_gaussian_1e2_diameter = 5 mm
clear_aperture = 15 mm
grid_N = 2048
focal_dx = 2.5 um
dx_doe = 4.457578125e-05 m
```

默认 DOE 相位：

```text
phase_refined.npy
```

默认 target、mask 和 config 从同一个第一阶段 baseline artifact 读取：

```text
E:\program\Point2P\rtad_mraf_gs_python\artifacts\fixed_baseline_bg0p9_initial_compare_20260429-175148
```

## Real Input Field Model

真实入射场在 DOE 平面建模。基本形式为：

```text
E_in(x,y) = A_real(x,y) * exp(i * phi_real(x,y))
E_after_DOE(x,y) = E_in(x,y) * exp(i * phase_refined(x,y))
```

其中 `A_real` 用于描述真实入射强度分布，`phi_real` 用于描述真实波前误差。

入射光先按未裁切高斯归一化为总功率 1，再经过 DOE 通光孔径裁切。这样可以
记录：

```text
aperture_throughput_percent
```

即真实入射光有多少能量真正打进 DOE / 通光孔径。

衍射效率仍使用：

```text
efficiency_e2_percent
```

表示通过 DOE 后的能量中，有多少进入 e^-2 目标区域。后续不再把
`efficiency_flat` 当作总衍射效率。

## First Batch Sweeps

第一批计划采用“两级扫描”：

1. 先用温和范围定位敏感项。
2. 再对明显变差的项目加宽范围，寻找失效边界。

### Defocus / 离焦

观察平面相对理想焦平面的 z 位移：

```text
defocus_mm = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
```

正负号需要在实施时写清楚：正值表示观察面位于理想焦平面之后，负值表示之前。

### Beam Position Offset / 入射光位置偏移

模拟光斑没有完全打到 DOE 中心：

```text
offset_x_mm = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
offset_y_mm = [-2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2]
```

这里的偏移发生在 DOE 平面，会同时改变入射振幅分布和孔径裁切情况。

### Beam Size / 入射光斑尺寸

模拟入射高斯光斑大小变化：

```text
diameter_1e2_mm = [3.5, 4, 4.5, 5, 5.5, 6, 6.5]
```

这里的尺寸是 1/e^2 intensity diameter。

### Divergence / 入射光发散或会聚

用 DOE 面处的波前曲率模拟发散/会聚。参数定义为 1/e^2 半径处的边缘角：

```text
divergence_edge_mrad = [-0.5, -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2, 0.5]
```

正值表示发散，负值表示会聚。

这个参数描述的是波前曲率，不等同于整束光倾斜。

### Pointing / 入射倾角

整束光倾斜会导致焦平面整体平移。第一版用等效焦平面平移表示：

```text
pointing_shift_x_um = [-100, -50, -25, 0, 25, 50, 100]
pointing_shift_y_um = [-100, -50, -25, 0, 25, 50, 100]
```

实施时也可以换算为入射角：

```text
theta_x = pointing_shift_x / focal_length
theta_y = pointing_shift_y / focal_length
```

### Clear Aperture / 有效通光孔径

模拟后续光学元件或 DOE 有效通光区域变小：

```text
clear_aperture_mm = [10, 12, 13, 14, 15]
```

### Ellipticity / 入射椭圆光斑

模拟 x/y 方向光斑尺寸不一致：

```text
(Dx, Dy) mm =
  (5.0, 5.0)
  (4.5, 5.0)
  (5.5, 5.0)
  (5.0, 4.5)
  (5.0, 5.5)
  (4.5, 5.5)
  (5.5, 4.5)
```

## Additional Simulations For Later

以下模拟建议后续再加入，第一版先不做：

- 波长漂移：模拟真实激光中心波长偏离 532 nm。
- 焦距误差：模拟实际聚焦模块焦距和设计值不一致。
- 透镜像差：加入 astigmatism、coma、spherical aberration 等 Zernike 相位。
- DOE 制造误差：相位深度缩放、相位量化、随机相位噪声、局部缺陷。
- DOE 旋转或安装角误差：相位图相对入射光旋转。
- DOE 横向安装偏移：DOE 相位图相对光束中心偏移。
- 光束截断形状误差：圆形孔径、方形孔径、轻微遮挡等。
- 输入光强非高斯：例如 super-Gaussian、椭圆 Gaussian、带热点或暗斑的入射场。

## Planned Outputs

每个 sweep 输出到：

```text
E:\program\Point2P\real_world_simulation\artifacts\<timestamp>_<sweep_name>
```

每个 sweep 保存：

```text
summary.csv
summary.json
center_profiles_overlay.png
metric_trends.png
simulation_report.pdf
representative_intensity_nominal.png
representative_intensity_worst.png
```

`summary.csv` 至少包含：

```text
sweep_name
sweep_parameter
sweep_value
size50_x_um
size50_y_um
size13p5_x_um
size13p5_y_um
transition_13p5_90_x_um
transition_13p5_90_y_um
rms_nonuniformity_percent
efficiency_e2_percent
aperture_throughput_percent
center_offset_x_um
center_offset_y_um
```

## Planned CLI

运行全部第一批 sweep：

```powershell
cd E:\program\Point2P\real_world_simulation
& 'D:\software\anaconda\envs\slmrtad\python.exe' .\run_real_world_sweep.py --sweep all
```

单独运行某一项：

```powershell
--sweep defocus
--sweep beam_offset_x
--sweep beam_offset_y
--sweep beam_size
--sweep divergence
--sweep pointing_x
--sweep pointing_y
--sweep aperture
--sweep ellipticity
```

## Test Plan

第一步先跑 nominal 单点，即所有误差为 0，确认结果接近第一阶段 baseline：

```text
size50_x/y ~= 330.18 / 123.62 um
rms_nonuniformity_percent ~= 1.86%
efficiency_e2_percent ~= 92.52%
```

然后用小网格 smoke test 验证代码路径：

- 不依赖完整 2048 数据也能生成强度图、CSV 和 PDF。
- linear tilt 的焦平面中心偏移方向正确。
- 大位置偏移或小孔径时 `aperture_throughput_percent` 下降。

完整 baseline 测试：

- 每个第一批 sweep 至少生成一个 `summary.csv` 和一张 overlay 图。
- PDF 中指标和 CSV 数值一致。
- `efficiency_e2_percent` 被标注为总衍射效率。
- 不使用 `efficiency_flat` 代替总衍射效率。

## Current Assumptions

- 第二阶段默认是固定 DOE 容差评估，不在误差条件下重新 WGS 优化。
- 第一版只做光路和入射误差，制造误差留到下一轮。
- “入射光发散角”按 DOE 面处波前曲率建模，不等同于整束光倾斜。
- 整束光倾斜单独用 `pointing_shift_x/y_um` 表示。
- 离焦按观察平面相对理想焦平面的 z 位移建模。
- 后续如果某一类误差非常敏感，再单独加密扫描范围。
