# Romero-Dickey MATLAB DOE 矩形平顶相位设计

本项目用于基于 Romero & Dickey, “Lossless laser beam shaping”, JOSA A, 1996 的解析思想，设计一个 phase-only DOE 初始相位，并通过 Fourier lens 焦平面传播仿真验证矩形平顶光斑。

本版本重点是“解析相位 + 物理传播 + 诊断”，不是盲目运行 GS/WGS/MRAF 优化器。后续可以在解析相位基础上加入轻量 refinement，但默认不启用。

## 如何运行

在 MATLAB R2025a 中进入项目目录：

```matlab
cd('E:\program\Point2P\DOE_ROMERO_DICKEY_MATLAB')
run_one
```

每次运行会生成：

```text
artifacts/YYYYMMDD-HHMMSS/
```

其中包含配置、beta 报告、metrics 报告、相位数据、焦平面图像和中心剖面 CSV。

运行时会分两步验证：

1. `ideal_` 前缀：ideal Fourier-lens FFT baseline，假设 DOE 面就是 Fourier transform pupil plane，用于先验证 Romero-Dickey separable phase 本身。
2. `full_` 前缀：full propagation model，加入 DOE 到 lens 的 `200 mm` 自由空间传播、透镜面圆形口径、薄透镜二次相位、lens 到焦平面的传播，用于评估真实光路偏离理想 Fourier transform 条件后的质量变化。

注意：如果 full model 使用无限透镜口径，在理想标量传播中 DOE 到 lens 的自由空间传播主要表现为频域相位因子，焦平面强度可能与 ideal FFT baseline 几乎相同。因此本项目默认在 full model 中加入可配置的 lens aperture。

同一 artifacts 文件夹内还会生成：

```text
comparison_ideal_vs_full.txt
```

用于直接对比 FWHM、transition width、uniformity、efficiency、side-lobe 和 overshoot。

## 参数在哪里改

主要参数集中在：

```text
config/default_config.m
```

常用可调参数包括：

- `lambda_m = 532e-9`：波长。
- `f_m = 429e-3`：场镜 / Fourier lens 焦距。
- `doe_to_lens_m = 200e-3`：DOE 到镜头距离；默认 Fourier 模型只记录该参数，不参与传播。
- `aperture_diameter_m = 15e-3`：DOE 圆形 clear aperture。
- `lens_aperture_diameter_m = 15e-3`：full propagation 中透镜面的圆形口径；默认等于 DOE clear aperture，用于体现 200 mm 传播后的实际截断影响。
- `apply_lens_aperture_in_full_model = true`：full propagation 是否应用透镜口径。
- `input_1e2_diameter_m = 5e-3`：入射高斯光强 `1/e^2` 直径。
- `target_size_x_m = 330e-6`、`target_size_y_m = 120e-6`：目标矩形 50% 强度边界尺寸。
- `N = 2048`：FFT 网格数。
- `requested_focal_dx_m = 2.5e-6`：期望焦平面采样。
- `phase_sign = 1`：相位符号约定，可改为 `-1` 对比 FFT/Fresnel 符号。
- `phase_scale_x`、`phase_scale_y`：x/y 方向相位缩放，用于后续调参。
- `target_edge = "hard"` 或 `"logistic"`：目标图样边缘模型。

## 重要采样约定

默认焦平面采样由以下关系确定：

```text
focal_dx = lambda * f / (N * dx_doe)
```

为了得到约 `2.5 um/pixel`，本项目默认 DOE 数值计算窗口约为 `91.2912 mm`，对应 `dx_doe ≈ 44.5758 um`。

这并不表示光斑或 DOE aperture 是 91 mm。实际物理含义是：

- `5 mm` 是入射高斯光斑 `1/e^2` 强度直径。
- `15 mm` 是 DOE clear aperture 圆形通光孔径。
- `91.2912 mm` 是 FFT 数值窗口，用于满足焦平面采样要求。

## 主要输出看什么

- `ideal_beta_report.txt`、`full_beta_report.txt`：检查 `ri`、`Ro_x/Ro_y`、`beta_x/beta_y` 和 beta 判断说明。
- `ideal_metrics.txt`、`full_metrics.txt`：分别检查两种模型的 FWHM、transition width、uniformity、efficiency、side-lobe、overshoot、undershoot。
- `ideal_phase_wrapped.png`、`full_phase_wrapped.png`：检查 DOE wrapped phase；两者相同，只是随各模型结果保存。
- `ideal_focal_intensity.png`、`full_focal_intensity.png`：两种传播模型的焦平面归一化强度。
- `ideal_focal_intensity_log.png`、`full_focal_intensity_log.png`：对数强度，用于看裙边和旁瓣。
- `ideal_roi_intensity.png`、`full_roi_intensity.png`：目标附近 ROI 强度。
- `ideal_center_profiles.png`、`full_center_profiles.png`：x/y 中心线、目标剖面、50%、13%、90% crossing 和目标边界。
- `ideal_x_profile.csv`、`full_x_profile.csv` 等：中心剖面数据。

## Romero-Dickey 与 GS/WGS/MRAF 的区别

Romero-Dickey 方法给出的是基于 stationary phase / geometrical optics 近似的解析 phase-only 初始设计。它强调能量映射和无损相位元件，而不是通过迭代强行拟合目标强度。

GS/WGS/MRAF 是迭代相位恢复 / 权重更新算法，通常更适合在给定采样、孔径、相位约束和目标 ROI 后做数值优化。但它们容易依赖初值、权重、ROI 策略和迭代参数。

本项目第一版先跑通解析设计，再用传播和诊断判断问题难度。若后续加入 WGS，应该以解析相位为初值，并保留解析结果用于对比。

## 为什么 beta 重要

论文中的核心无量纲参数为：

```text
beta = 2*pi*ri*Ro/(lambda*f)
```

其中 `ri` 是输入强度下降到 `1/e` peak 的半径，不是 `1/e^2` 强度半径。本项目输入为 `5 mm 1/e^2` 强度直径，所以：

```text
w = 2.5 mm
ri = w / sqrt(2)
```

默认矩形目标采用 separable approximation：

```text
Ro_x = 330 um / sqrt(pi)
Ro_y = 120 um / sqrt(pi)
```

因此默认约为：

```text
beta_x ≈ 9.061
beta_y ≈ 3.295
```

`beta_y` 明显更小，说明 y 方向更难获得锐利边缘和良好均匀性。

## 当前局限

- 这是解析 / separable 初始相位，不保证直接得到工业级完美矩形平顶。
- 有限 beta 会导致裙边变宽、边缘不够锐利。
- y 方向目标只有 `120 um`，`beta_y` 较小，边缘和均匀性更困难。
- 硬边 target 对频域高频要求高，容易产生 ringing。
- 15 mm aperture 截断会影响边缘 wiggle 和旁瓣。
- 默认传播假设 DOE 位于有效 pupil / Fourier transform plane。
- 当前不包含制造量化、相位深度、材料色散、SLM LUT、相机反馈或闭环实验。

## 下一步调参建议

优先尝试：

1. 改 `phase_sign = -1` 对比符号约定。
2. 小范围调 `phase_scale_x`、`phase_scale_y`，尤其 y 方向。
3. 将 `target_edge` 改为 `"logistic"` 并调整 `transition_width_x_m/y_m`。
4. 若解析结果合理但均匀性不足，再考虑加入少量 WGS refinement。
