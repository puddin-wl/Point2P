# 当前局限

## 解析初始相位不等于最终工业 DOE

本项目第一版实现的是 Romero-Dickey 解析 / 半解析 separable 初始相位。它适合作为理解物理尺度、beta 难度和初始相位结构的工具，但不保证直接得到完美工业矩形平顶。

若需要更高均匀性或更严格效率，通常还需要在此基础上加入 WGS、MRAF 或其他约束优化，并结合真实器件、相机和反馈闭环。

## 有限 beta 的影响

Romero-Dickey stationary-phase 近似在 beta 较大时更接近几何光学能量映射。有限 beta 会造成：

- flat-top 边缘变宽。
- 裙边增强。
- 目标内部均匀性下降。
- 对 aperture 和采样更敏感。

本项目默认 `beta_y ≈ 3.295`，显著小于 `beta_x ≈ 9.061`，所以 y 方向更容易出现边缘钝化、振铃和均匀性不足。

## 硬边目标会引入高频需求

`target_edge = "hard"` 表示理想硬矩形。硬边包含很强的高空间频率需求，而有限 aperture、有限 beta 和有限采样都无法无限支持高频，因此可能出现 ringing 和 side-lobe。

可尝试将目标改为：

```matlab
cfg.target_edge = "logistic";
```

并调节 `transition_width_x_m`、`transition_width_y_m`，观察边缘和均匀性的折中。

## Aperture 截断

15 mm DOE clear aperture 会对输入高斯远尾进行截断。虽然默认 5 mm 光斑远小于 15 mm aperture，但 aperture 的有限支持仍会影响焦平面旁瓣、边缘 wiggle 和能量分布。

## 离焦和偏心未默认建模

项目现在同时输出 ideal baseline 和 full propagation。ideal baseline 假设 DOE 位于有效 pupil / Fourier transform plane，并且目标平面为后焦平面。full propagation 加入 DOE 到 lens 的 200 mm 自由空间传播、薄透镜相位和 lens 到焦平面的传播。

即便如此，实际实验中的：

- DOE 到 lens 传播。
- 透镜相位。
- 离焦 `dz`。
- 横向偏心。
- 入射波前倾斜。

都可能导致 overshoot、skirt broaden、中心偏移或线性相位倾斜。当前 full model 仍不包含真实像差、偏心、倾斜、有限 lens clear aperture 或相机响应。

## 未包含实验与制造因素

当前只是仿真设计，不包含：

- DOE 制造量化。
- 最大相位深度限制。
- 材料色散。
- 多波长效应。
- SLM 像素填充因子。
- SLM LUT。
- 相机响应。
- 闭环实验反馈。
- 光学系统像差。

因此仿真输出应作为解析设计和参数诊断依据，而不是最终实验性能保证。
