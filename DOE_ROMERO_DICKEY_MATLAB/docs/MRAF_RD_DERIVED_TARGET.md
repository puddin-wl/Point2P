# MRAF RD-derived Target 设计说明

本项目的 MRAF refinement 默认不使用人工 hard rectangle 或 logistic rectangle target。Romero-Dickey / Point2P 初始相位已经给出一个接近真实可实现的圆角矩形平顶 envelope；MRAF 只应做轻量 refinement，而不是把焦斑硬压成数学矩形。

## 为什么不用硬 core/edge 拼接

旧版 target 采用分段逻辑：`core_mask` 内设常数，`edge_mask` 内设 `I_rd_norm^gamma`，`noise_mask` 内设 `NaN`。这种 core/edge 边界会产生振幅斜率不连续或 kink。相位-only MRAF 无法直接满足这种硬拼接，容易把误差变成平台边缘后的 shoulder：平台后先抬升/隆起，然后才下降。

因此新版默认 target 使用连续 smooth core-edge blend，优先消除 shoulder，再考虑 transition 是否变窄。

## Smooth RD-derived target

先用 RD phase 做 ideal FFT：

1. `U_rd = forward(input_amp .* exp(1i * phase_rd))`
2. `I_rd = abs(U_rd).^2`
3. 用目标核心区 mean intensity 归一化，得到 `I_rd_norm`

然后生成连续 blend map：

- `flat_high = 0.95`：`I_rd_norm >= flat_high` 时 `blend = 1`
- `flat_low = 0.75`：`I_rd_norm <= flat_low` 时 `blend = 0`
- 中间区域用 smoothstep：

```text
t = (I_rd_norm - flat_low) / (flat_high - flat_low)
blend = t^2 * (3 - 2*t)
```

目标强度定义为：

```text
I_env = min(I_rd_norm, 1)^gamma
I_flat = 1
I_target = blend * I_flat + (1 - blend) * I_env
A_target = sqrt(max(I_target, 0))
```

默认 `gamma = 1.0`，即先不 sharpen。第一目标是避免 shoulder 和边缘能量堆积，而不是追求更窄 transition。

## target 是 amplitude，不是 intensity

GS / WGS / MRAF 的 farfield projection 替换的是复光场振幅：

`Uf_projected = A_projected .* exp(1i * angle(Uf))`

因此主 target 是 `A_target`。输出的 `rd_derived_target_intensity.png` 只是 `A_target.^2` 的诊断图。

## NaN 表示 free/noise region

`A_target` 中的 `NaN` 不表示目标强度为 0，而表示该区域不强制投影到固定振幅。MRAF 在这些区域保留当前 farfield 的一部分：

`Uf_proj(noise_mask) = mraf_factor * Uf(noise_mask)`

新版默认 `mraf_factor = 0.7`，更保守；如果 shoulder 仍上升，可以尝试 `mraf_factor = 1.0`，即 free region 完全保留。

## mask 约定

mask 仍用于 metrics 和 projection 区域定义，但 target amplitude 不再用 core/edge 硬赋值：

- `core_mask`: `I_rd_norm >= 0.85`
- `signal_mask`: `I_rd_norm >= 0.50`
- `edge_mask`: `0.10 <= I_rd_norm < 0.85`
- `noise_mask`: `I_rd_norm < 0.05` 或 guard 外
- `guard_mask`: 包围主光斑的安全区域

## shoulder_peak 指标

沿中心 x/y profile，用 core mean 归一化后，在 50% 半宽附近的内侧边缘带统计肩部峰值：

- x/y 边缘带：`0.70*half50` 到 `1.05*half50`，左右/上下两侧都算
- `shoulder_peak_x = max(inner_edge_band_x) - 1`
- `shoulder_peak_y = max(inner_edge_band_y) - 1`

同时输出 left/right/top/bottom。MRAF 的 `shoulder_peak_x/y` 不应高于 RD baseline；如果 shoulder 明显上升，即使 transition 变窄也判为失败。

## FFT 和 full propagation

MRAF 闭环只使用 ideal FFT：

- forward: `Uf = fftshift(fft2(fftshift(U)))/N`
- backward: `U = ifftshift(ifft2(ifftshift(Uf)))*N`

焦平面采样：`dx_f = lambda*f/(N*dx_doe)`。当前 `N=2048`、`dx_f=2.5 um/pixel`，DOE 计算窗口由该关系反推。`15 mm` 是窗口内圆形 clear aperture，`5 mm` 是 Gaussian intensity `1/e^2` diameter。

Full propagation 不参与 MRAF 闭环，只用于验证：DOE 面 -> 自由空间 200 mm -> lens pupil/quadratic phase -> lens 到后焦面 429 mm 的单 FFT Fresnel 传播。报告中会明确写出 full propagation model；不允许把 ideal FFT 结果重复标成 full propagation。
## 当前主线 baseline：caseC_stable

当前 MRAF 主线已经冻结为 `caseC_stable`：

- `target_mode = "rd_derived"`
- `rd_target_gamma = 1.0`
- `mraf_factor = 1.0`
- `phase_blend = 0.25`
- `n_iter = 10`
- `method = "GS-MRAF"`
- `use_wgs = false` / `wgs_exponent = 0`

这个 baseline 的选择依据是：`mraf_factor = 1.0` 完全保留 free/noise region，显著降低 x 方向 shoulder，同时改善 core RMS，并保持 size50 基本稳定。后续不应默认使用 `gamma > 1`、WGS 或大迭代数。

## 三个主标注指标

图和 metrics 中统一显示以下三个主指标：

1. 匀化光斑尺寸（50%）：中心剖面在 0.50 intensity level 的左右 crossing 距离，变量为 `size50_x_um` / `size50_y_um`。
2. 匀化光斑尺寸（13.5%）：中心剖面在 0.135 intensity level 的左右 crossing 距离，变量为 `size13p5_x_um` / `size13p5_y_um`。
3. 传输区宽度（13.5%-90%）：左右两侧分别计算 90% crossing 到 13.5% crossing 的距离，再取平均，变量为 `transition_13p5_90_x_um` / `transition_13p5_90_y_um`。

13.5% 指标用于描述外层 skirt / envelope 尺寸。50% 尺寸更接近主匀化平台大小，而 13.5% 能反映边缘外侧能量包络是否扩张或收缩。

## 为什么先压 shoulder 而不是 transition

transition width 变窄不一定代表光斑更好。如果平台边缘出现 shoulder，说明能量在平台后方堆积，视觉形状和实际匀化质量都会变差。当前调参优先级是：先确保 `shoulder_peak_x/y` 不高于 RD baseline，再看 core RMS、size50、side lobe 和 transition。

## gamma probe 的边界

`gamma` 只作用在 RD envelope 项：`I_env = min(I_rd_norm, 1)^gamma`，再通过 smoothstep `blend_map` 与 `I_flat=1` 连续混合。probe 不允许恢复旧的 core/edge 硬拼接，也不允许重新启用 WGS 或 hard rectangle target。

当前只允许 `gamma <= 1.08` 的小幅试探；如果 shoulder 或 y 方向 size50 变差，即使 transition 变窄也不接受。

## pivot50 remap

普通 `gamma > 1` 会固定 `I=0` 和 `I=1`，但不会固定 `I=0.5`。因此它虽然能略微压窄 transition，却会把 50% crossing 往内拉，导致 `size50_x/y` 缩小。

`target_edge_mode="pivot50"` 用 50% 强度作为锚点改变边缘斜率：

```text
I < 0.5:   I_env = 0.5 * (2*I)^outer_power
I >= 0.5:  I_env = 0.5 + 0.5 * (2*I - 1)^inner_power
```

它保证 `I=0 -> 0`、`I=0.5 -> 0.5`、`I=1 -> 1`。`outer_power > 1` 压低 50% 外侧低强度 skirt；`inner_power < 1` 轻微抬高 50% 到平顶之间的内侧边缘。随后仍通过 `flat_low=0.75` 到 `flat_high=0.95` 的 smoothstep blend 平滑连接到核心平顶，不允许回到 core/edge 硬拼接。

## Hard rectangle + few-step MRAF 对照

`target_mode="hard_rectangle"` 是一个对照实验，不替代当前 RD-derived / pivot50 主线。它使用标准矩形信号区：`|x|<=165 um`、`|y|<=60 um`，矩形内 `target amplitude=1`，矩形外仍然是 MRAF free region (`NaN`)。因此它不是“外侧强制为 0”的全平面硬矩形 GS，而是“硬矩形信号区 + 自由外侧”的少步 MRAF。

这个实验的目的只是测试理想矩形目标是否能带来更窄 transition。判断时必须同时检查 shoulder、side lobe、size50 和 core RMS；不能只按 transition width 下结论。
