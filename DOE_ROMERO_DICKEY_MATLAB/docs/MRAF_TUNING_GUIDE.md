# MRAF 调参指南

当前优先级已经调整：先消除 shoulder，再考虑 transition width。不要为了让 13-90 transition 变窄而接受平台边缘能量堆积。

默认保守参数：

- `n_iter = 10`
- `target_mode = "rd_derived"`
- `rd_target_gamma = 1.0`
- `mraf_factor = 0.7`
- `phase_blend = 0.25`
- `method = "GS-MRAF"`
- `wgs_exponent = 0.2`，且默认不启用 WGS

## shoulder 出现或升高

shoulder 是平台边缘后紧接着出现的抬升/隆起，不是真正的边缘变锐。处理顺序：

1. 确认使用 smoothstep blend target，而不是 hard core/edge 拼接。
2. 保持 `gamma = 1.0`，不要 sharpen。
3. 增大 `mraf_factor` 到 `1.0`，让 free/noise region 完全保留。
4. 降低 `phase_blend` 到 `0.1-0.25`。
5. 减少迭代数，例如 `10 -> 5`。
6. 不要开启 WGS；如果必须开启，先用 `wgs_exponent <= 0.2`。

如果 `mraf_factor = 1.0` 明显更稳定，说明 shoulder 主要来自 free/noise 区约束过强。

## 为什么默认 gamma 改成 1.0

`gamma > 1` 会让 RD envelope 的边缘更陡，可能缩窄 transition，但也会提高边缘能量堆积风险。当前阶段先用 `gamma = 1.0` 保留 RD envelope，只在核心做平滑拉平。等 shoulder 不再增加后，再尝试 `gamma = 1.1` 或 `1.2`。

## smoothstep blend 定义

核心平坦区和 RD envelope 之间使用连续 blend：

```text
flat_low = 0.75
flat_high = 0.95
t = (I_rd_norm - flat_low) / (flat_high - flat_low)
blend = t^2 * (3 - 2*t)
I_target = blend * 1 + (1 - blend) * min(I_rd_norm,1)^gamma
```

这样 core 到 edge 的 target 连续，避免硬拼接 kink 被 MRAF 转换成 shoulder。

## 判断一组结果是否可接受

按优先级看：

1. `shoulder_peak_x/y` 不高于 RD baseline，最好下降。
2. `core_rms` 不恶化。
3. `size_50_x/y` 不明显偏离 `330 x 120 um`。
4. `side_lobe_peak_x/y_rel_to_core` 不明显升高。
5. `efficiency_inside_signal/guard` 不明显下降。
6. transition 可以稍微变窄，但不是第一优先级。

如果 shoulder 上升，即使 transition 变窄也不接受。

## 三个小 case 的解释

- Case A `gamma=1.0, mraf_factor=0.7`：保守无 sharpen，用于确认 smooth target 是否解决硬拼接问题。
- Case B `gamma=1.1, mraf_factor=0.7`：轻微 sharpen；如果 shoulder 升高，说明 gamma sharpen 诱发肩部。
- Case C `gamma=1.0, mraf_factor=1.0`：free/noise region 完全保留；如果最稳定，说明 free region 投影过强是主要问题。

## transition 变窄但 shoulder/overshoot 上升

不要继续增加 gamma 或迭代。应降低 sharpen、增大 `mraf_factor`、降低 `phase_blend`，必要时减少迭代。

## core 变花或 RMS 恶化

降低更新强度：

1. `phase_blend -> 0.25` 或更低
2. `n_iter -> 5`
3. 关闭 WGS
4. 降低 WGS exponent 到 `0.2`

## 几乎无变化但形状稳定

先确认 shoulder 不增加，然后再尝试：

1. `n_iter = 20`
2. `gamma = 1.1`
3. `phase_blend = 0.5`

一次只改一个参数。

## 不要追求硬边

当前 y 方向 `beta_y` 较小，硬边目标通常会被有限 aperture 和相位-only 约束转换成 overshoot、shoulder 或 side lobe。MRAF 的主线应是保留 RD 原始圆角矩形形状，在不增加 shoulder 的前提下轻微改善核心平坦度。
## caseC_stable 主线和 n_iter 扫描

当前主线 baseline 是 `caseC_stable`：`gamma=1.0`、`mraf_factor=1.0`、`phase_blend=0.25`、`n_iter=10`、`GS-MRAF`、WGS off。它不是 aggressive sharpen，而是保留 RD envelope，并用 smooth target 做轻量核心修正。

后续只建议做单变量小扫描。当前允许的第一组扫描是：

- `n_iter = 5`
- `n_iter = 10`
- `n_iter = 15`

其他参数固定不变。判断顺序：

1. `shoulder_peak_x/y` 不增加，最好下降。
2. `core_rms` 不恶化，最好下降。
3. `size50_x/y` 接近 `330 x 120 um`。
4. `size13p5_x/y` 不出现明显外扩，作为 skirt / envelope 尺寸观察项。
5. `side_lobe_peak_x/y_rel_to_core` 不明显升高。
6. `transition_13p5_90_x/y` 可以适度变窄，但不能作为唯一目标。

图中新增三个中文标注指标：匀化光斑尺寸（50%）、匀化光斑尺寸（13.5%）、传输区宽度（13.5%-90%）。13.5% level 用来观察外层 skirt / envelope，能补充 50% 主平台尺寸无法反映的边缘扩展信息。

## transition_probe gamma 小试探

在 `caseC_iter10_baseline` 固定后，只允许先做非常小的 gamma 单变量试探：

- `gamma = 1.00`
- `gamma = 1.02`
- `gamma = 1.04`
- `gamma = 1.06`
- `gamma = 1.08`

固定参数必须保持：`n_iter=10`、`mraf_factor=1.0`、`phase_blend=0.25`、`GS-MRAF`、WGS off。这个 probe 的目标只是试探 `transition_13p5_90_x/y` 是否能轻微下降；如果 shoulder 或 size50 约束失败，就不能接受。

安全约束相对于 `gamma=1.00` baseline：

- `shoulder_peak_x/y` 增加不超过 `0.01`
- `core_rms` 增加不超过 `10%`
- `size50_x` 在 `330 +/- 5 um` 内，`size50_y` 在 `120 +/- 5 um` 内
- `side_lobe_peak_x/y` 增加不超过 `0.02`
- `efficiency_inside_signal` 下降不超过 `2%`

如果 `gamma=1.02~1.08` 都不能在不增加 shoulder / 不破坏 size50 的情况下缩窄 transition，不要继续加大 gamma。下一轮只考虑：phase_blend 微扫、flat_low 微扫，或单独设计外边缘 sharpen。

## pivot50 probe 调参优先级

当普通 gamma 导致 `size50` 收缩时，应优先尝试 `target_edge_mode="pivot50"`。它以 50% 强度点为锚点，目标是在 `size50` 基本不动、shoulder 不回来的前提下，轻微压缩 `13.5%-90%` 传输区。

推荐只试保守范围：

- `inner_power = 0.98 / 0.95 / 0.92 / 0.90`
- `outer_power = 1.02 / 1.05 / 1.08 / 1.10`

不要超过 `inner_power < 0.90` 或 `outer_power > 1.10`。判断顺序：`size50` 不动、shoulder 不回来、transition 轻微变窄、side lobe 可轻微存在。

## hard_rectangle 对照实验

如果要回头测试硬矩形，必须保持少量迭代，例如 `n_iter=3/5/8`，并保持外侧为 MRAF free region，不要把矩形外强压成 0。硬矩形通常会显著压窄 transition，但也可能明显增加 shoulder 和 core RMS。

硬矩形结果只作为对照：如果 transition 变窄的代价是平台后 shoulder 明显升高、size50 跑偏或 core RMS 恶化，就不应替代 RD-derived / pivot50 主线。
