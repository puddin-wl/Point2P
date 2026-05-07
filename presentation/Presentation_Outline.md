# Point2P 项目汇报大纲

> 基于 `Point2P_Project_Summary.md` 精简，按 PPT 页面组织

---

## Slide 1: 封面

- **标题**: 矩形平顶 DOE 光束整形 — 从设计到容差评估
- **副标题**: 基于 Romero-Dickey 点对点法 + RTAD 目标 + WGS 优化
- **日期**: 2026 年 5 月

---

## Slide 2: 项目目标

设计一片衍射光学元件（DOE），将 532nm 高斯激光束整形成焦平面上的矩形平顶光斑。

| 参数 | 值 |
|---|---|
| 激光波长 | 532 nm |
| 入射光斑 (1/e² 直径) | 5 mm |
| 傅里叶透镜焦距 | 429 mm |
| 通光孔径 | 15 mm |
| 目标光斑 (50% 全宽) | **330 × 120 μm** |
| 计算网格 | 2048 × 2048 |

**核心流程**: Point2P 初始相位 → WGS 迭代优化 → 真实光路容差评估

---

## Slide 3: 理论基础 — Romero-Dickey 点对点法

- **参考文献**: Romero & Dickey, "Lossless laser beam shaping," *JOSA A*, 13(4), 751–760 (1996)
- **核心思想**: 菲涅尔近似 + 稳相法，解析求解高斯→平顶的 DOE 相位
- **关键无量纲参数 β**: β = 2π · rᵢ · Rₒ / (λ · f)，β 越大几何光学近似越精确

**可分离相位公式** (本项目使用):
```
φ(x,y) = φₓ(x) + φⵧ(y)
φₓ(xi) = βₓ · [xi·√π/2·erf(xi) + ½·exp(-xi²) - ½]
```
其中 xi = x / rᵢ

**本项目 β 值**: βₓ = 9.06, βⵧ = 3.30（y 方向为限制方向）

---

## Slide 4: 理论基础 — RTAD 目标函数

- **参考文献**: Chen et al., "Generation of high uniformity flat-top beams...," *Optics & Laser Technology*, 186, 112776 (2025)
- **核心思想**: 用升余弦下降边缘替代陡峭几何边界，抑制 FFT 频谱泄露

**升余弦边缘**:
```
C(u) = 1                         (u ≤ u₀)
     = ½[1 + cos(π(u-u₀)/Δ)]     (u₀ < u < u₁)
     = 0                         (u ≥ u₁)
I_full(x,y) = C(|x|) · C(|y|)
```

**截断 RTAD**: 只约束 I_full ≥ 13.5% 的信号区，低强度尾迹释放到自由区
- 避免强迫 DOE 跟随完整数学尾迹
- 文献报告: 仿真 98.8% 均匀性，实验 91.0%

---

## Slide 5: 理论基础 — WGS 优化

- **参考文献**: Alsaka et al., "Dynamic flat-topped laser beam shaping...," *Applied Physics B*, 128, 137 (2022)
- **算法参考**: slmsuite v0.4.1（GitHub 开源 SLM 控制库）

**MRAF 投影**（自由区松弛约束，不给背景加硬边界）:
```
信号区:  E' = A_target · exp(i·arg(E))
自由区:  E' = mraf_factor · E
远背景:  E' = bg_factor · E
```

**WGS 权重反馈**（仅在 mask_flat 内更新）:
```
weights ← weights · (mean(|E_flat|) / |E|)^exponent
weights = clip(weights, 0.5, 2.0)
weights /= mean(weights)
```

**策略**: 2D flat_local（同时作用于 x 和 y）

---

## Slide 6: 第一阶段 — 参数优化探索

**核心流程**: `phase0 (Point2P) → direct flat-local WGS`

共 11 轮参数扫描:

| 轮次 | 参数 | 结论 |
|---|---|---|
| 1-2 | WGS 权重范围 | min=0.5, max=1.5 |
| 3 | 迭代次数 | 200 次足够 |
| 4 | 反馈指数 | 0.6-0.7 以上饱和 |
| 5-6 | mraf_factor / release_level | mraf=1.0, rel=13.5% |
| **7** | **背景衰减 bg_factor** | **最大发现** |
| 8-11 | 无背景基线重扫 | 微调到最终参数 |

---

## Slide 7: 第一阶段 — 核心发现

**不约束远背景是最大单项改进。**

- 早期 bg_factor=0.05 时，平台被背景约束牵制
- bg_factor → 1.0：远背景不再被压低，平台均匀性显著改善
- 折衷：远背景能量增加，但均匀性改善远大于代价

**MRAF 预热没有必要。** 直接 WGS 优于 MRAF-only 和 MRAF→WGS 混合。

---

## Slide 8: 第一阶段 — 最终基线

**工作参数**:
```
method = wgs, strategy = flat_local
iters = 200, feedback_exponent = 0.8
weight_min = 0.5, weight_max = 2.0
mraf_factor = 0.8, bg_factor = 1.0
release_level = 0.135
```

**基线结果** (`fixed_baseline_bg0p9_initial_compare_20260429-175148`):

| 指标 | 值 |
|---|---|
| size50_x / y | 330.2 / 123.6 μm |
| RMS 非均匀性 | **1.86%** |
| e⁻² 衍射效率 | **92.5%** |
| transition width (x/y) | 20.3 / 21.7 μm |

---

## Slide 9: 第二阶段 — 真实光路容差评估

**目标**: 不重新优化 DOE，固定最终相位，评估 DOE 对真实光路误差的敏感度

**模拟的 7 类误差**:

| 误差类型 | 物理含义 |
|---|---|
| defocus | 观察面偏离焦平面 |
| beam_offset (x/y) | 入射光斑偏离 DOE 中心 |
| beam_size | 入射光斑尺寸偏差 |
| divergence | DOE 面波前曲率 |
| pointing | 光束倾斜 → 焦斑平移 |
| aperture | 有效通光孔径 |
| ellipticity | x/y 光斑尺寸不一致 |

---

## Slide 10: 第二阶段 — 5 种典型光斑误差及根因

| # | 问题 | 根因 | 机制 |
|---|---|---|---|
| ① | **长边中间内凹** | 椭圆度 Dx > Dy | x 过度照明抢走长边中心 y 向能量 |
| ② | 长边能量缺失+短边倾斜 | beam_offset | 光斑偏离造成不对称照明 |
| ③ | **四条边都内凹** | 负离焦 + 光斑偏大 | 观察面在焦点前 + 扩束过大 |
| ④ | 能量上下聚集 | beam_size | 入射光斑尺寸改变照明轮廓 |
| ⑤ | 长宽比不对 | 待进一步分析 | 可能与椭圆度相关 |

---

## Slide 11: 容差量级参考

| 误差参数 | 容差量级 (RMS<5%) |
|---|---|
| 发散角 | ~±0.005 mrad（极敏感） |
| 离焦 | ~±0.75 mm |
| 光束偏移 | ~±0.2 mm |
| 椭圆度 | Dx/Dy 偏差 < 0.2 mm |

---

## Slide 12: 项目架构

```
Point2P/
├── initial_phase_generation/  [MATLAB] Romero-Dickey 初始相位
├── rtad_mraf_gs_python/       [Python] MRAF/WGS 迭代优化
├── real_world_simulation/     [Python] 真实光路容差评估
├── result_diagnostics/        [MATLAB] 焦平面诊断
├── target/                    [MATLAB] RTAD 目标
└── text/                      参考文献 (×4)
```

---

## Slide 13: 总结

1. **Point2P 初始相位**（Romero-Dickey 解析解）提供了接近目标的物理初值
2. **RTAD 目标函数**（升余弦下降边缘）有效抑制了频谱泄露
3. **直接 WGS 优化**（不约束背景）是最简洁高效的匀化方案
4. **容差分析**揭示了 5 类光斑误差的根因，可用于实验反向纠错
5. 最终 RMS 非均匀性 **1.86%**，衍射效率 **92.5%**

---

## Slide 14: 参考文献

1. Romero, L.A. & Dickey, F.M. "Lossless laser beam shaping." *JOSA A* 13(4), 751–760 (1996).
2. Zhang, C. et al. "Optimized holographic femtosecond laser patterning..." *Scientific Reports* 6, 33281 (2016).
3. Alsaka, D.Y. et al. "Dynamic flat-topped laser beam shaping..." *Applied Physics B* 128, 137 (2022).
4. Chen, W. et al. "Generation of high uniformity flat-top beams..." *Optics & Laser Technology* 186, 112776 (2025).
5. slmsuite v0.4.1 — github.com/slmsuite/slmsuite
