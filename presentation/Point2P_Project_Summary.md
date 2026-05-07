# Point2P: 矩形平顶 DOE 光束整形 — 项目综述

## 1. 项目目标

设计一片**衍射光学元件（DOE）**，将532nm高斯激光束（1/e²直径5mm）整形成焦平面上的**矩形平顶光斑**（330μm × 120μm），并评估该DOE在真实光路误差条件下的容差能力。

核心流程：**Point2P初始相位 → MRAF/WGS迭代优化 → 真实光路容差评估**。

---

## 2. 理论基础

### 2.1 初始相位：Romero-Dickey 点对点法

**参考文献**: L.A. Romero & F.M. Dickey, "Lossless laser beam shaping," *J. Opt. Soc. Am. A*, Vol. 13, No. 4, pp. 751–760 (1996).

#### 核心思想

在菲涅尔近似下，用稳相法（stationary phase）解析求解将高斯光束转换为平顶光束所需的DOE相位。解的质量取决于无量纲参数 **β**：

```
β = 2π · rᵢ · Rₒ / (λ · f)
```

其中 rᵢ 为输入1/e幅度半径，Rₒ 为输出特征尺寸，λ 为波长，f 为焦距。β 越大，几何光学近似越精确。

#### 可分离相位公式（本项目使用）

对于将高斯 `exp(-xi²/2)` 转为平顶的一维情况，Romero-Dickey相位为：

```
φ(xi) = β · [xi · √π/2 · erf(xi) + 1/2 · exp(-xi²) - 1/2]
```

其中 `xi = x / rᵢ`，`βₓ = 2π · rᵢ · Roₓ / (λf)`。

二维可分离相位 `φ(x,y) = φₓ(x) + φⵧ(y)`，Roₓ 和 Roⵧ 分别由目标尺寸确定。

#### 本项目参数

| 参数 | 值 | 说明 |
|---|---|---|
| λ | 532 nm | 激光波长 |
| f | 429 mm | 傅里叶透镜焦距 |
| 入射光斑 (1/e²强度直径) | 5 mm | `rᵢ = 2.5mm / √2 = 1.768mm` |
| 通光孔径 | 15 mm | 圆形 |
| 目标尺寸 (50%全宽) | 330 × 120 μm | Roₓ = 330/√π μm, Roⵧ = 120/√π μm |
| 网格大小 N | 2048 × 2048 | |
| 焦面采样 | 2.5 μm/pixel | |
| βₓ | 9.061 | |
| βⵧ | 3.295 | y 方向为限制方向 |
| DOE采样 | 44.576 μm/pixel | |

### 2.2 目标函数：RTAD (Raised-Cosine Target Amplitude Distribution)

**参考文献**: W. Chen et al., "Generation of high uniformity flat-top beams by reconstructing the amplitude distribution at descending edges," *Optics & Laser Technology*, 186, 112776 (2025).

#### 核心思想

传统平顶目标使用陡峭边缘的几何形状，在FFT迭代中会产生频谱泄露（Gibbs现象）和散斑噪声。RTAD方法在目标平顶区与背景之间引入**升余弦下降边缘（raised-cosine descending edge）**，将振幅跳变转换为平滑过渡区，抑制高频分量。

#### 目标构建公式

```
I_full(x,y) = C(|x|; a₀, a₁) · C(|y|; b₀, b₁)
A_full = √(I_full)
```

其中 C(u; u₀, u₁) 为升余弦边缘函数：

```
C(u) = 1                    (u ≤ u₀)
     = ½[1 + cos(π(u-u₀)/(u₁-u₁))]  (u₀ < u < u₁)
     = 0                    (u ≥ u₁)
```

参数：a₀ = a₅₀ - Δx, a₁ = a₅₀ + Δx，默认 Δx = 15μm, Δy = 8μm。

#### 截断RTAD约束模式

本项目使用截断RTAD（truncated RTAD）：只约束 `I_full ≥ exp(-2) ≈ 13.5%` 的信号区域，低强度尾巴释放到MRAF自由区。这避免了强迫DOE跟随完整数学尾迹到零，保留物理上有用的旁瓣/光晕结构。

信号区域划分：
```
mask_signal = I_full ≥ 0.135  OR  mask_flat
mask_free   = 保护窗口 AND NOT mask_signal
mask_bg     = NOT (mask_signal OR mask_free)
```

文献报告RTAD方法在仿真中达到98.8%均匀性，实验中91.0%，显著优于传统几何目标法。

### 2.3 优化算法

#### 迭代框架：GS / MRAF / WGS

三种算法共享相同的 FFT 迭代框架，区别仅在输出平面的振幅替换方式。

**GS (Gerchberg-Saxton)**：信号区强制替换为目标振幅，其余归零。约束过强，收敛不稳定。

**MRAF (Mixed-Region Amplitude Freedom)**：
- **参考文献**: C. Zhang et al., "Optimized holographic femtosecond laser patterning...," *Scientific Reports*, 6, 33281 (2016).
- 将输出平面划分为信号区、自由区（噪声区）。信号区替换为目标振幅，自由区保留当前场乘以 mraf_factor。远背景可单独衰减或保留。
- 本项目中 MRAF 投影公式：
  ```
  信号区:  E' = A_target · exp(i·arg(E))
  自由区:  E' = mraf_factor · E
  远背景:  E' = bg_factor · E        (bg_factor→1.0 = 不约束)
  ```

**WGS (Weighted Gerchberg-Saxton)**：
- **参考文献**: D.Y. Alsaka et al., "Dynamic flat-topped laser beam shaping method using mixed region amplitude freedom algorithm," *Applied Physics B*, 128, 137 (2022).
- 在 MRAF/GS 投影基础上，对平顶核心区（mask_flat）引入自适应权重反馈：
  ```
  weights ← weights · (mean(|E_flat|) / |E|)^exponent
  weights = clip(weights, min, max)
  weights /= mean(weights)
  ```
- 权重只在 mask_flat（无边缘过渡区的纯平顶）内更新。边缘、自由区、远背景的权重保持 1.0。
- WGS 反馈是纯平顶均匀性校正，不把 size50、效率等诊断指标作为约束。

#### 算法参考

实现参考 **slmsuite v0.4.1**（GitHub 开源 SLM 控制库），特别是 `holography/algorithms/_hologram.py` 和 `_feedback.py`。关键语义：
- 目标数组为目标振幅，非强度
- 输入和目标振幅均 L2 归一化
- FFT 传播使用中心化 FFT + `norm="ortho"`
- 本项目用显式 mask 替代 slmsuite 的 `NaN` 标记

---

### 2.4 第一阶段优化探索

第一阶段通过系统性的参数扫描，确定最终的 DOE 优化管线。核心流程如下：

```
phase0 (Point2P解析解) → direct flat-local WGS (不经过MRAF预热)
```

#### 探索路线图

| 轮次 | 探索参数 | 关键发现 |
|---|---|---|
| 1 | WGS weight_max 扫描 (1.0~3.0) | max=1.0 太限制，≥1.5 之后差别不大。暂选 1.5。 |
| 2 | WGS weight_min 扫描 (0.5, 0.7, 0.8) | min 对 y 方向几乎无影响，暂选 0.5。 |
| 3 | 迭代次数扫描 (100~300) | 100 次之后变化很小。选 200。 |
| 4 | 反馈指数扫描 (0.2~0.9) | 0.6-0.7 以上改善饱和。暂选 0.7。 |
| 5 | mraf_factor 扫描 (0.2~1.2) | 越大越自由，但 1.2 过度放开。选 1.0。 |
| 6 | release_level 扫描 (0.025~0.135) | 影响不大，可回到 13.5%。 |
| 7 | **bg_factor 扫描 (0～1.0)** | **最大发现：不衰减背景是最大单项改进。** |
| 8 | 基于无背景基线的 mraf_factor 重扫 (0.8, 1.0, 1.2) | mraf_factor=0.8 剖面更好。 |
| 9 | 基于无背景基线的 release_level 重扫 | 回到原始 13.5%。 |
| 10 | 基于无背景基线的反馈指数重扫 | 0.8 的肩部最低。 |
| 11 | 基于无背景基线的 weight_max 重扫 (1.5~3.0) | max=2.0 为工作点。 |

#### 核心经验

**不强约束远背景是最大的单项改进。** 早期对远背景做较强衰减（bg_factor=0.05）时，平台被背景约束牵制；随着 bg_factor → 1.0，远背景基本不被压低，平台均匀性和中心剖面观感明显改善。折衷是远背景能量增加，但换来的均匀性改善远大于代价。

**MRAF 预热没有必要。** 直接 flat-local WGS（不经 MRAF 预热）给出了最好的平顶均匀性。MRAF-only 的结果较差，MRAF→WGS 混合也没有超越纯 WGS。

#### 最终工作基线

```
方法:       direct WGS (不经过 MRAF)
WGS 策略:   flat_local (2D局部反馈)
迭代次数:   200
反馈指数:   0.8
权重范围:   [0.5, 2.0]
mraf_factor: 0.8      (在 WGS 的自由区投影中)
bg_mode:    attenuate
bg_factor:  1.0       (实质 = 不管背景)
release_level: 0.135  (exp(-2))
```

**基线结果** (`fixed_baseline_bg0p9_initial_compare_20260429-175148`)：

| 指标 | 值 |
|---|---|
| size50_x / size50_y | 330.2 / 123.6 μm |
| size13.5_x / size13.5_y | 350.7 / 145.6 μm |
| transition_13.5_90 (x/y) | 20.3 / 21.7 μm |
| RMS 非均匀性 | **1.86%** |
| e⁻² 衍射效率 | **92.5%** |

---

## 3. 项目架构

```
Point2P/
├── initial_phase_generation/     # MATLAB: Romero-Dickey初始相位生成
│   ├── generate_initial_phase.m  # 核心: 稳相法解析相位
│   ├── run_initial_phase_generation.m
│   └── artifacts/<timestamp>/    # 基线phase0.mat (已冻结)
├── rtad_mraf_gs_python/          # Python: MRAF/WGS迭代优化
│   ├── src/mraf_gs.py            # 核心: GS/MRAF/WGS精化循环
│   ├── src/rtad_target.py        # RTAD升余弦矩形目标
│   ├── src/propagation.py        # FFT传播 + 高斯入射场
│   ├── src/metrics.py            # 诊断指标计算
│   ├── src/backend.py            # NumPy/CuPy后端
│   ├── run_rtad_mraf_gs_case.py  # CLI入口
│   └── artifacts/<timestamp>/    # phase_refined.npy + 诊断
├── real_world_simulation/        # Python: 真实光路容差评估
│   ├── src/field_models.py       # 真实入射场 (偏移/发散/倾角等)
│   ├── src/propagation.py        # 传播 + 离焦
│   ├── src/metrics.py            # 诊断 + 容差指标
│   ├── src/plotting.py           # 强度图/剖面图/PDF
│   ├── run_real_world_sweep.py   # CLI入口
│   └── artifacts/                # 容差扫描结果
├── result_diagnostics/           # MATLAB: 焦平面诊断 (备用)
├── target/                       # MATLAB: RTAD目标生成 (备用)
└── text/                         # 参考文献
    ├── josaa-13-4-751.pdf        # Romero-Dickey (1996)
    ├── srep33281.pdf             # MRAF (2016)
    ├── s00340-022-07860-5.pdf    # MMRAF动态整形 (2022)
    └── 1-s2.0-S0030399225003640-main.pdf  # RTAD (2025)
```

---

## 4. 真实光路误差 — 容差分析

固定设计好的DOE相位，改变光路/入射条件，评估当前DOE对误差的感度。

模拟的误差类型：

1. **离焦 (defocus)**: 观察面偏离理想焦平面
2. **光束偏移 (beam_offset)**: 入射光斑偏离DOE中心
3. **光束尺寸 (beam_size)**: 入射光斑1/e²直径偏差
4. **发散角 (divergence)**: DOE面处波前曲率（≠光束倾斜）
5. **光束倾角 (pointing)**: 整束光倾斜 → 焦斑平移
6. **通光孔径 (aperture)**: 有效通光区域变化
7. **椭圆度 (ellipticity)**: x/y方向光斑尺寸不一致

### 5种典型光斑误差 — 拟合结果

| 问题 | 根因 | 机制 |
|---|---|---|
| ① 长边中间内凹 | **椭圆度 Dx > Dy** | x方向过度照明抢走长边中心处y向能量 |
| ② 长边能量缺失+短边倾斜 | beam_offset（光斑偏离中心） | 不对称照明造成一侧能量缺失 |
| ③ 四条边都内凹 | **负离焦 + 光束偏大** | 观察面在焦点前 + 扩束过大 → 四边对称内缩 |
| ④ 能量上下聚集 | **beam_size（光束尺寸偏差）** | 入射光斑尺寸改变有效照明轮廓 |
| ⑤ 长宽比不对 | 待进一步分析 | 可能与椭圆度/发散角相关 |

---

## 5. 关键数值

| 指标 | 名义值 | 说明 |
|---|---|---|
| 目标光斑尺寸 (50%) | 330 × 120 μm | |
| WGS优化后RMS非均匀性 | 1.86% | `mask_flat`内 |
| e⁻²衍射效率 | 92.5% | 进入13.5%轮廓内的能量比例 |
| 发散角容差 (RMS<5%) | ~±0.005 mrad | 非常敏感 |
| 离焦容差 (RMS<5%) | ~±0.75 mm | |
| 光束偏移容差 (RMS<10%) | ~±0.2 mm | |
| 椭圆度容差 (Dx/Dy偏差) | < 0.2mm | |

---

## 6. 参考文献

1. Romero, L.A. & Dickey, F.M. "Lossless laser beam shaping." *J. Opt. Soc. Am. A* 13(4), 751–760 (1996).
2. Zhang, C. et al. "Optimized holographic femtosecond laser patterning method towards rapid integration of high-quality functional devices in microchannels." *Scientific Reports* 6, 33281 (2016).
3. Alsaka, D.Y. et al. "Dynamic flat-topped laser beam shaping method using mixed region amplitude freedom algorithm." *Applied Physics B* 128, 137 (2022).
4. Chen, W. et al. "Generation of high uniformity flat-top beams by reconstructing the amplitude distribution at descending edges." *Optics & Laser Technology* 186, 112776 (2025).
5. slmsuite v0.4.1 — https://github.com/slmsuite/slmsuite
