# Point2P — DOE 矩形平顶光斑设计：完整技术总结

> 将 532 nm 高斯激光束通过衍射光学元件（DOE）整形成焦平面上 330×120 μm 矩形平顶光斑。
> 核心管线：Romero-Dickey 初始相位 → RTAD 目标构建 → WGS 迭代精修 → 容差评估。

---

## 目录

1. [项目目标与物理参数](#1-项目目标与物理参数)
2. [项目架构](#2-项目架构)
3. [核心算法模块](#3-核心算法模块)
4. [理论文献](#4-理论文献)
5. [参数优化历程](#5-参数优化历程)
6. [关键数值结果](#6-关键数值结果)
7. [真实光路容差分析](#7-真实光路容差分析)
8. [核心经验教训](#8-核心经验教训)

---

## 1. 项目目标与物理参数

### 1.1 目标

设计一片**衍射光学元件（DOE）**，加载到空间光调制器（SLM）上，将 532 nm 高斯激光束整形成焦平面上的**矩形平顶光斑**，并评估该 DOE 在真实光路误差条件下的容差能力。

### 1.2 物理参数

| 参数 | 值 | 说明 |
|------|-----|------|
| 波长 λ | 532 nm | 绿光激光 |
| 焦距 f | 429 mm（正式）/ 100–300 mm（验证） | 傅里叶透镜 |
| 入射光斑 1/e² 直径 | 6 mm（2026-06 更新，原 5 mm） | 高斯光束 |
| 通光孔径 | 15 mm | 圆形硬边光阑 |
| 目标光斑尺寸 (50% 全宽) | **330 × 120 μm** | 矩形平顶 |
| 计算网格 N | 2048 × 2048 | 焦面采样 2.5 μm/pixel |
| SLM 规格 | 1024×1024 @ 17 μm 像素 | 目标 SLM |
| βₓ / β_y | 9.06 / 3.30 | y 方向为衍射限制方向 |

---

## 2. 项目架构

```
Point2P/
├── rtad_mraf_gs_python/          # ★ 正式管线 (f=429mm, Python/CuPy)
│   ├── src/
│   │   ├── mraf_gs.py            #   核心：GS / MRAF / WGS 迭代精修循环
│   │   ├── rtad_target.py        #   RTAD 升余弦边缘矩形目标生成
│   │   ├── propagation.py        #   中心化 FFT 传播 + 高斯入射场
│   │   ├── metrics.py            #   诊断指标（RMS, size50, 效率, 过渡宽度）
│   │   ├── diagnostics.py        #   综合诊断报告生成
│   │   ├── plotting.py           #   可视化（强度图/相位图/剖面对比）
│   │   ├── backend.py            #   NumPy / CuPy 后端抽象
│   │   ├── io_mat.py             #   MATLAB .mat 读写
│   │   └── utils.py              #   工具函数
│   ├── run_rtad_mraf_gs_case.py  #   CLI 入口（参数覆盖 config）
│   ├── run_pipeline.py           #   一键全流程（Stage 1→2→3）
│   ├── config_default.py         #   默认物理参数 + 精修参数
│   ├── make_phase0.py            #   纯 Python 生成 Romero-Dickey 初始相位
│   ├── convert_to_slm.py         #   2048² → SLM 原生分辨率转换（复振幅插值）
│   ├── add_blaze_grating.py      #   叠加闪耀光栅
│   └── artifacts/                #   输出：phase_refined.npy + 诊断图表
│
├── initial_phase_generation/     # Stage 1: Romero-Dickey 初始相位 (MATLAB)
│   ├── generate_initial_phase.m  #   核心：稳相法解析相位
│   ├── default_initial_phase_config.m
│   └── artifacts/20260428-141942/ #  基线 phase0.mat（已冻结）
│
├── real_world_simulation/        # Stage 3: 真实光路容差评估 (Python)
│   ├── src/
│   │   ├── field_models.py       #   真实入射场建模（偏移/发散/倾角等）
│   │   ├── propagation.py        #   传播 + 离焦
│   │   ├── metrics.py            #   容差诊断指标
│   │   └── plotting.py           #   强度图/剖面图/PDF
│   ├── run_real_world_sweep.py   #   CLI：7 类误差参数扫描
│   └── artifacts/                #   容差扫描结果
│
├── lab_test/                     # 实验室验证（非主管线）
│   ├── lab_test_f100mm/          #   f=100mm 焦距变体验证
│   ├── lab_test_f200mm/          #   f=200mm 标准验证
│   ├── lab_test_f300mm/          #   f=300mm 焦距变体验证
│   └── lab_wgs/                  #   实验反馈 WGS 在线优化 (MATLAB)
│
├── fig_analysis/                 # 捕获光斑图像分析
│   └── analyze_captured.py       #   梯度边缘法分析实验光斑照片
│
├── truncated_beam_phase/         # 截断光束数值相位计算
│   ├── compute_numerical_phase.py
│   └── run_wgs_refinement.py
│
├── result_diagnostics/           # 焦平面诊断 (MATLAB, 备用)
├── target/                       # RTAD 目标生成 (MATLAB, 备用)
│
├── presentation/                 # 项目汇报
│   ├── Point2P_Project_Report.pdf   # Beamer PDF
│   ├── Point2P_Project_Report.pptx  # PowerPoint
│   ├── Point2P_Project_Report.tex   # LaTeX 源码
│   └── Point2P_Project_Summary.md   # 项目综述
│
└── text/                         # 参考文献 (4篇)
    ├── josaa-13-4-751.pdf           # Romero & Dickey (1996)
    ├── srep33281.pdf                # Zhang et al. (2016) — MRAF
    ├── s00340-022-07860-5.pdf       # Alsaka et al. (2022) — WGS
    └── 1-s2.0-S0030399225003640-main.pdf  # Chen et al. (2025) — RTAD
```

### 2.1 核心管线流程

```
┌──────────────────────────────────────────────────────────────────┐
│                    Point2P 正式管线 (f = 429 mm)                    │
├──────────────────────────────────────────────────────────────────┤
│                                                                    │
│  Stage 1 (MATLAB)          Stage 2 (Python/CuPy)    Stage 3 (Python) │
│  ┌──────────────────┐     ┌────────────────────┐   ┌──────────────┐ │
│  │ Romero-Dickey     │     │ RTAD 目标构建       │   │ 7类误差扫描   │ │
│  │ 稳相法解析相位     │ ──→ │ + WGS 迭代精修      │ ─→│ 容差评估      │ │
│  │ → phase0.mat      │     │ → phase_refined.npy │   │ → 诊断报告    │ │
│  └──────────────────┘     └────────────────────┘   └──────────────┘ │
│                                                                    │
│  β_x=9.06, β_y=3.30      200 iters, flat_local    defocus, offset, │
│  基线已冻结                GPU ~3秒                 beam_size, ...  │
└──────────────────────────────────────────────────────────────────┘
```

---

## 3. 核心算法模块

### 3.1 初始相位生成（Romero-Dickey 点对点法）

**物理原理**：在菲涅尔近似下，用稳相法解析求解高斯→平顶能量守恒方程。

**一维相位公式**（可分离二维 = φₓ + φ_y）：
```
φ(xi) = β · [xi · √π/2 · erf(xi) + ½ · exp(−xi²) − ½]
其中 xi = x / rᵢ,  rᵢ = 光束 1/e 幅度半径
```

**β 参数**：`β = 2π · rᵢ · Rₒ / (λ · f)`

| β 范围 | 物理含义 |
|--------|---------|
| β > 20 | 几何光学区，稳相近似精确，初始相位接近真解 |
| 10 < β < 20 | 过渡区，WGS 小幅修正 |
| β < 10 | **衍射区**，WGS 明显改善均匀性 |

> 本项目 β_y = 3.30 处于衍射区，WGS 精修在 y 方向改善最显著。

### 3.2 RTAD 目标函数

**问题**：传统矩形目标用陡峭边缘 → FFT 频谱泄露（Gibbs 现象）→ 平台区波纹。

**RTAD 方法**：用**升余弦下降边缘（raised-cosine descending edge）**替代陡峭几何边界。

```
C(u) = 1                                         (u ≤ u₀)
     = ½[1 + cos(π(u−u₀)/(u₁−u₀))]               (u₀ < u < u₁)
     = 0                                         (u ≥ u₁)

I_full(x,y) = C(|x|; a₀, a₁) · C(|y|; b₀, b₁)
A_full = √(I_full)
```

**截断 RTAD（truncated RTAD）**：
- 只约束 `I_full ≥ exp(−2) ≈ 13.5%` 的信号区
- 低强度尾迹释放到 MRAF 自由区 → 避免强迫 DOE 跟随完整数学尾迹
- 参数：Δx = 15 μm, Δy = 8 μm, release_level = 0.135

**区域划分**：
```
mask_signal = I_full ≥ 13.5%  OR  mask_flat     ← 约束区
mask_free   = 保护窗口 AND NOT mask_signal       ← 自由区
mask_bg     = NOT (mask_signal OR mask_free)     ← 远背景
```

### 3.3 迭代精修算法：GS → MRAF → WGS

三种算法共享相同的 **FFT 迭代框架**（前向 FFT → 振幅替换 → 反向 FFT → 入射振幅约束），区别仅在输出平面的振幅替换方式。

#### GS（Gerchberg-Saxton）
```
信号区:  E' = A_target · exp(i·arg(E))
其余:    E' = 0
```
问题：约束过强，信号区外全部归零 → 散斑噪声严重（均匀性 ~37%）。

#### MRAF（Mixed-Region Amplitude Freedom）
```
信号区:  E' = A_target · exp(i·arg(E))
自由区:  E' = mraf_factor · E           ← 保留当前场，不强制归零
远背景:  E' = bg_factor · E
```
mraf_factor 控制自由区松弛程度：越小越自由，越大越接近 GS。

#### WGS（Weighted Gerchberg-Saxton）
在 MRAF 投影的基础上，对**平顶核心区（mask_flat）**引入自适应权重反馈：
```
weights ← weights · (mean(|E_flat|) / |E|)^exponent
weights = clip(weights, weight_min, weight_max)
weights /= mean(weights)                  ← 归一化防发散
```
权重仅在 mask_flat 内更新；边缘过渡区、自由区、远背景的权重保持 1.0。

**本项目策略**：`flat_local` — 2D 局部权重同时作用于 x 和 y，优于 `xy_then_x`（后者过早冻结 Y 方向）。

### 3.4 SLM 相位转换

将 2048×2048 计算相位转为 SLM 原生分辨率（1024×1024, 17 μm 像素）：

```python
# 正确：复振幅插值
c = np.exp(1j * phase_cropped)
phase_slm = np.arctan2(zoom(c.imag, ...), zoom(c.real, ...))

# 错误：直接对包裹相位插值 → 2π→0 跳变 → cubic spline 振铃
```

### 3.5 光斑图像分析（梯度边缘法）

实验光斑照片分析使用**梯度边缘法**定位平顶边界，用边界内中值强度做 flat-level。

> 避免了峰值阈值法因中心 hotspot 导致 flat-level 偏高的问题。

---

## 4. 理论文献

### 4.1 Romero & Dickey (1996) — 理论基石
**"Lossless laser beam shaping," *J. Opt. Soc. Am. A*, 13(4), 751–760**

- **贡献**：首次用稳相法解析求解高斯→平顶的 DOE 相位
- **关键公式**：一维 Romero-Dickey 相位积分 → 可分离二维相位
- **β 参数**：量化了几何光学近似的精度；β > 10 接近几何光学区
- **项目角色**：提供 Stage 1 初始相位 `phase0.mat`

### 4.2 Zhang et al. (2016) — MRAF 算法
**"Optimized holographic femtosecond laser patterning...," *Scientific Reports*, 6, 33281**

- **贡献**：提出 MRAF 算法——将输出平面划分为信号区和噪声区，噪声区保留当前场（×M）
- **与 GS 对比**：MRAF 均匀性 ~100% vs GS ~37%；能量利用率 4% vs 66%
- **M 参数**：M=1→GS, M=0.5→最优表面质量；在均匀性和能量效率之间折中
- **实验验证**：飞秒激光双光子聚合，单次曝光 240ms 制备 200×200μm 微结构阵列
- **项目角色**：提供自由区松弛概念（mraf_factor、bg_factor）

### 4.3 Alsaka et al. (2022) — WGS 权重反馈
**"Dynamic flat-topped laser beam shaping method using mixed region amplitude freedom algorithm," *Applied Physics B*, 128, 137**

- **贡献**：在 MRAF 基础上引入逐像素自适应 2D 权重反馈
- **权重公式**：`weights ← weights · (I_target/I)^exponent`，clip + normalize
- **物理含义**：暗区（|E| 小）→ 权重增大 → 下一轮给更多能量；亮区 → 权重减小
- **项目角色**：提供 WGS 精修核心算法（flat_local 策略）

### 4.4 Chen et al. (2025) — RTAD 目标
**"Generation of high uniformity flat-top beams by reconstructing the amplitude distribution at descending edges," *Optics & Laser Technology*, 186, 112776**

- **贡献**：用升余弦下降边缘替代陡峭几何边界，抑制 FFT 频谱泄露
- **截断 RTAD**：只约束 I ≥ 13.5% 信号区，尾迹释放到自由区
- **文献结果**：仿真 98.8% 均匀性，实验 91.0%
- **项目角色**：提供目标函数构建方法

### 四篇论文的协同关系

```
Romero-Dickey (1996)              Chen et al. (2025)
  稳相法解析初始相位  ────────────→  RTAD 升余弦目标
       │                                  │
       └──────────────┬───────────────────┘
                      ↓
               Zhang et al. (2016)
               MRAF 自由区松弛
                      │
                      ↓
               Alsaka et al. (2022)
               WGS 自适应权重反馈
                      │
                      ↓
               Point2P 完整管线
        phase0 → RTAD target → WGS refinement
```

> **一句话**：Romero-Dickey 给出物理上有意义的起点，RTAD 抑制了频谱泄露，MRAF 松弛了背景约束，WGS 在平顶核心区做逐像素匀化——四者叠加，从 RMS > 10% 收敛到 **1.86%**。

---

## 5. 参数优化历程

第一阶段通过 11 轮系统性参数扫描确定最终管线。核心发现：

### 5.1 扫描路线

| 轮次 | 参数 | 关键结论 |
|------|------|---------|
| 1 | WGS weight_max (1.0–3.0) | max=1.0 太限制，≥1.5 差别不大 |
| 2 | WGS weight_min (0.5, 0.7, 0.8) | 对 y 方向几乎无影响 |
| 3 | 迭代次数 (100–300) | 100 次后改进很小，选 200 |
| 4 | 反馈指数 (0.2–0.9) | 0.6–0.7 改善饱和 |
| 5 | mraf_factor (0.2–1.2) | 越大越自由，1.2 过度放开 |
| 6 | release_level (0.025–0.135) | 影响不大 |
| **7** | **bg_factor (0–1.0)** | **★ 最大发现：不衰减背景是最大单项改进** |
| 8–11 | 无背景基线重扫各参数 | 微调到最终参数 |

### 5.2 最终工作基线

```yaml
方法:         wgs（直接 WGS，不经 MRAF 预热）
WGS 策略:     flat_local（2D 局部反馈）
迭代次数:     200
反馈指数:     0.8
权重范围:     [0.5, 2.0]
mraf_factor:  0.8
bg_factor:    0.9（接近 1.0，实质不管背景）
release_level: 0.135 (exp(−2))
```

---

## 6. 关键数值结果

### 6.1 f=429mm 正式管线最佳结果

| 指标 | 值 | 说明 |
|------|-----|------|
| RMS 非均匀性 | **1.86%** | mask_flat 内 |
| size50_x / size50_y | 330.2 / 123.6 μm | 目标 330×120 |
| size13.5_x / size13.5_y | 350.7 / 145.6 μm | 13.5% 轮廓 |
| transition_13.5_90 (x/y) | 20.3 / 21.7 μm | 边缘过渡宽度 |
| e⁻² 衍射效率 | **92.5%** | 进入 13.5% 轮廓内能量比例 |

### 6.2 短焦距验证（f=100/200/300mm）

全部条件：λ=532 nm, N=2048, focal_dx=2.5 μm, 通光孔径 15 mm, 目标 330×120 μm。

| 焦距 | 光束直径 | β_x / β_y | RMS | size50_x/y | e⁻² 效率 |
|------|---------|-----------|-----|------------|---------|
| 100 mm | 5.0 mm | 38.9 / 14.1 | 0.27% | 330.8/116.9 | 99.1% |
| 200 mm | 5.0 mm | 19.4 / 7.1 | 0.12% | 328.2/117.8 | 96.5% |
| 200 mm | **7.0 mm** | **27.2 / 9.9** | **0.18%** | 327.1/118.3 | 97.6% |
| 300 mm | 6.0 mm | 13.0 / 5.0 | 1.22% | 329.3/118.1 | 94.9% |

> β 越大 → 初始相位越精确 → WGS 修正幅度越小。f=100mm 时 β_x=38.9 已进入几何光学区，WGS 几乎不需要修正。

---

## 7. 真实光路容差分析

### 7.1 模拟的 7 类误差

| 误差类型 | 物理含义 |
|---------|---------|
| defocus | 观察面偏离理想焦平面 |
| beam_offset (x/y) | 入射光斑偏离 DOE 中心 |
| beam_size | 入射光斑 1/e² 直径偏差 |
| divergence | DOE 面处波前曲率（≠ 光束倾斜） |
| pointing | 光束倾斜 → 焦斑平移 |
| aperture | 有效通光孔径变化 |
| ellipticity | x/y 方向光斑尺寸不一致 |

### 7.2 5 种典型光斑误差及根因

| # | 现象 | 根因 | 物理机制 |
|---|------|------|---------|
| ① | **长边中间内凹** | 椭圆度 Dx > Dy | x 方向过度照明抢走长边中心处 y 向能量 |
| ② | 长边能量缺失 + 短边倾斜 | beam_offset | 光斑偏离中心造成不对称照明 |
| ③ | **四条边都内凹** | 负离焦 + 光斑偏大 | 观察面在焦点前 + 扩束过大 |
| ④ | 能量上下聚集 | beam_size | 入射光斑尺寸改变有效照明轮廓 |
| ⑤ | 长宽比不对 | 待进一步分析 | 可能与椭圆度/发散角相关 |

### 7.3 容差量级

| 误差参数 | 容差量级 (RMS < 5%) |
|---------|---------------------|
| 发散角 | ~±0.005 mrad（**极敏感**） |
| 离焦 | ~±0.75 mm |
| 光束偏移 | ~±0.2 mm |
| 椭圆度 | Dx/Dy 偏差 < 0.2 mm |

---

## 8. 核心经验教训

1. **不强压背景是最大单项改进** — `bg_factor` ≈ 1.0 时平台更平；早期用 0.05 强压背景反而拖累均匀性
2. **MRAF 预热不需要** — 直接 WGS 的均匀性优于 MRAF→WGS 两段式；MRAF-only 结果更差
3. **X-only 不需要** — flat_local 的 2D 局部权重已足够；引入 X-only 会过早冻结 Y，损害均匀性
4. **复振幅插值是必须的** — SLM 相位转换必须对 `exp(iφ)` 插值，不能直接对包裹相位插值
5. **仿真 WGS 是实验 WGS 的必要热启动** — 不能从随机相位直接做实验反馈；散斑状态下拍不到平顶
6. **WGS 的 RMS 在首次权重更新后会跳升** — 这是预期行为，之后缓慢下降
7. **GPU 加速效果显著** — RTX 5070 Ti 上 200 轮约 3 秒（2048×2048 网格）
8. **发散角是最敏感的误差参数** — ±0.005 mrad 即可使 RMS 超过 5%，远超其他误差的容差量级

---

## 附录：运行命令速查

```bash
# 环境
conda activate slmrtad

# Stage 2: 仿真 WGS 精修（从已有 phase0）
cd E:\program\Point2P\rtad_mraf_gs_python
python run_rtad_mraf_gs_case.py \
    --phase-mat ../initial_phase_generation/artifacts/20260428-141942/phase0.mat \
    --phase-var phase0_wrapped_rad \
    --method wgs \
    --wgs-strategy flat_local \
    --iters 200 \
    --wgs-feedback-exponent 0.8 \
    --wgs-weight-min 0.5 \
    --wgs-weight-max 2.0 \
    --bg-factor 0.9

# Stage 3: 容差扫描
cd E:\program\Point2P\real_world_simulation
python run_real_world_sweep.py

# 光斑图像分析
python fig_analysis/analyze_captured.py <image.mat/.bmp/.png/.tif> --pixel-um 3.45

# 诊断已有结果
python rtad_mraf_gs_python/run_diagnostics_case.py rtad_mraf_gs_python/artifacts/<输出目录>

# 冒烟测试（验证环境，不需要 phase 文件）
python run_rtad_mraf_gs_case.py --iters 20 --smoke-shape 256
```

---

*文档生成时间：2026-06-08 | 基于 E:\program\Point2P 项目文件及四篇理论文献*
