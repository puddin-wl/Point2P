# truncated_beam_phase — 截断光束 DOE 相位设计

针对实际实验中「扩束 25mm → 5mm 孔径截断 → 10cm → DOE（f=300mm）」的光路，设计与之匹配的 DOE 相位。

## 物理问题

原始 DOE 相位是为**原生 5mm 高斯光束**设计的。在实际光路中，激光先经过 10× 扩束（扩束后约 25mm 直径），再用 5mm 硬边孔径截取中心部分，传播 10cm 后到达 DOE 面。此时 DOE 面的振幅分布与原生 5mm 高斯有本质区别：

- **孔径内强度接近均匀**（25mm 高斯中心 5mm 内变化 ~8%）
- **硬边截断产生衍射环**（经 100mm 传播后到达 DOE 面）
- **振幅呈圆形分布**（非可分离的 x/y 高斯）

用错误的输入场假设去跑 DOE 相位，焦面 RMS 非均匀性高达 **46.7%**。

## 解决方案

### 1. 数值初始相位 (`compute_numerical_phase.py`)

Romero-Dickey 解析公式 $\phi(\xi)=\beta[\xi\frac{\sqrt\pi}{2}\mathrm{erf}(\xi)+\frac12 e^{-\xi^2}-\frac12]$ 硬编码了高斯输入 $\exp(-\xi^2)$。对于截断光束，改用**数值能量守恒法**计算初始相位：

**数学原理**（稳相近似 + 能量守恒）：

$$\frac{d\phi}{dx} = \frac{2\pi}{\lambda f} \cdot u(x), \quad CDF_{in}(x) = CDF_{out}(u)$$

其中 $CDF(z) = \frac{\int_{-\infty}^z I(t)\,dt}{\int_{-\infty}^\infty I(t)\,dt}$

**步骤**：
1. 计算 DOE 面的截断光束复振幅（大高斯 → 5mm 硬边孔径 → 角谱传 100mm）
2. 提取 x/y 方向的**边际强度分布** $I_{\text{marg},x}(x) = \int I(x,y)\,dy$
3. 构建 RTAD 升余弦平顶输出目标
4. CDF 匹配 → 映射函数 $u(x)$ / $v(y)$
5. 数值积分 → $\phi_x(x)$ / $\phi_y(y)$
6. 构建可分离 2D 相位 $\phi(x,y) = \phi_x(x) + \phi_y(y)$

由边际分布计算初始相位后，焦面 RMS 从 46.7% 改善到 **39.6%**。

### 2. WGS 精化 (`run_wgs_refinement.py`)

将截断光束的真实复振幅 $|E_{\text{DOE}}| \cdot e^{i\phi_{\text{diff}}}$ 作为 `input_amp`，数值相位 $\phi_{\text{numerical}} + \phi_{\text{diff}}$ 作为 `phase0`，调用现有 WGS 管线进行 200 次迭代优化。

WGS 不依赖于可分离假设（使用 2D FFT），可以突破数值初始相位的限制，将 RMS 从 39.6% 进一步降低。

## 结果

| 阶段 | RMS 非均匀性 | 说明 |
|---|---|---|
| 原始相位（为高斯设计） | **46.7%** | 输入场假设完全错误 |
| 数值初始相位 | **39.6%** | 能量守恒可分离近似 |
| WGS 200 次精化 | **3.38%** | 2D FFT 迭代突破可分离限制 |

WGS 精化后焦面指标：

| 指标 | 值 | 目标 |
|---|---|---|
| RMS 非均匀性 | 3.38% | < 2% |
| size50_x | 325.0 μm | 330 μm |
| size50_y | 127.4 μm | 120 μm |
| e⁻² 效率 | 80.0% | — |

## 文件结构

```
truncated_beam_phase/
├── README.md                       # 本文件
├── compute_numerical_phase.py      # 步骤1: 边际分布 CDF 匹配 → 数值初始相位
├── run_wgs_refinement.py           # 步骤2: 截断场 + 数值相位 → WGS 精化
├── plot_diagnostics.py             # 诊断: 焦面强度图 + 中心剖面
├── artifacts/
│   ├── *_numerical_phase/          # 数值相位输出
│   │   ├── phase_numerical.npy     # 数值初始相位 (wrapped [0, 2π))
│   │   ├── amplitude_at_doe.npy    # DOE 面振幅 (截断光束)
│   │   ├── field_at_doe.npy        # DOE 面复振幅
│   │   └── ...
│   └── *_wgs_refined/              # WGS 精化输出
│       ├── phase_refined.npy       # 最终 DOE 相位
│       ├── reconstruction_refined.npy  # 焦面强度
│       ├── diagnostics_summary.png # 诊断图
│       └── ...
└── run_aperture_simulation.py      # 位于 ../lab_test_f300mm/: 截断光束前向模拟
```

## 运行

```bash
# 环境: conda activate slmrtad (CuPy 13.6 + NumPy 1.26, GPU RTX 5070 Ti)

# 步骤1: 计算数值初始相位
python compute_numerical_phase.py

# 步骤2: WGS 精化 (200 次迭代，GPU ~4秒)
python run_wgs_refinement.py --iters 200

# 诊断可视化
python plot_diagnostics.py
```

## 参数

| 参数 | 值 |
|---|---|
| 波长 λ | 532 nm |
| 焦距 f | 300 mm |
| 扩束后光斑 D_exp | 25 mm (1/e² 强度直径) |
| 孔径 | 5 mm (硬边圆形) |
| 孔径到 DOE 距离 z | 100 mm |
| 目标尺寸 | 330 × 120 μm |
| 网格 | 2048 × 2048 |
| DOE 像素 | 31.2 μm |
| 焦面像素 | 2.5 μm |

## 依赖

- `lab_test_f300mm/src/` — DOE 传播、WGS 优化、MAT 文件读写
- `real_world_simulation/src/` — 诊断指标、可视化
- `lab_test_f300mm/run_aperture_simulation.py` — 截断光束前向模拟（角谱传播 + 孔径截断函数）

## 作者

2026-05-08
