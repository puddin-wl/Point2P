# Point2P — DOE 平顶光斑设计管线

将高斯激光束整形成 330×120 μm 矩形平顶光斑。波长 532 nm，计算网格 2048×2048。

## 目录结构

```
Point2P/
├── initial_phase_generation/     # Stage 1: Romero-Dickey 初始相位生成 (MATLAB)
├── lab_test/
│   ├── lab_test_f200mm/          # 标准管线参考实现 (f=200mm)
│   ├── lab_test_f100mm/          # f=100mm 焦距变体
│   ├── lab_test_f300mm/          # f=300mm 焦距变体
│   └── lab_wgs/                  # Stage 5: 实验反馈 WGS 在线优化 (MATLAB)
├── rtad_mraf_gs_python/          # 原始基线 (f=429mm, Python)
├── real_world_simulation/        # 容差评估：固定相位、扫入射/光路误差
├── fig_analysis/                 # 捕获光斑图像分析（梯度边缘法）
├── result_diagnostics/           # 焦平面结果诊断 (MATLAB)
├── truncated_beam_phase/         # 截断光束相位模块
├── dmd_phase/                    # DMD 超像素相位转换 (192×192, 54.8µm)
├── presentation/                 # 项目汇报 Beamer PDF
└── target/                       # 目标图形定义
```

详细内容见各子目录的 README。

## 标准管线（f=200mm 已验证）

### Stage 1 — 生成初始相位

Romero-Dickey 解析法：在高斯光束与矩形平顶之间求解稳相能量守恒，得到可分离的 DOE 初始相位。

- MATLAB：`lab_test/lab_test_f200mm/matlab/run_initial_phase_generation`
- Python（无需 MATLAB）：`lab_test/lab_test_f200mm/make_phase0.py --beam 7`
- 输出：`phase0.mat`（2048×2048, [0, 2π)）

### Stage 2 — 仿真 WGS 精修

**直接从初始相位做 WGS，不需要 MRAF 预热，不需要 X-only 阶段。** 经过多轮参数扫描验证，`method=wgs` + `strategy=flat_local` 的简洁配置给出了最好的平台均匀性。

```bash
conda activate slmrtad
cd lab_test/lab_test_f200mm

python run_rtad_mraf_gs_case.py \
    --phase-mat make_phase0_output/phase0.mat \
    --phase-var phase0_wrapped_rad \
    --beam-diameter 7 \
    --method wgs \
    --wgs-strategy flat_local \
    --iters 200 \
    --wgs-feedback-exponent 0.8 \
    --wgs-weight-min 0.5 \
    --wgs-weight-max 1.5 \
    --bg-factor 0.9 \
    --no-swap-phase-xy
```

已验证的关键参数选择：
- **method=wgs**：直接 WGS，不加 MRAF 预热段。测试证明 MRAF→WGS 的两段式不带来额外收益
- **strategy=flat_local**：只对平顶核心区做 2D 局部权重反馈。不做 xy_then_x（先 XY 再 X-only），后者让 Y 方向过早冻结，反而损害均匀性
- **bg-factor=0.9**：极弱的背景衰减，实质接近"不管背景"。早期用 0.05 强压背景反而拖累平台均匀性；0.9 是多次扫描确认的最优值
- **feedback-exponent=0.8**：较强的逐轮修正力度，收敛更快
- **weight-max=1.5**：max≥2 没有进一步收益
- **no-swap-phase-xy**：Python 生成的 phase0 使用 numpy 惯例，不需要 x/y 转置

GPU（RTX 5070 Ti）上 200 轮约 3 秒。

### Stage 3 — SLM 相位转换

将 2048×2048 计算相位转为 SLM 原生分辨率（1920×1080, 6.4 μm 像素）。**关键教训**：必须用复振幅插值，不能直接对包裹相位做插值——2π→0 的跳变会导致 cubic spline 振铃。

```python
# 正确做法：对 exp(iφ) 的实部/虚部分别插值，再取 angle
c = np.exp(1j * phase_cropped)
phase_slm = np.arctan2(zoom(c.imag, ...), zoom(c.real, ...))
```

详见 `lab_test/lab_test_f200mm/README.md` Stage 3 章节。

### Stage 4 — 光斑图像分析

用梯度边缘法分析实验拍到的平顶光斑照片。梯度边缘法利用边缘处梯度峰值定位平顶边界，然后用边界内中值强度做 flat-level——避免了峰值阈值法因中心hotspot导致 flat-level 偏高的问题。

```bash
python fig_analysis/analyze_captured.py <image.mat/.bmp/.png/.tif> --pixel-um 3.45
```

### Stage 5 — 实验反馈 WGS（最新）

仿真 WGS 之后，用真实光路反馈跑实验 WGS，补偿仿真中未建模的误差（光路不对准、SLM 非线性、杂散光等）。这是纯 MATLAB 实现，放在 `lab_test/lab_wgs/` 中。

核心思路是混合场 WGS：前向传播走真实光路（SLM→相机→实测强度），反传走仿真 FFT（提供仿真相位 + 实测强度拼成焦面复振幅）。

**不能跳过仿真 WGS 直接从随机相位做实验 WGS**——散斑状态拍不到平顶，梯度边缘法找不到矩形亮区，WGS 没有有效反馈信号。仿真 WGS 是必须的热启动。

详见 `lab_test/lab_wgs/readme.md`。

## 各焦距结果汇总

全部条件：λ=532 nm, N=2048, focal_dx=2.5 μm, 通光孔径 15 mm, 目标 330×120 μm。

| 焦距 f | 光束 | β_x / β_y | RMS | size50_x | size50_y | e⁻² 效率 |
|--------|------|-----------|-----|----------|----------|----------|
| 100 mm | 5.0 mm | 38.9 / 14.1 | 0.27% | 330.8 | 116.9 | 99.1% |
| 200 mm | 4.5 mm | 17.5 / 6.4 | 0.25% | 328.9 | 118.1 | 96.8% |
| 200 mm | 5.0 mm | 19.4 / 7.1 | 0.12% | 328.2 | 117.8 | 96.5% |
| 200 mm | 5.5 mm | 21.4 / 7.8 | 0.13% | 326.5 | 118.1 | 95.7% |
| 200 mm | 6.0 mm | 23.3 / 8.5 | 0.09% | 326.9 | 118.1 | 95.7% |
| 200 mm | **7.0 mm** | **27.2 / 9.9** | **0.18%** | 327.1 | 118.3 | 97.6% |
| 300 mm | 6.0 mm | 13.0 / 5.0 | 1.22% | 329.3 | 118.1 | 94.9% |
| 429 mm | 5.0 mm | 9.1 / 3.3 | 1.86% | ~330 | ~124 | 92.5% |

β 越大，稳相近似越准确，初始相位越接近精确解。β > 10 进入几何光学区，WGS 修正幅度很小。β < 10 进入衍射区，WGS 明显改善均匀性。

## 重要经验

- **不要强压背景**：`bg-factor` 接近 1.0 时平台更平，远处背景的能量代价在可接受范围内
- **MRAF 不需要**：直接 WGS 的均匀性优于 MRAF→WGS 两段式，没必要用 MRAF
- **X-only 不需要**：flat_local 的 2D 局部权重已足够；引入 X-only 会过早冻结 Y，损害均匀性
- **复振幅插值**：SLM 相位转换必须对 `exp(iφ)` 插值，不能直接对包裹相位插值
- **曝光**：相机峰值 ~200（8-bit），留足余量，避免饱和像素

## DMD 参数

实验室 DMD（数字微镜器件）用于相位加载，参数如下：

| 参数 | 值 |
|------|-----|
| 像元尺寸 | 13.7 µm |
| 原始分辨率 | 1024×768 |
| 有效区域（方形） | 768×768 |
| 超像素 | 4×4 mirrors |
| 超像素分辨率 | **192×192** |
| 超像素尺寸 | **54.8 µm** (4×13.7) |
| 物理面积 | 10.52×10.52 mm |

与 SLM（1920×1080, 6.4 µm）不同，DMD 通过 4×4 超像素编码相位，每个超像素对应一个相位值。计算 DOE 相位（2048×2048, dx_doe）通过裁剪中心 10.52 mm 区域 + 复振幅插值映射到 192×192。

转换工具：`dmd_phase/convert_to_dmd.py`。

## 基线

基线初始相位：`initial_phase_generation/artifacts/20260428-141942/phase0.mat`（429 mm, β_x=9.1, β_y=3.3）。

详见 `BASELINES.md`。
