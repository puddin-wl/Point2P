# Point2P — DOE 平顶光斑设计管线

将高斯激光束整形成 330×120 μm 矩形平顶光斑。波长 532 nm，计算网格 2048×2048。

## 目录结构

```
Point2P/
├── rtad_mraf_gs_python/          # ★ 正式管线 (f=429mm, Python) — 最终应用目标
│   ├── run_rtad_mraf_gs_case.py  #   Stage 2: 仿真 WGS 精修
│   ├── run_diagnostics_case.py   #   诊断分析
│   └── src/                      #   核心算法库
├── initial_phase_generation/     # Stage 1: Romero-Dickey 初始相位生成 (MATLAB)
├── lab_test/                     # 实验室验证（非主管线）
│   ├── lab_test_f200mm/          #   验证: f=200mm 标准管线
│   ├── lab_test_f100mm/          #   验证: f=100mm 焦距变体
│   ├── lab_test_f300mm/          #   验证: f=300mm 焦距变体
│   └── lab_wgs/                  #   验证: 实验反馈 WGS 在线优化 (MATLAB)
├── real_world_simulation/        # 容差评估：固定相位、扫入射/光路误差
├── fig_analysis/                 # 捕获光斑图像分析（梯度边缘法）
├── result_diagnostics/           # 焦平面结果诊断 (MATLAB)
├── truncated_beam_phase/         # 截断光束相位模块
├── presentation/                 # 项目汇报 Beamer PDF
└── target/                       # 目标图形定义
```

详细内容见各子目录的 README。

## 正式管线（f=429mm，最终应用）

> 使用说明：`rtad_mraf_gs_python/USAGE.md`

物理参数：λ=532 nm, f=**429 mm**, 光束 1/e² 直径=5 mm, 通光孔径=15 mm, 目标 330×120 μm。

### Stage 1 — 生成初始相位

Romero-Dickey 解析法：在高斯光束与矩形平顶之间求解稳相能量守恒，得到可分离的 DOE 初始相位。

基线初始相位：`initial_phase_generation/artifacts/20260428-141942/phase0.mat`（β_x=9.1, β_y=3.3）。

### Stage 2 — 仿真 WGS 精修

**直接从初始相位做 WGS，不需要 MRAF 预热，不需要 X-only 阶段。** 经过多轮参数扫描验证，`method=wgs` + `strategy=flat_local` 的简洁配置给出了最好的平台均匀性。

```bash
conda activate slmrtad
cd rtad_mraf_gs_python

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
```

f=429mm 当前最佳结果：**RMS=1.86%, size50=~330×124 μm, e⁻² 效率=92.5%**。

已验证的关键参数选择：
- **method=wgs**：直接 WGS，不加 MRAF 预热段。测试证明 MRAF→WGS 的两段式不带来额外收益
- **strategy=flat_local**：只对平顶核心区做 2D 局部权重反馈。不做 xy_then_x（先 XY 再 X-only），后者让 Y 方向过早冻结，反而损害均匀性
- **bg-factor=0.9**：极弱的背景衰减，实质接近"不管背景"。早期用 0.05 强压背景反而拖累平台均匀性；0.9 是多次扫描确认的最优值
- **feedback-exponent=0.8**：较强的逐轮修正力度，收敛更快
- **weight-max=2.0**：max≥2.5 没有进一步收益

GPU（RTX 5070 Ti）上 200 轮约 3 秒。

### Stage 3 — SLM 相位转换

将 2048×2048 计算相位转为 SLM 原生分辨率（1920×1080, 6.4 μm 像素）。**关键教训**：必须用复振幅插值，不能直接对包裹相位做插值——2π→0 的跳变会导致 cubic spline 振铃。

```python
# 正确做法：对 exp(iφ) 的实部/虚部分别插值，再取 angle
c = np.exp(1j * phase_cropped)
phase_slm = np.arctan2(zoom(c.imag, ...), zoom(c.real, ...))
```

### Stage 4 — 光斑图像分析

用梯度边缘法分析实验拍到的平顶光斑照片。梯度边缘法利用边缘处梯度峰值定位平顶边界，然后用边界内中值强度做 flat-level——避免了峰值阈值法因中心hotspot导致 flat-level 偏高的问题。

```bash
python fig_analysis/analyze_captured.py <image.mat/.bmp/.png/.tif> --pixel-um 3.45
```

## 实验室验证

`lab_test/` 中的全部内容是在当前实验室环境下的方法验证，目的是确认管线各阶段可行后再应用到 f=429mm 正式管线。

### 仿真 WGS 验证（f=100/200/300mm）

在短焦距下验证 WGS 参数策略，利用更大的 β 值（更接近几何光学区）来确认算法行为。

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

β 越大，稳相近似越准确，初始相位越接近精确解。β > 10 进入几何光学区，WGS 修正幅度很小。β < 10 进入衍射区，WGS 明显改善均匀性。

### 实验反馈 WGS 验证（lab_wgs）

仿真 WGS 之后，用真实光路反馈跑实验 WGS，验证混合场方法能补偿仿真中未建模的误差（光路不对准、SLM 非线性、杂散光等）。纯 MATLAB 实现。

核心思路：前向传播走真实光路（SLM→相机→实测强度），反传走仿真 FFT（提供仿真相位 + 实测强度拼成焦面复振幅）。

**不能跳过仿真 WGS 直接从随机相位做实验 WGS**——散斑状态拍不到平顶，梯度边缘法找不到矩形亮区，WGS 没有有效反馈信号。仿真 WGS 是必须的热启动。

当前 `lab_wgs/` 配置为 f=200mm（验证环境）。应用到正式管线时需将配置改为 f=429mm。

详见 `lab_test/lab_wgs/readme.md`。

## 重要经验

- **不要强压背景**：`bg-factor` 接近 1.0 时平台更平，远处背景的能量代价在可接受范围内
- **MRAF 不需要**：直接 WGS 的均匀性优于 MRAF→WGS 两段式，没必要用 MRAF
- **X-only 不需要**：flat_local 的 2D 局部权重已足够；引入 X-only 会过早冻结 Y，损害均匀性
- **复振幅插值**：SLM 相位转换必须对 `exp(iφ)` 插值，不能直接对包裹相位插值
- **曝光**：相机峰值 ~200（8-bit），留足余量，避免饱和像素

## 基线

基线初始相位：`initial_phase_generation/artifacts/20260428-141942/phase0.mat`（429 mm, β_x=9.1, β_y=3.3）。

详见 `BASELINES.md`。
