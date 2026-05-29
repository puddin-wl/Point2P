# f=200mm, SLM 1024×1024 @ 17μm — 仿真验证管线

## 物理参数

| 参数 | 值 |
|------|-----|
| 波长 λ | 532 nm |
| 焦距 f | 200 mm |
| 计算网格 | 2048×2048 |
| 焦面采样 | 2.5 μm/pixel |
| 通光孔径 | 15 mm |
| 光束 1/e² 直径 | 5 / 6 / 7 mm |
| 目标尺寸 | 330×120 μm（50%强度） |
| SLM 分辨率 | 1024×1024 |
| SLM 像素 | 17 μm |
| SLM 物理尺寸 | 17.4×17.4 mm |
| β_x (beam=5mm) | 19.4 |
| β_y (beam=5mm) | 7.1 |

> β > 10 意味着进入了几何光学区域，稳相近似精度高，WGS 只需小幅度修正。

## 管线总览

```
Stage 1                    Stage 2                     Stage 3
make_phase0.py         run_rtad_mraf_gs_case.py    convert_to_slm.py
┌─────────────┐       ┌──────────────────┐       ┌──────────────┐
│ Romero-     │       │ WGS 仿真精修       │       │ 复振幅插值     │
│ Dickey      │ ────→ │ 2048×2048, 200轮  │ ────→ │ 2048→1024     │
│ 解析初始相位 │       │ GPU/CuPy 加速      │       │ @ 17μm        │
└─────────────┘       └──────────────────┘       └──────────────┘
                            ↓
                     run_all_beam_diameters.py  ← 一键跑完 5/6/7mm 三种光束
```

## 文件结构

```
lab_test_f200mm_slm_17um_1024/
├── config_default.py       # 默认配置（物理参数/网格/目标/精修/SLM）
├── make_phase0.py          # Stage 1: 生成 Romero-Dickey 初始相位
├── run_rtad_mraf_gs_case.py # Stage 2: WGS 仿真精修 (主入口, 50+ CLI 参数)
├── convert_to_slm.py       # Stage 3: 复振幅插值转换为 SLM 格式
├── run_all_beam_diameters.py # 批量: 自动串联 Stage 1→2→3, 跑 5/6/7mm
├── run_diagnostics_case.py  # 诊断: 分析已有结果 (profile crossing 等)
├── src/                    # 核心算法库
│   ├── backend.py           #   NumPy/CuPy 后端选择
│   ├── propagation.py       #   FFT 傅里叶透镜正反传播
│   ├── mraf_gs.py           #   GS/MRAF/WGS 迭代优化引擎
│   ├── rtad_target.py       #   升余弦矩形目标生成 + mask 区域划分
│   ├── metrics.py           #   迭代内指标计算 (RMS, size50, 效率)
│   ├── diagnostics.py       #   后处理诊断 (profile crossing, 导数旁瓣)
│   ├── io_mat.py            #   MATLAB .mat 读写 (v7 + v7.3 HDF5)
│   ├── plotting.py          #   所有 matplotlib 出图
│   └── utils.py             #   时间戳/目录/CSV/JSON 工具
└── artifacts/              # 输出目录 (每次运行一个时间戳子文件夹)
```

## 使用方法

### 单步运行

```bash
conda activate slmrtad
cd lab_test/lab_test_f200mm_slm_17um_1024

# Step 1: 生成初始相位
python make_phase0.py --beam 5 --out make_phase0_output/

# Step 2: WGS 精修
python run_rtad_mraf_gs_case.py \
    --phase-mat make_phase0_output/phase0.mat \
    --phase-var phase0_wrapped_rad \
    --beam-diameter 5 \
    --method wgs \
    --wgs-strategy flat_local \
    --iters 200 \
    --wgs-feedback-exponent 0.8 \
    --bg-factor 0.9

# Step 3: 转换为 SLM 格式
python convert_to_slm.py artifacts/<timestamp>/phase_refined.npy \
    --out artifacts/<timestamp>/slm_output/

# 可选: 诊断分析
python run_diagnostics_case.py artifacts/<timestamp>/
```

### 批量运行

```bash
python run_all_beam_diameters.py
```

自动生成 beam=5/6/7mm 三种初始相位 → WGS 精修 → SLM 转换，结果在 `artifacts/<timestamp>/slm_output/`。

## 关键设计决策

### 1. 复振幅插值（避免包裹相位振铃）

直接对 `[0, 2π)` 包裹相位做 cubic spline 插值时，2π→0 的跳变产生严重振铃。正确做法是对复振幅分量分别插值再取 angle：

```python
cfield = np.exp(1j * phase_crop)
real_z = zoom(cfield.real, zoom_ratio, order=3)
imag_z = zoom(cfield.imag, zoom_ratio, order=3)
phase_slm = np.arctan2(imag_z, real_z)
```

### 2. 截断 RTAD 目标

`constraint_mode='truncated_rtad'`：只有 `I_full ≥ release_level` 的像素才作为 signal 被 MRAF/WGS 约束。低强度的 full-template 尾部释放到 `mask_free`，不强制匹配数学模型。这样 GS/MRAF 有更多自由度，避免在边缘处浪费能量去拟合一个实际测不到的尾部。

区域划分：
- **mask_flat** — 平顶核心（140×49 μm，是 50% 以内收窄 δ 后的区域）
- **mask_edge_lock** — 边缘过渡区（signal 中去掉 flat 的部分）
- **mask_free** — guard_window 内、signal 外的自由区
- **mask_bg_far** — guard_window 外的远背景

### 3. WGS 权重更新

只在 `mask_flat` 上做 2D 局部反馈，权重逐轮更新：

```
w ← w × ( mean(|E_flat|) / |E_map| )^exponent
```

然后 clip 到 [wgs_weight_min, wgs_weight_max]（默认 0.5~1.5），再归一化。exponent 越大修正越激进，默认 0.8。

### 4. flat_local vs xy_then_x

- **flat_local**（当前默认）：直接在 mask_flat 上做 2D 权重反馈。X 和 Y 方向同时优化，不回过早冻结任一方向。经测试这是最优策略。
- **xy_then_x**：先 2D XY WGS，再切到 X-only WGS（Y 方向权重冻结）。早期测试发现 X-only 会使 Y 方向过早冻结，损害均匀性。

### 5. 不要强压背景

`bg_factor=0.9` 几乎保留背景原样，只极轻微衰减。早期用 `bg_factor=0.05` 强压背景反而拖累了平台均匀性。原因是背景区与信号区通过傅里叶变换耦合——强压背景会迫使能量再分布到平台区造成波纹。

## 与主管线的差异

| 项目 | 本管线 | 主管线 (rtad_mraf_gs_python) |
|------|--------|------------------------------|
| 焦距 f | 200 mm | 429 mm |
| 光束直径 | 5/6/7 mm | 5 mm |
| β_x / β_y | 19.4 / 7.1 | 9.1 / 3.3 |
| SLM 格式 | 1024×1024 @ 17μm | 1920×1080 @ 6.4μm |
| 用途 | 实验室方法验证 | 最终应用 |

β 越大越接近几何光学区，初始相位更精确。f=200mm 的 β 约为 f=429mm 的 2 倍，适合验证算法行为而不依赖 WGS 的大幅度修正。

## 当前结果

| 光束 | RMS | size50_x | size50_y | e⁻² 效率 |
|------|-----|----------|----------|----------|
| 5.0 mm | ~0.12% | ~328 | ~118 | ~96.5% |
| 6.0 mm | ~0.09% | ~327 | ~118 | ~95.7% |
| 7.0 mm | ~0.18% | ~327 | ~118 | ~97.6% |

GPU（RTX 5070 Ti）上 200 轮 WGS 约 3 秒。
