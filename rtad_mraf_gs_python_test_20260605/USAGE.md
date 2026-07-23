# 使用说明

## 环境

```bash
conda activate slmrtad
```

Python 路径：`D:\software\anaconda\envs\slmrtad\python.exe`

GPU：RTX 5070 Ti（CuPy 自动启用），200 轮约 3 秒。

## 快速开始（完整跑一遍）

### Step 1：生成初始相位（已有基线可跳过）

基线相位已生成好，直接用：

```
initial_phase_generation/artifacts/20260428-141942/phase0.mat
```

如需重新生成（修改了光束直径等参数时）：

```matlab
cd E:\program\Point2P\initial_phase_generation
run_initial_phase_generation
```

### Step 2：仿真 WGS 精修

```bash
conda activate slmrtad
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
```

输出在 `artifacts/<时间戳>_rtad_mraf_gs/` 下，关键文件：

| 文件 | 内容 |
|------|------|
| `phase_refined.npy` | 精修后的相位（2048×2048） |
| `phase_refined.mat` | 同上，MATLAB 可读 |
| `reconstruction_refined.npy` | 焦面强度图 |
| `phase_refined.png` | 相位可视化 |
| `refined_reconstruction_intensity.png` | 焦面强度可视化 |
| `center_profiles_compare.png` | 初始 vs 精修 中心剖面 |
| `convergence_metrics.png` | 收敛曲线 |
| `metrics.csv` | 每 10 轮的指标记录 |
| `report.txt` | 完整文本报告 |

### Step 3：SLM 相位转换

将 2048×2048 计算相位转为 SLM 可用格式（1920×1080, 8-bit PNG）。

**注意**：`save_slm_phase.py` 目前在 `lab_test/lab_test_f200mm/` 中，使用时要传入正确的 f=429mm 参数：

```bash
python ../lab_test/lab_test_f200mm/save_slm_phase.py \
    artifacts/<你的输出目录>/phase_refined.npy \
    --out slm_phase.png \
    --f 0.429 \
    --wavelength 532e-9 \
    --focal-dx 2.5
```

f=429mm 时，计算像素尺寸 dx_doe ≈ 44.6 μm，SLM 物理区域（12.288×6.912 mm）对应约 276×155 个计算像素。`save_slm_phase.py` 会自动裁剪中心区域、做复振幅插值、输出 8-bit 灰度 PNG。

### Step 4：诊断已有结果

```bash
python run_diagnostics_case.py artifacts/<输出目录>
```

## 冒烟测试（验证环境）

不需要 phase 文件，快速验证程序能跑通：

```bash
python run_rtad_mraf_gs_case.py --iters 20 --smoke-shape 256
```

## 关键参数说明

### 方法选择

| `--method` | 含义 |
|------------|------|
| `wgs` | **推荐**。直接 WGS，200 轮 |
| `mraf` | 仅 MRAF |
| `mraf_then_wgs` | MRAF 预热 + WGS，不推荐（测试证明不带来额外收益） |

### WGS 参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--wgs-strategy` | `flat_local` | 只用这个。`xy_then_x` 会过早冻结 Y，损害均匀性 |
| `--iters` | 200 | 100 轮后改进很小，200 是安全值 |
| `--wgs-feedback-exponent` | 0.8 | 逐轮修正力度，0.8 收敛快 |
| `--wgs-weight-min` | 0.5 | 权重下限 |
| `--wgs-weight-max` | 2.0 | 权重上限，≥2.5 无进一步收益 |

### 背景处理

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--bg-factor` | 0.9 | **不要用 0.05**。近 1 表示不强压背景，平台更平 |

### 相位加载

| 参数 | 说明 |
|------|------|
| `--phase-mat` | phase0.mat 路径 |
| `--phase-var` | mat 文件中的变量名，通常是 `phase0_wrapped_rad` |
| `--swap-phase-xy` | MATLAB 生成的 phase0 默认会交换 x/y（默认 True，不要改） |

## 如何读懂结果

看 `center_profiles_compare.png`：
- X 方向中心剖面应平顶宽度 ~330 μm
- Y 方向中心剖面应平顶宽度 ~120 μm
- 平台区越平越好，RMS < 2% 可接受

看 `convergence_metrics.png`：
- RMS 随迭代下降则正常
- WGS 的 RMS 在首次权重更新后会跳升，之后缓慢下降——这是预期行为

终端输出的最终指标：
```
Final metrics:
  rms_nonuniformity_percent = 1.86%
  size50_x/y = 330.2 / 123.6 um
  efficiency_e2_percent = 92.5%
```

## 修改物理参数

编辑 `config_default.py`，或通过命令行覆盖。常用的：

```bash
# 修改焦距（正式管线不要改，f=429mm 是固定的）
# 不支持命令行覆盖焦距，需改 config_default.py

# 圆形光束：同时覆盖 X/Y 的 1/e^2 强度直径（单位 mm）
python run_rtad_mraf_gs_case.py --beam-diameter 6.5

# 椭圆光束：分别覆盖 X/Y 的 1/e^2 强度直径（单位 mm）
python run_rtad_mraf_gs_case.py --beam-diameter-x 6.5 --beam-diameter-y 6.3
```

`--beam-diameter` 保留为圆形光束兼容接口；指定椭圆光束时应同时传入
`--beam-diameter-x` 和 `--beam-diameter-y`，并加载按相同 X/Y 尺寸生成的 phase0。

## 球差与离焦补偿相位

2026-07-23 起，使用 `export_zernike_compensated_slm.py` 在 V2 上叠加实验
补偿。当前首个基准为 `Z40=+0.10625`、`Z20=+0.25000` RMS waves，安装补偿
保持 `X+5/Y+5`，闪耀光栅保持焦面 `X+200/Y+200 µm`。

完整命令、处理顺序和输出说明见 `README_ZERNIKE_COMPENSATION.md`。
