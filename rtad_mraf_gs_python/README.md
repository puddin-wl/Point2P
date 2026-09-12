# RTAD/WGS 通用仿真管线

本目录是 Point2P 的通用 Python 仿真实现，默认焦距 429 mm。它负责读取或生成 Romero–Dickey 初始相位、建立矩形 RTAD 目标、执行 GS/MRAF/WGS 精修并输出诊断结果。

当前实验使用的 6.5 mm 冻结基线和 2026 年 7 月 V2 不在本目录继续修改，见 [`../rtad_mraf_gs_python_test_20260605/README.md`](../rtad_mraf_gs_python_test_20260605/README.md)。

## 默认参数

默认值以 [`config_default.py`](config_default.py) 为准：

| 参数 | 默认值 |
|---|---:|
| 波长 | 532 nm |
| 焦距 | 429 mm |
| 入射光 1/e² 强度直径 | 6.0 mm |
| 通光孔径 | 15 mm |
| 目标 W50×H50 | 330×120 μm |
| 计算网格 | 2048×2048 |
| 焦面采样 | 2.5 μm/像素 |
| SLM | 1024×1024，17 μm/像素 |
| 精修方法 | 直接 WGS，`flat_local`，200 轮 |

注意：这里的 6.0 mm 是通用目录默认值；2026-06-05 仿真参考基线和 7 月 V2 使用 6.5 mm。

## 文件说明

- `make_phase0.py`：生成 Romero–Dickey 初始相位；
- `run_rtad_mraf_gs_case.py`：单次精修入口；
- `run_pipeline.py`：组合管线；
- `run_diagnostics_case.py`：重新诊断已有结果；
- `convert_to_slm.py`：将计算相位转换到 SLM 网格；
- `shift_sweep_slm.py`：整幅相位安装偏移扫描；
- `add_blaze_grating.py`：叠加闪耀光栅；
- `src/`：传播、目标、WGS、指标、绘图和 MAT 输入输出；
- `artifacts/`：运行生成物，不进入 Git。

## 从头运行

```powershell
python make_phase0.py --beam 6.5 --out artifacts/phase0_beam6p5mm

python run_rtad_mraf_gs_case.py `
  --phase-mat artifacts/phase0_beam6p5mm/phase0.mat `
  --phase-var phase0_wrapped_rad `
  --beam-diameter 6.5 `
  --method wgs --wgs-strategy flat_local --iters 200 `
  --wgs-feedback-exponent 0.8 `
  --wgs-weight-min 0.5 --wgs-weight-max 1.5 `
  --bg-factor 0.9 --no-swap-phase-xy
```

由 Python 生成的初始相位使用 `--no-swap-phase-xy`。只有读取方向未经校正的旧 MATLAB 相位时才使用轴交换。

## 输出

每次运行在 `artifacts/<运行名>/` 下保存：

- `phase_refined.npy/.mat`：WGS 精修相位；
- `reconstruction_refined.npy`：重建强度；
- `config_used.json`：实际参数；
- `metrics.csv`、`report.txt`：迭代和结果记录；
- `diagnostics_python/`：统一诊断指标与图；
- 目标、掩膜、相位、收敛曲线和中心剖面图。

## 算法约定

- 目标先按强度构建，再取平方根得到目标振幅；
- WGS 只更新平顶核心区的局部权重；
- `bg_factor` 接近 1 表示不强压远背景；
- 复振幅插值必须先计算 `exp(iφ)`，分别插值实部/虚部后再取 `atan2`；
- 诊断指标是后处理，不参与优化约束。

详细命令见 [`USAGE.md`](USAGE.md)，当前项目状态见 [`../docs/PROJECT_STATUS.md`](../docs/PROJECT_STATUS.md)。
