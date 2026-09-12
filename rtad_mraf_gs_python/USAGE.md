# 通用 RTAD/WGS 管线使用说明

## 环境

从仓库根目录安装依赖：

```powershell
python -m pip install -r requirements.txt
python -m pip install -r requirements-gpu.txt
```

## 从头运行

```powershell
cd rtad_mraf_gs_python
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

## 诊断

```powershell
python run_diagnostics_case.py artifacts/<运行目录>
```

## 冒烟测试

```powershell
python run_rtad_mraf_gs_case.py --iters 2 --smoke-shape 256 --no-cupy
```

## 重要参数

| 参数 | 说明 |
|---|---|
| `--method wgs` | 直接运行 WGS |
| `--wgs-strategy flat_local` | 平顶核心区二维局部反馈 |
| `--wgs-feedback-exponent` | 权重反馈强度 |
| `--wgs-weight-min/max` | 权重截断范围 |
| `--bg-factor` | 远背景保留比例，当前常用 0.9 |
| `--beam-diameter` | 圆形 Gaussian 的 1/e² 强度直径，单位 mm |
| `--swap-phase-xy` | 仅用于需要交换轴方向的旧 MATLAB 输入 |
| `--no-swap-phase-xy` | Python 生成或方向已确认的输入 |

## SLM 转换

本目录默认目标 SLM 为 1024×1024、17 μm/像素。优先使用本目录 `convert_to_slm.py`，或使用测试目录中带 manifest 和回归校验的导出器。不要再调用 f=200 mm 目录下针对 1920×1080、6.4 μm SLM 的旧默认参数。

## 输出解释

- `phase_refined.npy`：纯 WGS 相位；
- `phase_refined.mat`：MATLAB 兼容副本；
- `reconstruction_refined.npy`：焦面重建；
- `config_used.json`：本次真实参数；
- `metrics.csv`：迭代趋势；
- `diagnostics_python/`：统一验收指标。
