# f=200 mm 实验室验证管线

该目录验证 532 nm、f=200 mm、不同入射光直径和 330×120 μm 目标下的完整流程：初始相位、WGS 精修、1920×1080 SLM 转换和焦面诊断。

> 此处 SLM 是 1920×1080、6.4 μm/像素的实验设备；不要与项目当前 1024×1024、17 μm/像素的目标 SLM 混用。

## 主要入口

- `make_phase0.py`：Python 初始相位；
- `matlab/`：MATLAB 初始相位；
- `run_rtad_mraf_gs_case.py`：单次 WGS；
- `run_all_beam_diameters.py`：多光束直径批处理；
- `save_slm_phase.py`：1920×1080 SLM 导出；
- `run_diagnostics_case.py`：结果诊断；
- `artifacts/`：运行输出，不进入 Git。

## 生成初始相位

Python 示例：

```powershell
python make_phase0.py --beam 5 --out make_phase0_output_5mm
python make_phase0.py --beam 6 --out make_phase0_output_6mm
python make_phase0.py --beam 7 --out make_phase0_output_7mm
```

MATLAB 示例：

```matlab
cd('E:\program\Point2P\lab_test\lab_test_f200mm\matlab')
run_all_beam_diameters
```

## WGS 精修

```powershell
python run_rtad_mraf_gs_case.py `
  --phase-mat make_phase0_output_7mm/phase0.mat `
  --phase-var phase0_wrapped_rad --beam-diameter 7 `
  --method wgs --wgs-strategy flat_local --iters 200 `
  --wgs-feedback-exponent 0.8 `
  --wgs-weight-min 0.5 --wgs-weight-max 1.5 `
  --bg-factor 0.9 --no-swap-phase-xy
```

## SLM 导出

```powershell
python save_slm_phase.py artifacts/<运行目录>/phase_refined.npy `
  --out artifacts/<运行目录>/slm_phase_f200mm.png
```

程序按 DOE 面实际物理尺寸裁切后，对 `exp(iφ)` 的实部/虚部分别插值，再取相角并量化成 8-bit。禁止直接对包裹相位做三次插值。

## 历史结果

| 入射光直径 | RMS 非均匀性 | size50_x | size50_y | e⁻² 效率 |
|---:|---:|---:|---:|---:|
| 4.5 mm | 0.25% | 328.9 μm | 118.1 μm | 96.8% |
| 5.0 mm | 0.12% | 328.2 μm | 117.8 μm | 96.5% |
| 5.5 mm | 0.13% | 326.5 μm | 118.1 μm | 95.7% |
| 6.0 mm | 0.09% | 326.9 μm | 118.1 μm | 95.7% |
| 7.0 mm | 0.18% | 327.1 μm | 118.3 μm | 97.6% |

这些数据只用于方法验证，不应替换 f=429 mm 的项目状态。
