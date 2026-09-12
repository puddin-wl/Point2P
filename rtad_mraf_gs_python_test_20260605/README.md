# 6.5 mm 冻结基线、V2 实验与补偿工具

本目录从 2026-06-05 的 6.5 mm 仿真基线演进而来，用于保存不会覆盖通用管线的实验版本。当前包含尺寸补偿 V2、SLM 安装偏移、输入光失配分析、中心空洞分析接口和 Zernike 补偿导出。

## 目录定位

- 通用算法开发：使用 `../rtad_mraf_gs_python/`；
- 重现 6.5 mm 冻结基线：使用本目录默认配置；
- 查找 7 月实验最终文件：先看 `log/README.md`；
- 生成球差/离焦补偿候选：看 `README_ZERNIKE_COMPENSATION.md`；
- 做安装偏移扫描：看 `README_SLM_SHIFT.md`。

## 冻结基线

本目录跟踪了最小可复现输入：

```text
artifacts/phase0_beam6p5mm/phase0.mat
artifacts/phase0_beam6p5mm/config_snapshot.mat
```

默认参数：532 nm、f=429 mm、6.5 mm 圆形 Gaussian、15 mm 孔径、目标 330×120 μm、2048² 网格、直接 WGS 200 轮。

```powershell
python run_rtad_mraf_gs_case.py
```

结果写入新的 `artifacts/<运行名>/`，不会覆盖冻结输入。

## 当前实验 V2

当前实验基线目录：

```text
artifacts/run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm
```

- 纯 WGS 相位：`phase_refined.npy`；
- 设计目标：W50=375.02 μm、H50=137 μm；
- 安装补偿：X+5、Y+5；
- 闪耀光栅：焦面 X+200 μm、Y+200 μm；
- 现场加载文件：`SLM_LOAD_20260715/phase_20260715_targetExpX330_v2_shiftX+5_Y+5.bmp`。

2026-07-23 生成的 Z40/Z20 文件是候选和扫描组，尚未取代该 V2 基线。

## 主要工具

- `run_rtad_mraf_gs_case.py`：支持圆形或椭圆 Gaussian 的 WGS 运行；
- `shift_sweep_slm.py`、`verify_shift_sweep.py`：安装偏移生成与仿真筛选；
- `analyze_beam_diameter_mismatch.py`：输入直径失配；
- `analyze_elliptical_beam_mismatch.py`：椭圆输入失配；
- `analyze_simulation_percent_energy.py`：仿真能量尺寸；
- `export_zernike_compensated_slm.py`：固定光轴的 Zernike 补偿导出；
- `export_zernike_compensated_y_sweep.py`：整幅 V2+Zernike 相位 Y 扫描。

## 测试

从仓库根目录运行：

```powershell
python -m unittest discover -s rtad_mraf_gs_python_test_20260605 -p "test_*.py"
```

## 文档

- [`FROZEN_BASELINE_20260605.md`](FROZEN_BASELINE_20260605.md)：冻结基线说明；
- [`USAGE.md`](USAGE.md)：参数和运行命令；
- [`README_SLM_SHIFT.md`](README_SLM_SHIFT.md)：安装偏移；
- [`README_ZERNIKE_COMPENSATION.md`](README_ZERNIKE_COMPENSATION.md)：球差/离焦补偿；
- [`log/README.md`](log/README.md)：7 月实验日志索引；
- [`../docs/PROJECT_STATUS.md`](../docs/PROJECT_STATUS.md)：跨目录当前状态。
