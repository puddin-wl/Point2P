# Point2P 基线说明

本文件保留基线定义。最新使用状态见 [`docs/PROJECT_STATUS.md`](docs/PROJECT_STATUS.md)。

## 初始相位历史基线：2026-04-28

历史目录：

```text
initial_phase_generation/artifacts/20260428-141942
```

参数：

- βx=9.060975545，βy=3.294900198；
- 目标 330×120 μm；
- 波长 532 nm；
- 焦距 429 mm；
- 入射光 1/e² 强度直径 5 mm；
- 通光孔径 15 mm。

该目录是早期 5 mm Romero–Dickey 初始相位基线，不代表当前 6.5 mm 实验输入。目录被 Git 忽略，因此新环境应运行生成脚本，不应假定该本地路径存在。

## 仿真参考基线：2026-06-05

```text
rtad_mraf_gs_python/artifacts/20260605-144020_rtad_mraf_gs_truncI0135
```

- 入射光直径：6.5 mm；
- 目标：330×120 μm；
- RMS 非均匀性：1.40310042%；
- size50：329.5557×119.8249 μm；
- e⁻² 效率：96.2935%。

## 实验基线：2026-07-15 V2

```text
rtad_mraf_gs_python_test_20260605/artifacts/
run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm
```

纯 WGS 相位是 `phase_refined.npy`；当前确认的现场加载图是 `SLM_LOAD_20260715/phase_20260715_targetExpX330_v2_shiftX+5_Y+5.bmp`。

## 基线管理规则

- 不覆盖已命名基线；新结果使用日期和用途命名；
- 每个基线保存配置、命令、指标、核心文件和 SHA256；
- “候选”“扫描”“现场确认”必须显式区分；
- 大数组放外部归档或 Git LFS，仓库至少保留 manifest；
- 更新实验基线时同步更新 `docs/PROJECT_STATUS.md` 和日志索引。

