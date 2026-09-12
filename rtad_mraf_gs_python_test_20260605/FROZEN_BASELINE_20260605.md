# 2026-06-05 冻结基线

本目录用于重现并扩展 2026-06-05 的 f=429 mm、6.5 mm 入射光基线，同时避免覆盖 `rtad_mraf_gs_python` 中的通用实现。

## 冻结输入

```text
artifacts/phase0_beam6p5mm/phase0.mat
artifacts/phase0_beam6p5mm/config_snapshot.mat
```

## 参数

- 波长：532 nm；
- 焦距：429 mm；
- 入射 Gaussian 1/e² 强度直径：6.5 mm；
- 通光孔径：15 mm；
- 目标：330×120 μm；
- 网格：2048×2048；
- 方法：直接 WGS，`flat_local`，200 轮；
- `mraf_factor=0.8`，`bg_factor=0.9`；
- 权重范围：0.5–1.5；
- `swap_phase_xy=False`。

## 重现

```powershell
python run_rtad_mraf_gs_case.py
```

所有新运行必须写入新的 `artifacts/<运行名>/`。不要修改冻结输入；需要测试新参数时使用命令行覆盖或复制到新实验目录。
