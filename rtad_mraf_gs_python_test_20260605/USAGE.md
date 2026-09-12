# 6.5 mm 冻结目录使用说明

## 重现冻结基线

```powershell
cd rtad_mraf_gs_python_test_20260605
python run_rtad_mraf_gs_case.py
```

默认读取仓库已跟踪的 `artifacts/phase0_beam6p5mm/phase0.mat`，并使用 6.5 mm 圆形 Gaussian、f=429 mm、目标 330×120 μm、直接 WGS 200 轮。

## 覆盖入射光

```powershell
# 圆形光束
python run_rtad_mraf_gs_case.py --beam-diameter 6.5

# 椭圆光束；必须加载与 X/Y 尺寸一致的 phase0
python run_rtad_mraf_gs_case.py `
  --beam-diameter-x 6.5 --beam-diameter-y 6.3 `
  --phase-mat artifacts/phase0_library/beamX6p5mm_Y6p3mm/phase0.mat `
  --no-swap-phase-xy
```

## 当前 V2

```text
artifacts/run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm
```

- `phase_refined.npy`：无安装偏移、无闪耀光栅的纯 WGS 相位；
- `SLM_LOAD_20260715/...X+5_Y+5.bmp`：当前现场确认加载图；
- Zernike 输出目录：历史候选；PBS 损伤根因确认后已停用，不替代当前 V2。

## 安装偏移

先仿真筛选：

```powershell
python verify_shift_sweep.py <phase_refined.npy> `
  --shifts-y -20,-15,-10,-5,0,5,10,15,20 `
  --shifts-x 0 --out <输出目录>
```

再生成 SLM 文件：

```powershell
python shift_sweep_slm.py <phase_refined.npy> `
  --out <输出目录> `
  --shifts-y -20,-15,-10,-5,0,5,10,15,20 `
  --shifts-x 0 --blaze-x 200 --blaze-y 200
```

## Zernike 补偿工具（历史保留）

后续已确认中心空洞来自 PBS 在极高功率下损坏，因此这些工具和参数不再用于
修正该问题，只保留复现和通用导出能力。

固定光轴导出使用 `export_zernike_compensated_slm.py`；整幅 V2+Zernike 相位共同移动的 Y 扫描使用 `export_zernike_compensated_y_sweep.py`。两者处理顺序不同，详见 `README_ZERNIKE_COMPENSATION.md`。

## 测试

从仓库根目录：

```powershell
python -m unittest discover -s rtad_mraf_gs_python_test_20260605 -p "test_*.py"
```
