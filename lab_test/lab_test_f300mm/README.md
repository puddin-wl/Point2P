# f=300 mm 短焦距验证

该目录验证 532 nm、f=300 mm、5 mm 入射光和 330×120 μm 目标下的 WGS 行为。它位于 f=100/200 mm 几何光学区验证和 f=429 mm 正式场景之间。

## 运行

```matlab
cd('E:\program\Point2P\lab_test\lab_test_f300mm\matlab')
run_initial_phase_generation
```

```powershell
python run_rtad_mraf_gs_case.py `
  --phase-mat matlab/artifacts/<时间戳>/phase0.mat `
  --phase-var phase0_wrapped_rad `
  --method wgs --wgs-strategy flat_local --iters 200 `
  --wgs-feedback-exponent 0.8 `
  --wgs-weight-min 0.5 --wgs-weight-max 2.0 `
  --bg-factor 0.9
```

CPU 冒烟测试：

```powershell
python run_rtad_mraf_gs_case.py --smoke-shape 256 --iters 2 --no-cupy
```

## 历史结论

f=300 mm 时 βx≈13、βy≈5，Y 方向进入明显衍射区。历史结果 RMS 非均匀性约 1.22%，size50 约 329.3×118.1 μm，e⁻² 效率约 94.9%。

运行输出位于 `artifacts/`，默认不进入 Git。
