# f=100 mm 短焦距验证

该目录验证 532 nm、f=100 mm、5 mm 入射光和 330×120 μm 目标下的 Romero–Dickey 初始相位与 WGS 行为。它是方法验证，不是当前实验生产管线。

## 运行

先在 MATLAB 中生成初始相位：

```matlab
cd('E:\program\Point2P\lab_test\lab_test_f100mm\matlab')
run_initial_phase_generation
```

再在本目录运行精修：

```powershell
python run_rtad_mraf_gs_case.py `
  --phase-mat matlab/artifacts/<时间戳>/phase0.mat `
  --phase-var phase0_wrapped_rad --method wgs `
  --wgs-strategy flat_local --iters 200
```

无 MATLAB 时可做 CPU 冒烟测试：

```powershell
python run_rtad_mraf_gs_case.py --smoke-shape 256 --iters 2 --no-cupy
```

## 结论

f=100 mm 时 βx≈38.9、βy≈14.1，稳相近似已接近几何光学区，初始相位本身较好。历史 WGS 结果 RMS 非均匀性约 0.27%，size50 约 330.8×116.9 μm，e⁻² 效率约 99.1%。

所有运行输出位于 `artifacts/`，默认不进入 Git。
