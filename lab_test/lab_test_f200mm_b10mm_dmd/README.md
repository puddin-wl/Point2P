# lab_test_f200mm_b10mm_dmd — f=200mm, beam=10mm, DMD 输出

物理参数：λ=532 nm, f=200 mm, 光束 1/e² 直径=10 mm, 通光孔径=15 mm, 目标 330×120 μm。

β_x ≈ 38.9, β_y ≈ 14.1（几何光学区），初始相位质量高，WGS 修正量小。

## 管线

```
Stage 1: make_phase0.py --beam 10
         → phase0.mat (2048×2048)

Stage 2: run_rtad_mraf_gs_case.py --beam-diameter 10 --method wgs ...
         → phase_refined.npy/.mat (2048×2048)

Stage 3: convert_to_dmd.py artifacts/<ts>/phase_refined.npy --out dmd_output/ --f 0.2
         → phase_dmd_192x192.mat  (DMD 相位, 变量 phase_dmd_rad)
         → phase_dmd_192x192.png  (8-bit 灰度预览)
```

## 使用

```bash
conda activate slmrtad
cd E:\program\Point2P\lab_test\lab_test_f200mm_b10mm_dmd
```

### Stage 1 — 生成初始相位

```bash
python make_phase0.py --beam 10 --out make_phase0_output/
```

输出 `phase0.mat`（2048×2048, β_x≈38.9, β_y≈14.1）。

### Stage 2 — WGS 精修

```bash
python run_rtad_mraf_gs_case.py \
    --phase-mat make_phase0_output/phase0.mat \
    --phase-var phase0_wrapped_rad \
    --beam-diameter 10 \
    --method wgs \
    --wgs-strategy flat_local \
    --iters 200 \
    --wgs-feedback-exponent 0.8 \
    --wgs-weight-min 0.5 \
    --wgs-weight-max 1.5 \
    --bg-factor 0.9 \
    --no-swap-phase-xy
```

输出在 `artifacts/<时间戳>/`，关键文件：`phase_refined.npy`、`phase_refined.mat`、`phase_refined.png`。

### Stage 3 — DMD 转换

```bash
python convert_to_dmd.py artifacts/<时间戳>/phase_refined.npy \
    --out dmd_output/ --f 0.2
```

输出：
- `dmd_output/phase_dmd_192x192.mat` — DMD 相位，变量 `phase_dmd_rad`，192×192，[0, 2π)
- `dmd_output/phase_dmd_192x192.png` — 8-bit 灰度预览

## DMD 参数

| 参数 | 值 |
|------|-----|
| 原始分辨率 | 1024×768 |
| 超像素 | 4×4 mirrors |
| 超像素分辨率 | 192×192 |
| 像元尺寸 | 13.7 µm |
| 超像素尺寸 | 54.8 µm |
| 物理面积 | 10.52×10.52 mm |

转换过程：从 2048×2048 计算网格裁剪中心 ~506×506 像素（对应 DMD 物理面积 10.52mm），复振幅插值到 192×192。

## 文件清单

```
lab_test_f200mm_b10mm_dmd/
├── config_default.py          # beam=10mm
├── make_phase0.py             # Stage 1: Python 初始相位生成
├── run_rtad_mraf_gs_case.py   # Stage 2: WGS 精修
├── run_diagnostics_case.py    # 诊断已有结果
├── convert_to_dmd.py          # Stage 3: 2048→192 DMD 转换
├── src/                       # 核心算法库
└── artifacts/                 # 运行输出
```
