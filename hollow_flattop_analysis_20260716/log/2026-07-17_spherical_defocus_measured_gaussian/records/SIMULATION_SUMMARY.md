# 仿真结果摘要

## 固定输入

- V2 artifact：
  `E:\program\Point2P\rtad_mraf_gs_python_test_20260605\artifacts\run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm`
- V2 相位：`phase_refined.npy`
- 正向传播：原项目 `src.propagation.forward_fft`
- 实测入射光：
  `E:\program\Point2P\real_test\20260716\G-光斑-1.bgData`
- Zernike pupil：15 mm clear aperture
- 系数单位：RMS waves

## 扫描规模

标准 Zernike 球差 + 离焦共完成 1047 个案例：

| 扫描 | 案例数 |
|---|---:|
| 大范围粗扫 | 297 |
| 近零细扫 | 561 |
| 关键区域局部细化 | 189 |

## 主案例

```text
Z40 = -0.10625 RMS waves
Z20 = -0.25000 RMS waves
```

联合波前：

```text
W/lambda =
  -0.10625 * sqrt(5) * (6*rho^4 - 6*rho^2 + 1)
  -0.25000 * sqrt(3) * (2*rho^2 - 1)
```

## 结果对比

| 指标 | 实验 | 理想 Gaussian + Z40/Z20 | 实测 Gaussian + Z40/Z20 | 实测 Gaussian + 平面波前 |
|---|---:|---:|---:|---:|
| 中心/核心均值 | 0.77359 | 0.77338 | 0.75287 | 0.97274 |
| 中部/两侧 | 0.89813 | 0.91947 | 0.88799 | 0.97692 |
| 核心 RMS | 25.01% | 18.08% | 22.66% | 5.86% |
| 50% footprint IoU | — | 0.99265 | 0.98249 | 0.98580 |
| 50% 宽×高 | 约 339×118 µm | 370×135 µm | 375×137.5 µm | — |

## 数值解释

- 实测 Gaussian 单独配合平面波前不会产生空洞。
- 加入相同 Z40/Z20 后，实测 Gaussian 案例的空洞比理想 Gaussian 案例略深。
- 实测 Gaussian 案例的中部/两侧和核心 RMS 更接近实验。
- 矩形 IoU 仍为 `0.98249`，说明实测振幅没有破坏主要矩形外轮廓。

## 符号说明

`Z40` 和 `Z20` 的正负号遵循当前 Python 正向传播的相位约定。将参数映射到
实验中的机械调焦方向、SLM 相位或补偿相位时，需要用已知离焦量做一次正负方向
标定，不能直接只按符号判断实际移动方向。

## 关键文件

- `../figures/03_MEASURED_GAUSSIAN_ZERNIKE_COMPARISON.png`
- `../data/measured_gaussian_zernike_summary.json`
- `../data/measured_gaussian_zernike_focal_intensity.npy`
- `../data/KEY_ZERNIKE_CASES.json`
- `../data/best_defocus_for_each_spherical.csv`
- `../data/local_refine_all_cases.csv`
