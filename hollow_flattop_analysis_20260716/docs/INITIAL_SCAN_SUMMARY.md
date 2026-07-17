# V2 平顶光中心空洞第一轮正向传播扫描

- V2 相位：`E:\program\Point2P\rtad_mraf_gs_python_test_20260605\artifacts\run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm\phase_refined.npy`
- 正向传播：复用项目 `src.propagation.forward_fft`（居中、正交 FFT）。
- 实验目标：`20260716-1.bmData` 的 16 帧平均。
- 扫描耗时：58.5 s。

## 基准复现

- 与保存的 `reconstruction_refined.npy` 相对 L2 误差：`9.99543e-08`。
- 基准中心窗口/核心均值：`0.9989`。

## 实验目标

- 物理核心框：`339.48 × 118.08 um`。
- 中心窗口：`36.90 × 36.90 um`。
- 中心窗口/核心均值：`0.7736`。
- 中间三分之一区域/左右区域：`0.8981`。
- 核心 RMS：`25.01%`。

## 当前最佳候选

- 名称：`horizontal_input_dip_d0.600_s1.750mm`
- 参数：`{"amplitude_dip_family": "horizontal_input_dip", "amplitude_dip_depth": 0.6, "amplitude_dip_sigma_mm": 1.75}`
- 综合差异分数：`2.7071`（越低越接近实验）。
- 中心窗口/核心均值：`0.7819`。
- 中间三分之一区域/左右区域：`0.9973`。
- 核心 RMS：`22.26%`。

## 每类物理因素的最佳结果

| 因素 | 最佳候选 | 分数 | 中心比 | 中部/两侧 | RMS |
|---|---|---:|---:|---:|---:|
| horizontal_input_dip | `horizontal_input_dip_d0.600_s1.750mm` | 2.707 | 0.782 | 0.997 | 22.3% |
| HG02_mix | `HG02_mix_+0.200` | 3.020 | 0.793 | 0.999 | 18.2% |
| center_dip | `center_dip_d0.300_s1.750mm` | 3.278 | 0.819 | 0.848 | 12.1% |
| vertical_input_dip | `vertical_input_dip_d0.300_s1.750mm` | 4.004 | 0.835 | 0.785 | 15.7% |
| astigmatism | `astigmatism_-0.200waves` | 4.238 | 0.882 | 0.993 | 33.0% |
| HG20_mix | `HG20_mix_+0.100` | 4.293 | 0.895 | 0.855 | 12.5% |
| beam_offset_x | `beam_offset_x_+0.600mm` | 5.069 | 0.952 | 0.932 | 31.6% |
| beam_offset_y | `beam_offset_y_-0.600mm` | 5.789 | 0.953 | 1.000 | 31.3% |
| defocus | `defocus_+0.200waves` | 5.985 | 0.902 | 1.091 | 26.9% |
| centered_aperture | `centered_aperture_9.500mm` | 6.642 | 0.977 | 0.997 | 9.8% |
| baseline | `baseline` | 8.034 | 0.999 | 1.000 | 0.8% |
| spherical | `spherical_+0.000waves` | 8.034 | 0.999 | 1.000 | 0.8% |
| measured_ellipse | `measured_Gaussian_X6p40_Y6p31mm` | 8.649 | 1.032 | 1.016 | 2.3% |
| aperture_offset_x | `aperture7mm_offset_x_-0.500mm` | 13.594 | 1.249 | 1.088 | 32.2% |
| aperture_offset_y | `aperture7mm_offset_y_-1.500mm` | 14.840 | 1.261 | 1.094 | 51.1% |

## 说明

- 这一轮固定 V2 相位，只改变入射振幅、孔径或附加波前。
- 闪耀光栅只平移焦斑，不影响中心形态，因此没有加入。
- 安装补偿等价于相位和入射光的相对横向偏移，本轮通过 beam offset 扫描。
- 第一轮结果用于锁定能够产生同类中心暗带的因素，不能单凭拟合分数认定真实物理原因。
