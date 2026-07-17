# 二维椭圆输入凹陷细化扫描

- 案例数：294
- 耗时：93.7 s
- V2 相位保持不变；正向传播仍复用项目 `forward_fft`。

## 实验目标

- 中心窗口/核心均值：`0.7736`
- 中间三分之一/左右两侧：`0.8981`
- 核心 RMS：`25.01%`

## 最佳仿真

- 名称：`elliptical_dip_d0.600_sx4.000_sy1.500mm`
- 参数：`{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 4.0, "amplitude_dip_sigma_y_mm": 1.5}`
- 分数：`1.3438`
- 中心窗口/核心均值：`0.7609`
- 中间三分之一/左右两侧：`0.8997`
- 核心 RMS：`20.70%`

## 前十名

| 排名 | 参数 | 分数 | 中心比 | 中部/两侧 | RMS |
|---:|---|---:|---:|---:|---:|
| 1 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 4.0, "amplitude_dip_sigma_y_mm": 1.5}` | 1.344 | 0.761 | 0.900 | 20.7% |
| 2 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.7000000000000001, "amplitude_dip_sigma_x_mm": 4.0, "amplitude_dip_sigma_y_mm": 1.2}` | 1.538 | 0.808 | 0.902 | 23.8% |
| 3 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 6.0, "amplitude_dip_sigma_y_mm": 1.8}` | 1.693 | 0.753 | 0.942 | 21.3% |
| 4 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.7000000000000001, "amplitude_dip_sigma_x_mm": 6.0, "amplitude_dip_sigma_y_mm": 1.5}` | 1.702 | 0.745 | 0.938 | 26.8% |
| 5 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 6.0, "amplitude_dip_sigma_y_mm": 2.1}` | 1.741 | 0.743 | 0.934 | 20.6% |
| 6 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 6.0, "amplitude_dip_sigma_y_mm": 2.4}` | 1.751 | 0.746 | 0.929 | 19.4% |
| 7 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 6.0, "amplitude_dip_sigma_y_mm": 1.5}` | 1.803 | 0.785 | 0.950 | 21.0% |
| 8 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.6000000000000001, "amplitude_dip_sigma_x_mm": 4.0, "amplitude_dip_sigma_y_mm": 1.8}` | 1.806 | 0.724 | 0.883 | 21.3% |
| 9 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.8, "amplitude_dip_sigma_x_mm": 4.0, "amplitude_dip_sigma_y_mm": 1.2}` | 1.837 | 0.798 | 0.887 | 29.1% |
| 10 | `{"amplitude_dip_family": "elliptical_input_dip", "amplitude_dip_depth": 0.7000000000000001, "amplitude_dip_sigma_x_mm": 3.0, "amplitude_dip_sigma_y_mm": 1.2}` | 1.916 | 0.781 | 0.846 | 23.8% |
