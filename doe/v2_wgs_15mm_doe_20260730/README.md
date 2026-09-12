# V2 WGS 15 mm 方形 DOE

## 相位来源

- 来源：`phase_refined.npy`，即 2026-07-15 最终 V2 的 WGS 精修相位。
- 原始数组：2048 × 2048。
- 原始计算网格：44.57578125 µm/像素，整幅约 91.2912 mm。
- 本输出取中心 15.000 mm × 15.000 mm 的物理区域，并对
  `exp(i*phase)` 做三次插值到 1024 × 1024，避免直接插值包裹相位产生跳变。
- 输出 DOE：15.000 mm × 15.000 mm，14.6484375 µm/像素。
- 未施加 X/Y 安装平移，未叠加闪耀光栅，未使用圆形掩膜。

## 主要文件

- `doe_v2_wgs_15mm_square_4mask.gds`：可由 KLayout 打开的四掩膜 GDS。
- `doe_phase_v2_wgs_15mm_1024.mat`：MATLAB 数据，变量 `doe_phase`、
  `doe_phase_rad` 和 `doe_level`。
- `doe_phase_v2_wgs_15mm_1024.npy`：Python 浮点相位，单位 rad。
- `doe_levels_v2_wgs_15mm_1024.npy`：量化台阶编号。
- `doe_phase_v2_wgs_15mm_16bit.png`：16-bit 包裹相位图。
- `doe_levels_v2_wgs_15mm.png`：量化台阶预览。
- `mask_bit_8.png`、`mask_bit_4.png`、`mask_bit_2.png`、
  `mask_bit_1.png`：四张二值掩膜预览。
- `doe_v2_wgs_15mm_preview.png`：总览图。
- `metadata.json`：尺寸、材料参数、量化统计及 SHA256。

## GDS 层定义

沿用原 `kalyout_doemake.m` 的 532 nm、折射率 1.458/1.00029 和
77.5 nm 基础刻蚀深度：

| GDS layer | 二进制权重 | 刻蚀深度增量 |
|---:|---:|---:|
| 1 | 8 | 620.0 nm |
| 2 | 4 | 310.0 nm |
| 3 | 2 | 155.0 nm |
| 4 | 1 | 77.5 nm |
| 100 | — | 15 mm 方形参考边界，不是刻蚀层 |

GDS 坐标单位为 1 µm，数据库精度为 1 nm。刻蚀深度通过 layer 定义，
不写入 GDS 的 `dbunit`。

## 量化说明

为与原程序逐式一致，本结果使用 `floor(height/77.5 nm)`。在给定折射率下，
2π 深度为 1162.308 nm，等于 14.9975 个基础深度，因此实际台阶编号为
0–14；四张二进制掩膜仍完整输出。若加工方要求严格的 16 个等相位级
（0–15 均匀出现），应先确认是否改为 `2π深度/16` 的基础刻蚀深度，不能
在未确认工艺参数时直接替换当前 77.5 nm 规则。
