# 2026-06-08 Real Test 数据处理记录

## 数据位置

- 目录：`E:\program\Point2P\real_test\20260605-1`
- 扩束后光斑：`光斑.bgData`
- 最终矩形平顶光：`20260605-5.bmData`
- 运行环境：`D:\software\anaconda\envs\slmrtad\python.exe`

## 1. 确认 bmData/bgData 格式

- `.bmData` 和 `.bgData` 实际是 HDF5 文件。
- 像素数据路径：`/BG_DATA/1/DATA`
- 图像尺寸路径：
  - `/BG_DATA/1/RAWFRAME/WIDTH`
  - `/BG_DATA/1/RAWFRAME/HEIGHT`
- 像素标定路径：
  - `/BG_DATA/1/RAWFRAME/PIXELSCALEXUM`
  - `/BG_DATA/1/RAWFRAME/PIXELSCALEYUM`
- 额外需要考虑 binning：
  - `/BG_DATA/1/RAWFRAME/BINNINGX`
  - `/BG_DATA/1/RAWFRAME/BINNINGY`

## 2. 提取图片和原始数据

脚本/操作生成了：

- 输出目录：`real_test/20260605-1/extracted_images`
- 每个文件输出：
  - `.png`：预览图
  - `.npy`：原始一维 `int32 DATA`
  - `_metadata.json`：尺寸、曝光、增益、像素标定等元数据

重要结果：

- `20260605-5.bmData`
  - HDF5 图像尺寸：`964 x 724`
  - `PIXELSCALE = 3.69 um`
  - `BINNING = 2 x 2`
  - 有效像素尺寸：`7.38 um x 7.38 um`
  - 曝光：`457.123...`
  - 增益：`0.0`
- `光斑.bgData`
  - HDF5 图像尺寸：`5120 x 5120`
  - `PIXELSCALE = 4.5 um`
  - `BINNING = 1 x 1`
  - 有效像素尺寸：`4.5 um x 4.5 um`
  - 曝光：`44.59`
  - 增益：`0.0`

## 3. 尺寸分析脚本

新增两个独立脚本：

- `real_test/analyze_expanded_beam_size.py`
  - 分析扩束后光斑。
  - 输出：`real_test/20260605-1/beam_size_analysis`
- `real_test/analyze_rect_flattop_size.py`
  - 分析 `20260605-5.bmData` 矩形平顶光。
  - 输出：`real_test/20260605-1/flattop_size_analysis`

当前计算结果：

- 扩束后光斑：
  - D4σ：`6.373 mm x 6.578 mm`
  - 剖面 `1/e^2` 宽度：`5.920 mm x 6.236 mm`
  - FWHM：`3.423 mm x 3.453 mm`
- 平顶光：
  - size50：`328.6 um x 131.2 um`
  - size90：`277.1 um x 79.1 um`
  - size13.5：`375.2 um x 184.8 um`

## 4. 可视化

新增脚本：

- `real_test/make_real_test_visualizations.py`

输出：

- `real_test/20260605-1/visualizations/20260605-5_flattop_crop_uniformity.png`
  - 平顶光 ROI、相对强度、均匀性偏差图。
- `real_test/20260605-1/visualizations/expanded_beam_colormap_log.png`
  - 扩束光斑线性伪彩色和 log 强度图。
- `real_test/20260605-1/visualizations/20260605-5_flattop_raw_profiles.png`
  - X/Y raw profile。
- `real_test/20260605-1/visualizations/20260605-5_flattop_raw_profile_sizes.json`
  - raw profile 阈值宽度和目标尺寸反求结果。

## 5. Raw profile 上的阈值宽度

归一化定义：

- `100% = median(flat)`，即平顶核心区中位强度。
- 边界过冲会显示为 `>100%`。

当前 raw profile 外包络宽度：

| 阈值 | X size | Y size |
|---|---:|---:|
| 90% | `292.3 um` | `89.5 um` |
| 86.5% | `293.7 um` | `92.7 um` |
| 50% | `308.1 um` | `109.6 um` |
| 13.5% | `328.5 um` | `129.0 um` |

反求目标尺寸对应的强度位置：

| 目标尺寸 | 对应 profile 强度位置 |
|---|---:|
| X = `330 um` | `11.56% of median(flat)` |
| Y = `120 um` | `27.46% of median(flat)` |

这个结果说明：当前实验边缘在 X/Y 方向不满足“同一个强度阈值同时定义出 `330 x 120 um`”。

## 6. 尺寸定义备注

记录文件：

- `real_test/size_definition_notes.md`

其中记录了三类后续可能需要比较的定义：

- raw profile 强度阈值宽度
- Percent Energy 宽度
- raised-cosine edge fit 宽度

## 7. 平顶光偏斜确认

新增脚本：

- `real_test/make_flattop_skew_figure.py`

输出：

- `real_test/20260605-1/visualizations/20260605-5_skew_explanation.png`
- `real_test/20260605-1/visualizations/20260605-5_skew_metrics.json`

当前基于 50% 等值线的量化：

- 顶边拟合倾角：约 `8.8 deg`
- 底边拟合倾角：约 `4.5 deg`
- PCA 主轴角：约 `-4.2 deg`
- 说明：顶边和底边拟合角度不同，说明图像不是一个单纯旋转后的理想矩形；它还带有边缘畸变/剪切成分。

解释：

- 相机旋转只能对图像坐标施加刚性旋转。
- 振镜 + 场镜系统若存在非正入射、振镜轴未共轭、扫描平面与相机平面不平行、场镜非远心或系统存在仿射剪切/梯形畸变，则矩形会发生 shear/keystone 类型变形。
- 这类畸变不是单纯旋转，因此不管怎么转动相机，都只能改变整体角度，不能把剪切或梯形畸变消掉。

## 8. 明天可继续的方向

1. 尝试复现软件里的 `86.5%` 尺寸定义，尤其确认它是不是 Percent Energy。
2. 针对平顶光边沿做 raised-cosine fit，但需要处理边界过冲。
3. 用平顶偏斜图判断是否需要引入 affine 校正或四点透视/仿射标定。
4. 若后续要做实验纠偏，优先区分：
   - DOE/SLM 图案本身导致的形变；
   - 振镜/场镜映射导致的几何畸变；
   - 相机安装角度导致的全局旋转。
