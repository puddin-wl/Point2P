# 平顶光尺寸定义记录

本文只记录后续可能要比较的几种尺寸定义；当前优先任务是：在 raw profile 上反求 `330 um x 120 um` 对应的相对强度位置。

## 1. Raw Profile 强度阈值宽度

- 定义：沿 X/Y 中心线取 profile，背景扣除后用 `median(flat)` 归一化。
- 宽度：在某个强度比例上，例如 `50%`、`86.5%`、`13.5%`，找左右交点或外包络跨度。
- 优点：直观，能直接显示边界过冲和平台起伏。
- 风险：边界过冲、局部热点、平台内部凹陷会影响交点；不同“从中心找交点”或“外包络找交点”的算法会给出不同结果。

## 2. Percent Energy 宽度

- 定义：找一个区域/宽度，使其包含总能量的某个比例，例如 `86.5%`。
- 背景：`86.5% = 1 - exp(-2)`，常用于高斯光束的 `1/e^2` 能量等价定义。
- 优点：可能与 Spiricon/BeamGage 等软件的 beam width 设置一致。
- 风险：用于矩形平顶光时物理意义不如高斯光束直接；需要明确是 X/Y 方向一维累计能量，还是二维区域累计能量。

## 3. Raised-Cosine Edge Fit

- 定义：用设计里的升余弦边沿模型拟合左右/上下边缘，再从拟合曲线读取 flat core、transition width、50%、13.5% 等尺寸。
- 优点：与 DOE 设计模型一致，对噪声更稳，能分离平顶核心和边沿过渡。
- 风险：实验边界过冲会拉偏拟合；需要对过冲点降权，或者只拟合中低强度下降沿。

## 当前采用

- 图：`real_test/20260605-1/visualizations/20260605-5_flattop_raw_profiles.png`
- JSON：`real_test/20260605-1/visualizations/20260605-5_flattop_raw_profile_sizes.json`
- 当前新增项：`target_size_equivalent_levels`
  - X 方向反求 `330 um` 对应的 `I / median(flat)`。
  - Y 方向反求 `120 um` 对应的 `I / median(flat)`。
