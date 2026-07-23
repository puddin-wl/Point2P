# 保留图片索引

## KEEP_measured_pupil_12mm_z1m_comparison.png

- 状态：明确保留，后续讨论截断影响时使用。
- 内容：实测 `G-光斑-1.bgData` 振幅叠加固定 V2 相位，自由传播 1 m 后，对比无入瞳截断与 12 mm 圆形入瞳截断的场镜后焦面。
- 条件：不加闪耀光栅，不加球差、离焦或其他 Zernike 项。
- 关键结果：12 mm 入瞳截去约 0.3085% 功率，中心比由 0.9849 变为 0.9903，没有产生中心空洞。
- 原始生成结果：
  `../results/12_measured_pupil_12mm_1m/measured_pupil_12mm_z1m_comparison.png`

## KEEP_measured_input_forward_diagnostics.png

- 状态：明确保留，后续汇报/讨论需要使用。
- 内容：
  - 左上：`G-光斑-1.bgData` 中框选后的实测 Gaussian 入射光；
  - 右上：重采样到 V2 DOE 网格后的实测强度；
  - 左下：2026-07-16 实验平顶光 16 帧平均；
  - 右下：实测入射振幅、平面波前和固定 V2 相位的正向传播结果。
- 关键用途：直观看出实测 Gaussian 振幅本身不会产生实验中的明显中心空洞。
- 原始生成结果：
  `../results/03_measured_input_forward/measured_input_forward_diagnostics.png`
