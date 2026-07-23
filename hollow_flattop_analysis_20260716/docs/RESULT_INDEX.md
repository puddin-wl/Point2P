# 结果索引

## 最终收口

`FINAL_CONCLUSION_20260717.md`

用途：记录2026-07-17排查结束后的最终判断、排除项、球差主案例和后续实验方向。

`ADDENDUM_20260721.md`

用途：记录真实16帧尺寸复算、12 mm入瞳、Zernike符号复核，以及转入SLM
正球差补偿的决定。

## 本轮归档日志

`../log/2026-07-17_spherical_defocus_measured_gaussian`

该目录集中保存球差 + 离焦讨论、用户判断、关键仿真图、JSON/CSV/NPY 数据和
脚本快照，后续回顾本轮结论时优先从该目录进入。

## 后续优先使用

### 已标记的实测 Gaussian 诊断图

`../reference_figures/KEEP_measured_input_forward_diagnostics.png`

用途：说明 `G-光斑-1.bgData` 的实测 Gaussian 振幅本身不会在固定 V2 相位下
产生实验空洞。

### 球差 + 离焦关键对比图

`../results/07_zernike_spherical_defocus_local_refine/KEY_ZERNIKE_COMPARISON.png`

用途：同时比较实验、理想 V2、外轮廓优先案例和内部形貌最接近案例。

### 球差 + 离焦关键参数

`../results/07_zernike_spherical_defocus_local_refine/KEY_ZERNIKE_CASES.json`

用途：保存两个关键案例的 Z40/Z20 系数、矩形指标和空洞指标，供下一轮仿真或
实验参数换算使用。

### 实测 Gaussian + 关键球差/离焦

`../results/08_measured_input_zernike_key_case/MEASURED_GAUSSIAN_ZERNIKE_COMPARISON.png`

用途：在相同 `Z40=-0.10625`、`Z20=-0.25000` 下，对比实验、理想 Gaussian、
实测 Gaussian，以及实测 Gaussian 平面波前基准。

机器可读指标：

`../results/08_measured_input_zernike_key_case/measured_gaussian_zernike_summary.json`

### SLM到场镜自由空间传播与15 mm入瞳

`../results/09_slm_to_field_lens_free_space/slm_to_field_lens_trends.png`

用途：显示原始整幅实测Gaussian与理想6.5 mm Gaussian在0–1 m传播中的
15 mm入瞳外功率、D99和D99.9；无闪耀、无Zernike像差。

传播面图：

`../results/09_slm_to_field_lens_free_space/slm_to_field_lens_snapshots.png`

完整数据：

- `../results/09_slm_to_field_lens_free_space/slm_to_field_lens_scan.csv`
- `../results/09_slm_to_field_lens_free_space/slm_to_field_lens_scan.json`
- `../results/09_slm_to_field_lens_free_space/SUMMARY_SLM_TO_FIELD_LENS.md`

### 15 mm圆孔截断后的傅里叶面

`../results/10_field_lens_pupil_fourier/field_lens_pupil_fourier_montage.png`

用途：并列比较0、0.5、1.0 m处实测/理想Gaussian在“不截断”和“乘15 mm
圆孔”后的429 mm场镜后焦面。

趋势图：

`../results/10_field_lens_pupil_fourier/field_lens_pupil_fourier_trends.png`

完整数据：

- `../results/10_field_lens_pupil_fourier/field_lens_pupil_fourier_scan.csv`
- `../results/10_field_lens_pupil_fourier/field_lens_pupil_fourier_scan.json`
- `../results/10_field_lens_pupil_fourier/SUMMARY_FIELD_LENS_PUPIL_FOURIER.md`

### 实测Gaussian的14 mm与13 mm入瞳

`../results/11_measured_pupil_13_14mm/measured_pupil_13_14mm_montage.png`

用途：只用原始整幅实测Gaussian，并列比较无圆孔、14 mm圆孔和13 mm圆孔
在0、0.5、1.0 m场镜距离下的后焦面。

趋势图：

`../results/11_measured_pupil_13_14mm/measured_pupil_13_14mm_trends.png`

完整数据：

- `../results/11_measured_pupil_13_14mm/measured_pupil_13_14mm_scan.csv`
- `../results/11_measured_pupil_13_14mm/measured_pupil_13_14mm_scan.json`
- `../results/11_measured_pupil_13_14mm/SUMMARY_MEASURED_PUPIL_13_14MM.md`

### 实测Gaussian的12 mm入瞳

`../results/12_measured_pupil_12mm_1m/measured_pupil_12mm_z1m_comparison.png`

用途：在1 m场镜距离比较无圆孔与12 mm圆孔。12 mm圆孔截去约`0.3085%`
功率，中心比由`0.9849`变为`0.9903`，仍未产生中心空洞。

### Zernike反号案例

`../results/13_measured_input_positive_zernike_key_case/EXPERIMENT_VS_MEASURED_GAUSSIAN_ZERNIKE.png`

用途：复核`Z40=+0.10625`、`Z20=+0.25000` RMS waves会把负号案例的中心
凹陷变为中心隆起，为后续正球差SLM补偿提供符号依据。

### 像散和彗差探索（暂停）

`../results/exploratory_14_astigmatism_coma/ASTIGMATISM_COMA_EXPLORATORY_MONTAGE.png`

用途：保存低阶非对称像差的独立探索。当前不进入第一轮SLM补偿。

## 球差 + 离焦扫描

- `05_zernike_spherical_defocus`：大范围粗扫，只用于确定数量级；大球差会破坏
  矩形，不应把该目录的 overview 当作最终结果。
- `06_zernike_spherical_defocus_fine`：`Z40=-1…+1` 的近零细扫。
- `07_zernike_spherical_defocus_local_refine`：关键区域局部细化，是当前主要
  结果目录。
- `08_measured_input_zernike_key_case`：把第三列关键 Z40/Z20 原样应用到实测
  Gaussian 振幅。
- `09_slm_to_field_lens_free_space`：实际SLM物理网格到场镜的0–1 m自由空间
  传播与15 mm入瞳检查；不加入Zernike像差。
- `10_field_lens_pupil_fourier`：在各场镜距离实际施加15 mm圆孔，再计算
  429 mm场镜后焦面；用于直接判断截断是否会制造中心空洞。
- `11_measured_pupil_13_14mm`：只用实测Gaussian测试14 mm和13 mm圆孔；
  两者均未产生实验量级的中心空洞。
- `12_measured_pupil_12mm_1m`：只用实测Gaussian测试1 m处12 mm圆孔；仍未
  产生中心空洞。
- `13_measured_input_positive_zernike_key_case`：主Z40/Z20案例同时反号的符号
  复核结果。
- `exploratory_14_astigmatism_coma`：像散/彗差独立探索；已暂停，不属于当前
  补偿主线。

每个细扫目录中：

- `all_cases.csv`：所有组合；
- `best_defocus_for_each_spherical.csv`：每个球差对应的最佳矩形离焦；
- `fine_scan_results.json`：完整机器可读结果；
- `fine_scan_trend.png`：球差—最佳离焦及内部指标趋势；
- `fine_scan_montage.png`：代表案例图。

## 旧方向与反例

- `01_initial_scan`：初始扰动探索；
- `02_elliptical_input_counterexample`：人为严重双瓣振幅反例，不作为实验归因；
- `03_measured_input_forward`：实测 Gaussian 振幅正向传播；
- `04_measured_input_low_order_wavefront`：早期非标准归一化低阶波前扫描。
