# 文件清单

## figures

- `01_KEEP_measured_input_forward_diagnostics.png`
  - 实测 Gaussian、DOE 重采样、实验平顶光和实测振幅平面波前传播诊断图。
- `02_KEY_ZERNIKE_COMPARISON.png`
  - 实验、理想 V2、第三列主案例和内部形貌最接近案例。
- `03_MEASURED_GAUSSIAN_ZERNIKE_COMPARISON.png`
  - 实验、理想 Gaussian 主案例、实测 Gaussian 主案例和平面波前基准。
- `04_SLM_TO_FIELD_LENS_TRENDS.png`
  - 修正矩形ROI问题后，两组输入在0–1 m传播中的15 mm入瞳外功率、
    D99、D99.9和质心偏移。
- `05_SLM_TO_FIELD_LENS_SNAPSHOTS.png`
  - 原始整幅实测Gaussian与理想Gaussian在0、0.5、1.0 m传播面强度及
    15 mm圆形入瞳。
- `06_FIELD_LENS_PUPIL_FOURIER_TRENDS.png`
  - 15 mm圆孔截去功率及焦面变化随SLM到场镜距离的趋势。
- `07_FIELD_LENS_PUPIL_FOURIER_MONTAGE.png`
  - 两组输入在0、0.5、1.0 m处截断前后的场镜后焦面对比。
- `08_MEASURED_PUPIL_13_14MM_TRENDS.png`
  - 仅实测Gaussian下，14 mm和13 mm圆孔影响随场镜距离的趋势。
- `09_MEASURED_PUPIL_13_14MM_MONTAGE.png`
  - 无圆孔、14 mm和13 mm圆孔在0、0.5、1.0 m处的后焦面对比。

## data

说明：大型 `.npy` 数值阵列保留在本地，但不纳入Git；它们均可由归档脚本
重新生成。Git中保留CSV/JSON指标、关键PNG和脚本。

- `KEY_ZERNIKE_CASES.json`：理想 Gaussian 下的两个关键 Z40/Z20 案例。
- `best_defocus_for_each_spherical.csv`：局部细扫中每个球差对应的最佳矩形离焦。
- `local_refine_all_cases.csv`：局部细扫的全部 189 个案例。
- `measured_gaussian_zernike_summary.json`：实测 Gaussian 主案例完整指标。
- `measured_gaussian_zernike_focal_intensity.npy`：实测 Gaussian 主案例焦面原始强度。
- `slm_to_field_lens_scan.csv`：修正矩形ROI后两个输入分支、每0.1 m的入瞳统计。
- `slm_to_field_lens_scan.json`：自由空间传播模型和完整机器可读结果。
- `field_lens_pupil_fourier_scan.csv`：15 mm圆孔截断前后焦面指标。
- `field_lens_pupil_fourier_scan.json`：圆孔、场镜和传播模型的完整机器可读结果。
- `measured_G_z1m_focal_with_15mm_pupil.npy`：实测输入在1 m场镜距离、15 mm圆孔后的焦面强度。
- `ideal_6p5mm_G_z1m_focal_with_15mm_pupil.npy`：理想输入对应的焦面强度。
- `measured_pupil_13_14mm_scan.csv`：14 mm和13 mm入瞳逐距离焦面指标。
- `measured_pupil_13_14mm_scan.json`：对应模型与完整机器可读结果。
- `measured_G_z1m_focal_with_14mm_pupil.npy`：14 mm圆孔、1 m距离的焦面强度。
- `measured_G_z1m_focal_with_13mm_pupil.npy`：13 mm圆孔、1 m距离的焦面强度。

## scripts

- `refine_zernike_spherical_defocus.py`：球差 + 离焦局部细扫脚本快照。
- `simulate_measured_input_zernike_key_case.py`：实测 Gaussian 主案例脚本快照。
- `scan_slm_to_field_lens_free_space.py`：实际SLM网格角谱传播脚本快照。
- `simulate_field_lens_pupil_fourier.py`：15 mm圆孔后傅里叶面对比脚本快照。
- `simulate_measured_pupil_13_14mm.py`：实测Gaussian的14/13 mm圆孔测试脚本快照。

脚本快照用于记录当时实现；正式重跑优先使用分析根目录 `scripts` 下的同名脚本，
因为那里保留了正确的项目导入路径。

## records

- `DISCUSSION_AND_JUDGMENT.md`：本轮讨论、用户判断与当前共同结论。
- `SIMULATION_SUMMARY.md`：参数、模型和数值结果。
- `CURRENT_SCAN_RECORD_SNAPSHOT.md`：截至本次15 mm圆孔截断检查的全量记录快照。
- `SLM_TO_FIELD_LENS_SUMMARY.md`：15 mm场镜入瞳检查摘要。
- `FIELD_LENS_PUPIL_FOURIER_SUMMARY.md`：15 mm圆孔截断后焦面变化摘要。
- `MEASURED_PUPIL_13_14MM_SUMMARY.md`：14 mm与13 mm圆孔测试摘要。
- `FINAL_CONCLUSION_20260717.md`：本轮结束时的最终结论与后续实验方向。
