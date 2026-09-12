# Point2P 文档导航

文档按“当前状态、操作说明、技术原理、历史记录”四层组织。

## 当前状态

- [`PROJECT_STATUS.md`](PROJECT_STATUS.md)：当前仿真、实验、补偿和 DOE 状态；发生冲突时以这里为入口。
- [`DATA_POLICY.md`](DATA_POLICY.md)：哪些文件进入 Git，哪些保留在本机或外部归档。
- [`../baselines/README.md`](../baselines/README.md)：稳定结果的 SHA256 清单。

## 操作说明

- [`../rtad_mraf_gs_python/README.md`](../rtad_mraf_gs_python/README.md)：通用 f=429 mm 仿真管线。
- [`../rtad_mraf_gs_python_test_20260605/README.md`](../rtad_mraf_gs_python_test_20260605/README.md)：6.5 mm 冻结基线、V2 和补偿工具。
- [`../initial_phase_generation/README.md`](../initial_phase_generation/README.md)：MATLAB 初始相位生成。
- [`../real_world_simulation/USAGE.md`](../real_world_simulation/USAGE.md)：容差扫描命令。
- [`../result_diagnostics/README.md`](../result_diagnostics/README.md)：MATLAB 诊断接口。

## 技术原理

- [`../Point2P_Complete_Summary.md`](../Point2P_Complete_Summary.md)：算法、指标和项目架构说明。该文件保留历史参数，最新运行状态仍以 `PROJECT_STATUS.md` 为准。
- [`../lab_test/lab_test_f200mm_slm_17um_1024/CODE_GUIDE.md`](../lab_test/lab_test_f200mm_slm_17um_1024/CODE_GUIDE.md)：关键算法代码导读。

## 历史记录

- [`history/RTAD_PARAMETER_TRIALS_20260428_20260602.md`](history/RTAD_PARAMETER_TRIALS_20260428_20260602.md)：原两份重复 README 中的早期参数扫描全文。
- [`../rtad_mraf_gs_python_test_20260605/log/README.md`](../rtad_mraf_gs_python_test_20260605/log/README.md)：2026 年 7 月实验日志索引。
- [`../hollow_flattop_analysis_20260716/docs/ROOT_CAUSE_UPDATE_20260912.md`](../hollow_flattop_analysis_20260716/docs/ROOT_CAUSE_UPDATE_20260912.md)：中心空洞的 PBS 损伤根因更正。
- [`../hollow_flattop_analysis_20260716/README.md`](../hollow_flattop_analysis_20260716/README.md)：已归档的中心空洞仿真排查。
- [`../real_test/`](../real_test/)：实验图像与尺寸定义记录。
- [`../presentation/`](../presentation/)：阶段汇报材料。

命令、路径或指标若与当前代码不一致，应先更新 `PROJECT_STATUS.md`，再同步到对应操作文档，避免多份“最终版本”并存。
