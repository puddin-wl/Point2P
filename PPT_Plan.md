# Point2P 项目报告 PPT 制作计划

基于 `Point2P_Complete_Summary.md` 和最佳结果 artifacts 生成 18 页项目汇报 PPT。
核心要求：**思路重要，代码不重要**。

## 叙事线

问题 → 理论四部曲 → 参数探索 → 结果呈现 → 设计洞察

## 18 页结构

| P# | 标题 | 核心信息 | 图片 | 节奏 |
|----|------|---------|------|------|
| P01 | 封面 | 从高斯到平顶：532nm DOE矩形光束整形 | refined_reconstruction 低透明度背景 | breathing |
| P02 | 问题定义 + 物理参数 | 设计DOE将532nm高斯→330×120μm矩形平顶；β_y=3.30在衍射区，需要迭代 | initial_reconstruction（小图） | breathing |
| P03 | 四篇论文总览 | Romero-Dickey→RTAD→MRAF→WGS 构成完整链路，从RMS>10%到1.86% | 协同关系图（SVG手绘） | dense |
| P04 | Romero-Dickey 初始相位 (1996) | 稳相法解析求解；β参数决定精度；β_y=3.30 → 需要迭代 | phase0.png | breathing |
| P05 | RTAD 升余弦目标 (2025) | 升余弦边缘替代陡峭边界→抑制Gibbs现象；截断RTAD：只约束I≥13.5% | target_profiles + target_intensity | breathing |
| P06 | 区域划分 + MRAF 自由区松弛 (2016) | 四区划分(signal/free/bg/flat)是均匀性的关键；自由区保留当前场，不强制归零 | masks.png（大图） | breathing |
| P07 | WGS 自适应权重反馈 (2022) | 暗区增大权重、亮区减小→像自动调平；flat_local：2D同时作用XY | wgs_weights_final.png | breathing |
| P08 | 理论协同总结 | Romero-Dickey给起点 + RTAD抑制泄露 + MRAF松弛约束 + WGS逐像素匀化 | phase0 vs phase_refined 并排 | dense |
| P09 | 11轮参数扫描 | 系统性扫参路线图；bg_factor→1.0是最大发现；MRAF预热/X-only被证明不需要 | 彩色表格 | dense |
| P10 | 最终工作基线 | method=wgs, flat_local, 200轮, feedback=0.8, weight=[0.5,2.0], bg=0.9 | wgs_weight_hist.png | dense |
| P11 | 蜕变：优化前后对比 | 从粗糙phase0到WGS精修后高均匀矩形平顶——RMS从>10%降至1.86% | initial + refined 并排（Hero） | breathing |
| P12 | X/Y剖面对比验证 | 优化后X/Y方向都呈现高度均匀平顶；Y方向改善最显著 | center_profiles_compare（大幅） | breathing |
| P13 | 收敛行为 | RMS首次权重更新后跳升是预期行为；200轮充分收敛 | convergence_metrics（大幅） | breathing |
| P14 | 最终DOE相位 | WGS优化是"微调"而非"重构"——相位保持清晰物理结构 | phase_refined（大幅居中） | breathing |
| P15 | 反直觉的三大发现 | ①不强压背景 ②MRAF预热没必要 ③X-only反而有害 | 纯排版卡片 | dense |
| P16 | 工程关键细节 | 复振幅插值；仿真WGS是实验热启动；RMS跳升是预期；GPU 3秒/200轮 | 纯排版 | dense |
| P17 | 方法论的普遍意义 | 四齿轮范式可迁移到任何光束整形问题；物理推理 > 盲目优化 | 齿轮图（SVG） | breathing |
| P18 | Thanks | 从一个物理起点出发，四篇论文的理论协同，最终RMS 1.86% | refined_reconstruction（居中） | breathing |

## 配色方案（深色主题，呼应532nm绿光激光）

| 角色 | HEX |
|------|-----|
| 背景 | #0D1117 |
| 次级背景 | #1A1D2E |
| 主色（激光绿） | #00E676 |
| 辅色（冷蓝） | #64B5F6 |
| 强调（琥珀） | #FFAB40 |
| 正文 | #E8E8EC |
| 辅助文字 | #8B8B96 |
| 分隔线 | #2A2D3A |

## 字体

"Microsoft YaHei", "PingFang SC", Arial, sans-serif（正文 18-20px，标题 36-40px）

## 图片素材来源

E:\program\Point2P\rtad_mraf_gs_python\artifacts\20260605-144020_rtad_mraf_gs_truncI0135\

| 文件 | 使用页面 |
|------|---------|
| initial_reconstruction_intensity.png | P02, P11 |
| refined_reconstruction_intensity.png | P01(bg), P11, P18 |
| center_profiles_compare.png | P12 |
| convergence_metrics.png | P13 |
| phase_refined.png | P08, P14 |
| phase0.png | P04, P08 |
| masks.png | P06 |
| target_profiles.png | P05 |
| target_intensity.png | P05 |
| wgs_weights_final.png | P07 |
| wgs_weight_hist.png | P10 |
