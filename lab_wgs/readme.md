# lab_wgs — 实验反馈 WGS 在线优化

## 概述

仿真 WGS（Stage 2）用 FFT 模拟焦面强度来更新相位。实验 WGS 用**真实的相机图像**
替代仿真：每次迭代把相位加载到 SLM 上，用相机拍下实际焦面光斑，分析图像得到
平顶均匀性，再反传更新相位。这是一个**闭环实验优化**。

```
┌─────────────┐     ┌──────────┐     ┌──────────┐     ┌──────────────┐
│ Python 计算  │ ──→ │ MATLAB   │ ──→ │   SLM    │ ──→ │  光学系统    │
│ WGS 权重更新 │     │ 加载相位 │     │ 显示相位 │     │  透镜+DOE    │
│ 反传+相位提取│     └──────────┘     └──────────┘     └──────┬───────┘
└──────┬──────┘                                              │
       │                                              ┌──────▼───────┐
       │                                              │   相机成像   │
       │                                              └──────┬───────┘
       │                                              ┌──────▼───────┐
       │                                              │ MATLAB snap  │
       └──────────────────────────────────────────────│ 保存 .mat    │
                                                      └──────────────┘
```

## 与仿真 WGS 的区别

| | 仿真 WGS（Stage 2） | 实验 WGS（Stage 5） |
|---|---|---|
| 前向传播 | FFT（计算） | 真实光路 + 相机 |
| 后向传播 | FFT（计算） | FFT（计算，复用现有代码） |
| 输入振幅 | 理想高斯模型 | 实际光束（含像差、不对准等） |
| 收敛判断 | 仿真 RMS | 相机图像 RMS |
| 迭代速度 | ~0.015s/iter（GPU） | 受限于 SLM 刷新+相机曝光（~0.5-2s/iter） |

关键优势：实验 WGS 自动补偿所有**仿真中未建模的误差**（光路不对准、光束像差、
SLM 非线性响应、光学元件瑕疵等）。

## 硬件接口

当前硬件控制代码位于 `basic/`：

| 文件 | 功能 |
|------|------|
| `basic/CallDemon.m` | SLM 控制：通过 `SecondDll.dll` 在 1920×1080 窗口显示图像 |
| `basic/snap.mlx` | 相机采集：触发相机拍摄，保存图像数据 |

SLM 关键 API（来自 `SecondDll.dll`）：
```matlab
% 显示单张图片
calllib('SecondDll', 'saShowImageFromFilePath', 'phase.png', 0, 0, 0, 1920, 1080, 1);

% 带超时控制的窗口（适合自动化循环）
calllib('SecondDll', 'Timeout_CreateWindow');
calllib('SecondDll', 'Timeout_ShowWindow', 1);
calllib('SecondDll', 'Timeout_ShowImageFromFilePath', 'phase.png', 0, 0, 0, 1920, 1080, 1, 2000);
calllib('SecondDll', 'Timeout_CloseWindow');
```

## 算法流程

每次迭代：

```
1. Python 保存当前相位为 1920×1080 PNG（复用 save_slm_phase.py）
2. 调用 MATLAB 脚本：
   a. 加载 PNG 到 SLM（CallDemon）
   b. 等待 SLM 稳定（~100ms）
   c. 触发相机拍照（snap）
   d. 保存图像为 captured.mat
3. Python 加载 captured.mat
4. 用梯度边缘法分析图像（复用 fig_analysis/analyze_captured.py）：
   - 找平顶区域中心
   - 提取中心剖面
   - 计算 size50, RMS 均匀性
5. 如果收敛（RMS < 阈值 或 迭代耗尽）→ 结束
6. 用实测强度替代仿真强度，计算 WGS 权重更新
7. 反传（FFT）→ 提取新相位 → 回到步骤 1
```

## 需要新建的文件

```
lab_wgs/
  basic/                         # 硬件控制（用户已提供）
    CallDemon.m                  # SLM DLL 封装
    snap.mlx                     # 相机采集
  
  run_experimental_wgs.py        # 主循环：协调 SLM/相机/WGS 迭代
  matlab_bridge.py               # Python ↔ MATLAB 桥接
  camera_analysis.py             # 相机图像→平顶强度分布 的转换
  experimental_wgs_config.py     # 实验 WGS 配置参数
  display_phase.m                # MATLAB：接收相位 PNG，加载到 SLM
  capture_image.m                # MATLAB：触发相机，保存 .mat
  
  artifacts/                     # 运行时输出（gitignored）
```

### `matlab_bridge.py` — 核心桥接逻辑

```python
# 方案 A：通过命令行调用 MATLAB
subprocess.run([
    'matlab', '-batch',
    f"cd {PROJECT_DIR}; display_phase('{png_path}'); capture_image('{out_mat}')"
])

# 方案 B：Python 写 PNG → MATLAB 监控文件夹 → 自动执行
# （适合 MATLAB 已打开 GUI 的情况）
```

具体方案取决于实验电脑上 MATLAB 的使用方式（命令行 vs GUI）。

### `camera_analysis.py` — 从相机图像提取平顶强度分布

关键步骤：
1. **平顶定位**：梯度边缘法找到矩形亮区（复用 `fig_analysis/analyze_captured.py`）
2. **背景扣除**：用图像四角估计背景，减去
3. **归一化**：以平顶区域中值为 1 做归一化
4. **坐标映射**：相机像素坐标 → 焦面物理坐标（μm）
5. **裁剪**：提取 target 对应的焦面区域，作为 WGS 的 "测量强度"
6. **遮挡处理**：相机图像中被 SLM 边框或其他遮挡物影响的区域需标记

## 关键挑战与对策

### 1. 坐标对齐
相机图像中的平顶位置每轮可能漂移。对策：
- 每轮自动检测平顶中心（梯度边缘法已支持）
- 以检测到的中心为基准裁剪 target 区域
- 如果漂移过大（>20μm），报警提示重新对准

### 2. 曝光控制
过曝图像无法用于 WGS 更新。对策：
- 每轮检查饱和像素比例
- 如果饱和 >1%，自动降低相机曝光或等待用户调整激光功率
- 初始几轮可先用低功率"对准模式"

### 3. 散斑噪声
相机图像有散斑噪声，直接用于 WGS 会导致权重震荡。对策：
- 每个相位拍 3-5 张取平均（如果 SLM 和相机支持）
- 在 WGS 权重更新前对图像做轻度高斯平滑（σ=1-2 px）

### 4. SLM 相位响应
SLM 的灰度-相位映射可能不是严格线性的。对策：
- 实验 WGS 不需要知道绝对相位值——相位更新是相对的（从当前相位出发）
- 非线性响应会被 WGS 迭代自动补偿

### 5. 收敛判据
仿真 RMS 和实验 RMS 数值不可直接比较。对策：
- 以实验 RMS 的相对变化为判据：连续 5 轮 RMS 变化 < 0.01% → 收敛
- 设置最大迭代次数（建议 30-50 轮，受限于实验时间）

## 实施步骤

### Step 1：验证硬件链路
```matlab
% MATLAB 中手动测试
loadlibrary('SecondDll.dll', 'SecondDll.h');
calllib('SecondDll', 'Timeout_CreateWindow');
calllib('SecondDll', 'Timeout_ShowWindow', 1);
calllib('SecondDll', 'Timeout_ShowImageFromFilePath', 'slm_phase_f200mm_d5mm.png', ...
    0, 0, 0, 1920, 1080, 1, 2000);
% 肉眼确认 SLM 上显示了正确的相位图
% 运行 snap 确认相机能拍到焦面光斑
```

### Step 2：建立 Python→MATLAB 通信
确认 Python 能调用 MATLAB 执行 SLM 显示和相机采集。

### Step 3：单次闭环测试
跑一轮：加载相位 → 拍照 → 分析 → 更新相位 → 对比新旧相位。

### Step 4：多轮迭代
跑完整的实验 WGS 循环，观察 RMS 收敛曲线。

### Step 5：对比仿真结果
将实验 WGS 的最终相位与仿真 WGS 的最终相位对比，分析差异来源。

## 与现有管线的衔接

实验 WGS 的输入是仿真 WGS 的输出（`phase_refined.npy`）。
实验 WGS 的输出是一个经过真实光路验证/修正的相位。

```bash
# Step 1: 仿真 WGS（已有）
cd lab_test_f200mm
python run_rtad_mraf_gs_case.py ...  # 得到 phase_refined.npy

# Step 2: 转为 SLM PNG
python save_slm_phase.py artifacts/.../phase_refined.npy --out slm_initial.png

# Step 3: 实验 WGS（新建）
cd lab_wgs
python run_experimental_wgs.py --phase slm_initial.png --beam 5.5 --iters 30
```
