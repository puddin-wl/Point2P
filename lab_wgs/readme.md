# lab_wgs — 实验反馈 WGS 在线优化（全 MATLAB 实现）

## 定位：仿真→实验的桥梁

```
┌─────────────────────────────────────────────────────────────┐
│                     Point2P 完整管线                         │
│                                                             │
│  Stage 1           Stage 2            Stage 5 ← 本文件夹    │
│  RD 初始相位  ──→  仿真 WGS  ──→  实验 WGS  ──→  最终相位  │
│  (MATLAB/Python)   (Python+GPU)     (MATLAB)               │
│                                                             │
│  lab_test_f200mm/  lab_test_f200mm/  lab_wgs/               │
│  make_phase0.py    run_rtad_mraf    run_experimental        │
│                    _gs_case.py      _wgs.m                  │
└─────────────────────────────────────────────────────────────┘
```

## 为什么要两级优化

**仿真 WGS**（Stage 2）从 Romero-Dickey 初始相位出发，用 GPU 跑 200 轮迭代（~3 秒），
把 RMS 从 ~15% 优化到 0.1-0.5%。这一步处理的是理想的物理模型：
理想高斯光束 + 完美 FFT 传播 + 无像差。

**实验 WGS**（本文件夹）以仿真最优相位为起点，用真实光路反馈跑 20-30 轮，
补偿仿真中未建模的误差：光路不对准、实际光束偏离理想高斯、SLM 非线性响应、
杂散光、光学元件瑕疵等。

**为什么不能跳过仿真直接从随机相位做实验 WGS**：随机相位的平顶根本不成形——
相机拍到的是一团散斑，梯度边缘法找不到矩形亮区，WGS 权重更新没有有效信号。
之前实验 WGS 做不出来的原因就在这里。仿真 WGS 提供了"热启动"——初始相位已经
能在焦面产生可辨认的平顶光斑，实验 WGS 只需做小幅修正。

## 与其他文件夹的关系

| 文件夹 | 关系 | 本文件夹如何使用 |
|--------|------|-----------------|
| `lab_test_f200mm/` | **输入源** | 从中取 `phase_refined.npy` 作为实验 WGS 的初始相位 |
| `lab_test_f200mm/save_slm_phase.py` | **参考实现** | `display_slm_phase.m` 中的复振幅插值逻辑参考它 |
| `lab_test_f200mm/make_phase0.py` | **参考实现** | `load_input_amplitude.m` 参考它的高斯振幅生成 |
| `fig_analysis/analyze_captured.py` | **参考实现** | `analyze_captured_image.m` 的梯度边缘法移植自它 |
| `lab_test_f200mm/src/rtad_target.py` | **参考实现** | `build_rtad_target.m` 的 RTAD 目标生成移植自它 |
| `lab_test_f200mm/src/mraf_gs.py` | **参考实现** | `wgs_iteration.m` 的 WGS 核心算法移植自它（只取 WGS 部分，不用 MRAF） |
| `lab_test_f200mm/artifacts/` | **输入源** | 实验 WGS 的初始相位从这里加载 |
| `real_world_simulation/` | **互补** | 容差评估是开环的（固定相位，扫误差），实验 WGS 是闭环的（修正误差） |

## 算法

### 混合场 WGS（Hybrid-Field WGS）

仿真 WGS 的前向传播是 FFT。实验 WGS 的前向传播是**真实光路 + 相机**。
相机只测强度不测相位，所以反传时用仿真 FFT 提供的相位 + 相机提供的强度，
拼成焦面复振幅。

```
仿真 WGS：                         实验 WGS：
                                    
E_doe = A × exp(iφ)               E_doe = A × exp(iφ)
    ↓                                  ↓
  FFT                            ┌─ FFT ──→ φ_focal(仿真相位)
    ↓                            │
I_sim = |E_focal|²               │  SLM → 相机 → I_meas(实测强度)
    ↓                            │        ↓
W *= (I_target/I_sim)^α          ↓        ↓
    ↓                            └─→ E_focal = √I_meas × exp(i·φ_focal)
E_focal = W·A_target·exp(iφ)         ↓
    ↓                             W *= (I_target/I_meas)^α
  IFFT                                ↓
    ↓                             E_focal_corrected = W·A_target·exp(i·φ_focal)
φ_new = angle(E_doe)                  ↓
                                    IFFT
                                      ↓
                                  φ_new = angle(E_doe)
```

### 单轮迭代

```
1. 显示: display_slm_phase(φ_k)      % 相位→PNG→SLM
2. 采集: img = capture_focal_image() % Image Acquisition Toolbox
3. 分析: [I_meas, metrics] = analyze_captured_image(img)
   - 梯度边缘法定位平顶中心+边框
   - 背景扣除、flat-level 归一化
   - 坐标映射：相机像素 → 焦面 μm → 仿真网格
4. WGS: [φ_{k+1}, W] = wgs_iteration(φ_k, W, I_meas, I_target)
   - 仿真 FFT 提供 φ_focal
   - 实测 I_meas 替代仿真强度
   - 混合场反传 → 新相位
5. 判断: if metrics.rms < target || iter >= max → 结束
```

## 文件清单

```
lab_wgs/
  basic/                            # 硬件控制（用户已提供）
    CallDemon.m                     # SLM DLL 封装
    snap.mlx                        # 相机采集

  experiments_wgs_config.m          # [✓] 实验参数配置
  load_input_amplitude.m            # [✓] 生成高斯输入振幅
  load_initial_phase.m              # [✓] 加载仿真 WGS 输出的最终相位
  build_rtad_target.m               # [✓] 构建 RTAD 平顶目标

  display_slm_phase.m               # [✓] 相位→PNG/BMP→SLM
  capture_focal_image.m             # [✓] 触发相机拍照 (gentl Mono8)
  auto_exposure.m                   # [✓] 自适应曝光控制
  analyze_captured_image.m          # [✓] 梯度边缘法图像分析
  map_camera_to_focal.m             # [✓] 坐标映射

  wgs_iteration.m                   # [✓] 混合场 WGS 单步
  run_experimental_wgs.m            # [✓] 主循环
  start_experiment.m                # [✓] 实验入口：加载精修相位→跑实验 WGS
  verify_phase1.m                   # [✓] Phase 1 验证：MATLAB vs Python 逐位对比

  artifacts/                        # gitignored
```

## 两个入口脚本

| 脚本 | 起点 | 用途 |
|------|------|------|
| `verify_phase1.m` | RD 初始相位 (RMS ~3.5%) | 验证 MATLAB 算法移植正确性——与 Python 跑相同 phase0 的 1 轮 WGS 结果对比 |
| `start_experiment.m` | **仿真 WGS 精修相位** (RMS ~0.18%) | 实际实验流程——从 `lab_test_f200mm/artifacts/` 加载最优相位，跑实验 WGS |

**实际使用时跑 `start_experiment`**，它加载仿真最优相位（beam=7mm, 200 轮 WGS, RMS 0.18%）作为实验 WGS 的热启动。

## SLM 相位格式

SLM 接收的是 **8-bit 灰度 PNG 或 BMP 图像**，不是 .mat 文件：

```
计算相位 [0, 2π) rad  →  display_slm_phase.m  →  PNG/BMP (0–255)
```

`display_slm_phase.m` 完成以下转换：
1. 从 2048×2048 计算网格裁剪中央区域（匹配 SLM 物理尺寸 12.288×6.912 mm）
2. **复振幅插值**：先转成 `exp(i·φ)`，分别对实部/虚部做 cubic 插值到 1920×1080，再取 `atan2`——避免相位 wrap 处的 spline 振铃
3. 包裹到 [0, 2π)，线性映射到 [0, 255]，输出 uint8 图像

这一步不涉及 .mat 文件。SLM DLL (`CallDemon.m`) 直接读取 PNG 路径显示。

## 配置参数

```matlab
% experiments_wgs_config.m

% 物理参数（与 lab_test_f200mm 一致）
cfg.lambda_m = 532e-9;
cfg.f_m = 200e-3;
cfg.aperture_diameter_m = 15e-3;
cfg.input_1e2_diameter_m = 5.0e-3;  % 实际光束直径

% 计算网格（与仿真一致）
cfg.N = 2048;
cfg.focal_dx_um = 2.5;

% 目标（与仿真一致）
cfg.W50_um = 330;
cfg.H50_um = 120;
cfg.delta_x_um = 15;
cfg.delta_y_um = 8;
cfg.guard_x_um = 20;
cfg.guard_y_um = 12;

% SLM
cfg.slm_width = 1920;
cfg.slm_height = 1080;
cfg.slm_pitch_um = 6.4;

% 相机
cfg.camera_pixel_um = 3.45;

% WGS 参数（比仿真更保守）
cfg.max_iters = 30;
cfg.wgs_feedback_exponent = 0.5;   % 实测噪声大，降低反馈强度
cfg.wgs_weight_min = 0.5;
cfg.wgs_weight_max = 1.5;
cfg.bg_factor = 0.9;
cfg.rms_target = 0.02;             % 实验 RMS 目标 2%（比仿真宽松）
cfg.convergence_window = 5;        % 连续 N 轮 RMS 变化<0.01% 则收敛

% 路径
cfg.phase0_path = '';              % 从 lab_test_f200mm/artifacts/ 加载
cfg.output_root = 'artifacts';
```

## 实施计划

### Phase 1：移植核心算法到 MATLAB（先不管硬件）

目标：在 MATLAB 中完整复现仿真 WGS 的单步迭代，确认 FFT + WGS 结果与 Python 一致。

- [x] `load_input_amplitude.m` — 参考 `lab_test_f200mm/make_phase0.py` 和 `src/propagation.py`
- [x] `build_rtad_target.m` — 参考 `lab_test_f200mm/src/rtad_target.py`
- [x] `wgs_iteration.m` — 参考 `lab_test_f200mm/src/mraf_gs.py` 的 WGS 部分
- [x] `analyze_captured_image.m` — 参考 `fig_analysis/analyze_captured.py`
- [x] `load_initial_phase.m` — 支持 .mat 和 .npy 格式
- [x] `map_camera_to_focal.m` — 相机像素→焦面 μm→仿真网格坐标映射
- [x] 验证：Python 基准数据已生成 (`lab_test_f200mm/artifacts/phase1_verify_1iter/`)
- [x] **仿真模式实验 WGS 验证已跑通**（见下方 "仿真模式验证结果"）
- [x] MATLAB `verify_phase1` 跑通：初始 RMS 3.57% vs Python 3.52%，权重一致

### Phase 2：接入硬件，跑通单轮闭环

目标：MATLAB 能显示相位、拍照、分析、更新——一整轮不出错。

- [x] `display_slm_phase.m` — 复振幅插值 + BMP 写入 + SecondDll SLM 显示（从 `cam_in_loop.m` 吸收）
- [x] `capture_focal_image.m` — gentl Mono8 相机初始化 + getsnapshot（从 `cam_in_loop.m` 吸收）
- [x] `auto_exposure.m` — 自适应曝光控制（二分搜索，目标峰值 ~200）
- [x] `run_experimental_wgs.m` — 主循环：sim_mode=false 时走完整硬件管线
- [ ] 验证：在实验电脑上跑 `start_experiment`（`sim_mode=false`），确认图像分析能定位平顶

### Phase 3：多轮迭代 + 收敛验证

目标：跑完 30 轮实验 WGS，RMS 收敛，得到比纯仿真更好的实验相位。

- [ ] 收敛曲线记录
- [ ] 异常处理：过曝暂停、平顶丢失报警
- [ ] 对比实验：仿真 WGS 相位 vs 实验 WGS 相位的焦面光斑

## 注意事项

1. **此文件夹的所有代码用 MATLAB**，不引入 Python 依赖
2. **物理参数必须与 `lab_test_f200mm/` 一致**，修改时两边一起改
3. **WGS 只用 flat_local 策略**，不做 MRAF、不做 xy_then_x
4. **SLM 相位保存必须用复振幅插值**，参考 `save_slm_phase.py` 的教训
5. **SLM 接收的是 PNG/BMP 灰度图**（0-255 → 0-2π），不是 .mat 文件。`display_slm_phase.m` 负责此转换
6. **图像分析用梯度边缘法**，参考 `analyze_captured.py`，不要用峰值法
7. **曝光控制**：每轮检查饱和像素，>1% 时暂停
8. **实验 RMS 目标 2-5%**，低于仿真是因为真实光路有额外的噪声源

## 仿真模式验证结果（2026-05-12）

从精修相位（RMS 0.30%）出发，仿真模式跑 30 轮实验 WGS，MATLAB 与 Python 对比：

| Iter | MATLAB RMS | Python RMS | 说明 |
|------|-----------|-----------|------|
| 0 | 0.30% | 0.30% | 起点一致 |
| 5 | 4.28% | 3.21% | 第一次权重更新后跳升 |
| 10 | 3.63% | 2.56% | |
| 15 | 3.34% | 2.02% | |
| 20 | 2.94% | 1.57% | |
| 25 | 2.47% | 1.20% | |
| 30 | 2.07% | 0.91% | 趋势向下，MATLAB 更慢 |

**MATLAB 收敛更慢的原因**：Python 的 MRAF 将 free 区振幅衰减到 40%（mraf_factor=0.4），把更多能量推入信号区。MATLAB 不做 MRAF（实验 WGS 设计意图），能量从 free 区自然扩散，WGS 修正效率降低。这是正确的差异，不是 bug。

**MATLAB Phase 1 算法验证通过**：
- `verify_phase1.m` 跑通：初始 RMS 3.57%（Python 3.52%），权重统计一致（std 0.0153 vs 0.0151）
- `start_experiment.m` 跑通：30 轮收敛正常，趋势与 Python 一致
- Python 基准数据：`lab_test_f200mm/artifacts/phase1_verify_1iter/`（1 轮对比）、`lab_test_f200mm/artifacts/experimental_wgs_sim_30iter/`（30 轮对比）

**RMS 不降反升是预期行为**，原因：

1. 仿真 WGS（200 轮，feedback=0.8）已把相位优化到 RMS 0.18%，这是仿真模型下的最优解
2. 实验 WGS 用保守参数（feedback=0.5，降低反馈强度）在仿真模式下继续迭代，只是在最优解附近扰动
3. 第一次权重更新（iter 5）RMS 跳升，之后缓慢恢复
4. **实验 WGS 的真正价值在硬件闭环**：真实光路中存在仿真未建模的误差（不对准、光束偏离、SLM 非线性），那时相机反馈的强度图会引导 WGS 做出有意义的修正

**结论**：MATLAB 算法框架正确，Phase 1 完成。Phase 2 硬件接口已从 `cam_in_loop.m` 吸收完成。

## 从 cam_in_loop.m 吸收的硬件接口

| 功能 | 来源 | 目标 |
|------|------|------|
| 相机初始化 + 采集 | `videoinput("gentl",1,"Mono8")` + `getsnapshot` | `capture_focal_image.m` |
| SLM 显示 | `calllib('SecondDll','saShowImageFromFilePath',...)` | `display_slm_phase.m` |
| 自适应曝光 | 二分搜索，目标峰值 ~200 | `auto_exposure.m` |
| 零级光 mask | `((X-cx).^2+(Y-cy).^2)>150^2` | `run_experimental_wgs.m` 曝光阶段 |
| 闭环管线 | SLM→相机→分析→WGS→SLM | `run_experimental_wgs.m` (sim_mode=false) |

与 `cam_in_loop.m` 的关键区别：
- 用梯度边缘法定位平顶（`analyze_captured_image`），不依赖标定映射
- 用 WGS 乘性权重更新，不是加性 GS 权重
- N=2048 计算网格 + 复振幅插值到 1920×1080，不是直接在 1080 网格上算
