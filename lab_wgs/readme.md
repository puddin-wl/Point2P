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

## 文件清单（待实现）

```
lab_wgs/
  basic/                            # 硬件控制（用户已提供）
    CallDemon.m                     # SLM DLL 封装
    snap.mlx                        # 相机采集

  experiments_wgs_config.m          # [待实现] 实验参数配置
  load_input_amplitude.m            # [待实现] 生成高斯输入振幅
  load_initial_phase.m              # [待实现] 加载仿真 WGS 输出的最终相位
  build_rtad_target.m               # [待实现] 构建 RTAD 平顶目标

  display_slm_phase.m               # [待实现] 相位→PNG→SLM
  capture_focal_image.m             # [待实现] 触发相机拍照
  analyze_captured_image.m          # [待实现] 梯度边缘法图像分析
  map_camera_to_focal.m             # [待实现] 坐标映射

  wgs_iteration.m                   # [待实现] 混合场 WGS 单步
  run_experimental_wgs.m            # [待实现] 主循环

  artifacts/                        # gitignored
```

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

- [ ] `load_input_amplitude.m` — 参考 `lab_test_f200mm/make_phase0.py`
- [ ] `build_rtad_target.m` — 参考 `lab_test_f200mm/src/rtad_target.py`
- [ ] `wgs_iteration.m` — 参考 `lab_test_f200mm/src/mraf_gs.py` 的 WGS 部分
- [ ] `analyze_captured_image.m` — 参考 `fig_analysis/analyze_captured.py`
- [ ] 验证：用 Python 的 phase_refined 做一轮 WGS，MATLAB 和 Python 结果一致

### Phase 2：接入硬件，跑通单轮闭环

目标：MATLAB 能显示相位、拍照、分析、更新——一整轮不出错。

- [ ] `display_slm_phase.m` — 参考 `lab_test_f200mm/save_slm_phase.py` 的复振幅插值
- [ ] `capture_focal_image.m` — 基于 `basic/snap.mlx`
- [ ] `run_experimental_wgs.m` — 主循环框架
- [ ] 验证：手动跑一轮，确认图像分析能定位平顶

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
5. **图像分析用梯度边缘法**，参考 `analyze_captured.py`，不要用峰值法
6. **曝光控制**：每轮检查饱和像素，>1% 时暂停
7. **实验 RMS 目标 2-5%**，低于仿真是因为真实光路有额外的噪声源
