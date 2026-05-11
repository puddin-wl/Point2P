# lab_wgs — 实验反馈 WGS 在线优化（全 MATLAB 实现）

## 概述

仿真 WGS 用 FFT 模拟焦面强度。实验 WGS 把相位加载到 SLM 上，用相机拍下**真实的**
焦面光斑，分析实测强度来更新 WGS 权重，再通过 FFT 反传更新相位。全流程用 MATLAB。

```
┌──────────┐     ┌──────────┐     ┌───────────┐
│ MATLAB   │ ──→ │   SLM    │ ──→ │  光学系统  │
│ 相位→PNG │     │ 显示相位 │     │  透镜+DOE  │
└────┬─────┘     └──────────┘     └─────┬─────┘
     │                                  │
     │                            ┌─────▼─────┐
     │                            │  相机成像  │
     │                            └─────┬─────┘
     │                            ┌─────▼─────┐
     │                            │ Image Acq │
     │                            │ getsnapshot│
     └────────────────────────────│  .mat     │
                                  └──────┬────┘
                                         │
┌─────────────┐                           │
│ MATLAB WGS  │◄──────────────────────────┘
│ 权重更新    │
│ FFT 反传    │
│ 相位提取    │
└─────────────┘
```

## 为什么全用 MATLAB

- **SLM**：`SecondDll.dll` 的 `calllib` 接口是 MATLAB 原生支持的
- **相机**：使用 MATLAB Image Acquisition Toolbox（`videoinput` / `getsnapshot`），Python 没有等价库
- **FFT**：MATLAB 的 `fft2` / `ifft2` 与 Python CuPy 版本等价，虽然慢一些但足够用
- **一致性**：硬件控制、图像处理、WGS 计算都在同一个环境，避免语言间数据传递的麻烦

## 算法

### 与仿真 WGS 的对应关系

仿真 WGS（Python）：
```
E_doe = A_input × exp(iφ)  →  FFT  →  E_focal  →  |E_focal|² = I_sim
W_new = W_old × (I_target / I_sim)^α
E_doe_new = IFFT( W × A_target × exp(i×angle(E_focal)) )
φ_new = angle(E_doe_new)
```

实验 WGS（MATLAB）：
```
E_doe = A_input × exp(iφ)  →  FFT  →  E_focal  →  angle(E_focal)  ← 保留相位
                                                              ↓
                           SLM显示  →  相机拍照  →  I_meas  ← 实测强度
                                                              ↓
W_new = W_old × (I_target / I_meas)^α                        ↓
E_doe_new = IFFT( W × A_target × exp(i×angle(E_focal)) )  ←──┘
φ_new = angle(E_doe_new)
```

关键：**实测强度替代仿真强度，FFT 相位保留使用**。相机只能测强度不能测相位，
所以仿真 FFT 仍然需要跑——它提供焦面的相位信息用于反传。

### 单次迭代流程

```
1. E_doe = A_input × exp(i × φ_k)          % DOE 面复振幅
2. E_focal_sim = fftshift(fft2(ifftshift(E_doe)))  % FFT 得焦面复振幅（含相位）
3. φ_focal = angle(E_focal_sim)            % 保留仿真相位

4. 把 φ_k 转为 PNG，加载到 SLM              % 硬件显示
5. 相机拍照 → img_raw                       % 硬件采集
6. I_meas = analyze_image(img_raw)          % 梯度边缘法提取平顶区域强度
   - 找平顶中心
   - 背景扣除、归一化
   - 坐标映射（相机像素 → 焦面 μm）
   - 裁剪到 target 对应的焦面区域

7. W_k = W_{k-1} × (I_target ./ I_meas).^α  % WGS 权重更新
   W_k = clip(W_k, w_min, w_max)            % 限幅

8. E_focal_corrected = W_k .* A_target .* exp(i × φ_focal)  % 加权目标
9. E_doe_new = fftshift(ifft2(ifftshift(E_focal_corrected))) % 反传
10. φ_{k+1} = mod(angle(E_doe_new), 2π)     % 提取新相位，包裹到 [0, 2π)
```

步骤 6 是实验 WGS 最核心的新增功能。步骤 3 的仿真相位只参与步骤 8 的反传，
实测数据只贡献强度——两者结合的混合场才是物理上正确的反传输入。

## 文件结构

```
lab_wgs/
  basic/                            # 硬件控制（用户已提供）
    CallDemon.m                     # SLM DLL 封装
    snap.mlx                        # 相机采集（Image Acquisition）

  experiments_wgs_config.m          # 配置参数：波长、焦距、光束直径、target、WGS超参

  % ---- 初始化 ----
  load_input_amplitude.m            # 生成输入高斯振幅 A_input（复用 RD 初始相位逻辑）
  load_initial_phase.m              # 从 .mat 加载 phase0 或 phase_refined

  % ---- 硬件交互 ----
  display_slm_phase.m               # 相位数组 → PNG → CallDemon 显示到 SLM
  capture_focal_image.m             # 调用 Image Acquisition 拍照，返回 uint8 图像

  % ---- 图像分析（核心新增） ----
  analyze_captured_image.m          # 梯度边缘法定位平顶、提取实测强度 I_meas
  map_camera_to_focal.m             # 相机像素坐标 → 焦面物理坐标（μm）

  % ---- WGS 迭代 ----
  wgs_iteration.m                   # 一次迭代：仿真FFT + 权重更新 + 反传
  run_experimental_wgs.m            # 主循环：显示→拍照→分析→WGS→重复

  % ---- 工具 ----
  save_phase_png.m                  # 相位数组 → 1920×1080 8-bit PNG（等同 save_slm_phase.py）
  load_phase_mat.m                  # 读取 phase_refined.mat / phase0.mat

  artifacts/                        # 运行时输出（gitignored）
    <timestamp>_expwgs/
      phase_001.png ... phase_N.png # 每轮 SLM 相位
      captured_001.mat ...          # 每轮相机图像
      I_meas_001.mat ...            # 提取的实测强度
      metrics.csv                   # 每轮 RMS / size50
      phase_final.mat               # 最终相位
```

## 各模块详细设计

### `analyze_captured_image.m` — 图像→强度分布

输入：相机原始图像 (2048×2448 uint8)
输出：I_meas — 映射到焦面坐标的归一化强度分布

步骤：
1. **背景估计**：取图像四角 30×30 区域，中值 = 背景
2. **寻找亮区**：img > bg + 5*std → 连通域 → 取最大 → 亮度加权质心
3. **梯度边缘定位**：
   - 过质心取 x/y 中心剖面，平滑 (窗口=7)
   - 计算梯度，找中心两侧最强梯度峰 → 边缘位置 lx, rx, ly, ry
   - 平顶参考电平 = 边缘间剖面中值
4. **尺寸测量**：以平顶电平为参考，线性插值找 50%/90%/13.5% 交叉点
5. **坐标映射**：
   - 已知靶面尺寸 (330×120 μm)，size50 的像素值 → μm/px 标定
   - 将图像裁剪/重采样到仿真焦面网格 (N=2048, dx=2.5μm)
6. **强度归一化**：I_meas = (img - bg) / flat_level → 平顶区域均值≈1

```matlab
function [I_meas, metrics] = analyze_captured_image(img, target_params)
    % target_params: W50_um=330, H50_um=120, pixel_um=3.45
    % I_meas: 归一化焦面强度，尺寸匹配仿真网格
    % metrics: size50, RMS, hotspot_ratio 等
```

### `display_slm_phase.m` — 相位显示到 SLM

```matlab
function display_slm_phase(phase_rad, slm_width, slm_height)
    % phase_rad: 2048×2048, [0, 2π)
    % 1. 裁剪中心区域 → 重采样到 slm_width×slm_height
    % 2. 复振幅插值（避免包裹跳变）
    % 3. 量化到 8-bit → 保存为临时 PNG
    % 4. CallDemon 加载 PNG 到 SLM
```

### `wgs_iteration.m` — 核心 WGS 一步

```matlab
function [phase_new, weights_new, metrics] = wgs_iteration(phase, weights, ...
    A_input, I_target, A_target, I_meas, cfg)
    % phase: 当前 DOE 相位
    % I_meas: 相机实测的焦面强度（已配准到仿真网格）
    % 
    % 1. E_focal_sim = FFT(A_input .* exp(1i*phase))  → φ_focal
    % 2. weights_new = weights .* (I_target ./ I_meas).^α, clip
    % 3. E_corrected = weights_new .* A_target .* exp(1i*φ_focal)
    % 4. phase_new = angle(IFFT(E_corrected)), wrap to [0,2π)
```

### `run_experimental_wgs.m` — 主循环

```matlab
function run_experimental_wgs(phase0_path, cfg)
    % 初始化
    phase = load_initial_phase(phase0_path);
    A_input = load_input_amplitude(cfg);
    [I_target, A_target] = build_rtad_target(cfg);  % 复用现有 target 逻辑
    weights = ones(size(I_target));
    
    for iter = 1:cfg.max_iters
        % 显示
        display_slm_phase(phase, 1920, 1080);
        pause(0.2);  % SLM 稳定
        
        % 采集
        img = capture_focal_image();
        
        % 分析
        [I_meas, m] = analyze_captured_image(img, cfg);
        
        % 记录
        fprintf('[%d] RMS=%.2f%%  size50=%.0fx%.0f um\n', ...
            iter, m.rms_pct, m.size50_x_um, m.size50_y_um);
        
        % 收敛判断
        if m.rms_pct < cfg.rms_target, break; end
        
        % WGS 迭代
        [phase, weights] = wgs_iteration(phase, weights, ...
            A_input, I_target, A_target, I_meas, cfg);
    end
    
    save(fullfile(cfg.out_dir, 'phase_final.mat'), 'phase', 'weights');
end
```

## WGS 参数建议

沿用仿真 WGS 的经验：

| 参数 | 值 | 说明 |
|------|-----|------|
| wgs_feedback_exponent | 0.8 | 实测噪声比仿真大，可以稍微降一点（0.5-0.8） |
| wgs_weight_min / max | 0.5 / 1.5 | 沿用仿真设置 |
| bg_factor | 0.9 | 轻度背景衰减 |
| max_iters | 30-50 | 实验迭代慢，不可能跑 200 轮 |
| rms_target | 2% | 实验 RMS 通常会比仿真高 |

注意：实验 WGS 的 RMS 和仿真 WGS 的 RMS 不可直接对比——实验 RMS 包含了
相机噪声、散斑、杂散光等仿真中没有的因素。合理的实验 RMS 目标可能是 2-5%。

## 实施步骤

### Step 1：基础验证
- 在 MATLAB 中手动加载一张 SLM 相位 PNG，拍照，确认能看到平顶光斑
- 验证 `analyze_captured_image` 能正确定位和测量

### Step 2：单轮闭环
- 跑一轮完整的：显示相位 → 拍照 → 分析 → WGS 更新 → 得到新相位
- 对比新旧相位的差异

### Step 3：多轮迭代
- 跑 30 轮，观察 RMS 是否收敛
- 如果震荡，降低 `wgs_feedback_exponent`
- 如果不动，检查坐标映射是否正确

### Step 4：对比评估
- 将实验 WGS 最终相位与仿真 WGS 最终相位对比
- 将两者都加载到 SLM 上拍照，直接对比焦面光斑

## 关键注意事项

1. **曝光**：每轮检查饱和像素，饱和 >1% 时暂停并提示调整激光/曝光
2. **坐标对齐**：每轮自动检测平顶位置，如果漂移 >20μm 报警
3. **散斑**：如果条件允许，每个相位拍 3 张取平均；或加大分析时的平滑窗口
4. **SLM 灰度响应**：不需要知道绝对映射关系，WGS 权重更新是相对的
5. **FFT 复用**：反传用的 FFT 可以用 CPU（MATLAB 内置），不需要 GPU。2048×2048
   的 FFT 在 MATLAB 中约 0.05-0.1 秒，不是瓶颈（硬件时间远大于此）
