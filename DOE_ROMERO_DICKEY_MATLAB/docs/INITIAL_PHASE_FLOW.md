# Initial Phase Flow

一句话总结：当前 MATLAB 主流程中的 `phase0` 是 **解析式 separable Romero-Dickey Gaussian-to-rectangular-flat-top 初始相位**，由 x/y 两个 1D RD 相位相加得到；它不是当前代码里的迭代 point-to-point 映射结果。现在干净入口已分离到 `E:\program\Point2P\initial_phase_generation\generate_initial_phase.m`。

## 调用链

### 不跑 MRAF，只看初始相位

```text
run_initial_phase_only.m
  -> default_config.m
  -> initial_phase_generation/generate_initial_phase.m
       -> internal grid generation
       -> internal Gaussian input/aperture generation
       -> internal separable Romero-Dickey phase generation
       -> optional internal ideal FFT forward propagation
  -> artifacts/initial_phase_only/<timestamp>/
```

### 原始 RD baseline

```text
run_one.m
  -> default_config.m
  -> initial_phase_generation/generate_initial_phase.m
  -> fresnel_to_focal_fft.m
  -> full_lens_propagate.m
  -> save_artifacts.m
```

### MRAF 使用 phase0

```text
run_mraf_one.m
  -> default_config.m
  -> initial_phase_generation/generate_initial_phase.m
  -> run_mraf_refinement(input_amplitude, phase_rd_unwrapped, focal_x_m, focal_y_m, cfg)
       -> mraf_forward_fft(input_amplitude .* exp(1i * phase0))
       -> MRAF target/projection iterations
```

## 核心文件列表

- `initial_phase_generation/generate_initial_phase.m`：当前推荐的独立 phase0 生成入口；返回网格、输入振幅、孔径、phase0 和可选初始焦平面强度。
- `initial_phase_generation/README.md`：说明外部/下游脚本如何调用独立 phase0 模块。
- `config/default_config.m`：定义波长、焦距、DOE 网格、孔径、入射光斑、目标尺寸、FFT 采样，以及 RD beta 相关参数。
- `src/make_grid.m`：生成 DOE/SLM 面的物理坐标 `X_m/Y_m`。
- `src/gaussian_input_field.m`：生成 DOE 面 5 mm 1/e² 高斯入射振幅和 15 mm clear aperture mask。
- `src/build_separable_phase_2d.m`：初始相位 `phase0` 的直接构造函数；把 x/y 两个 1D RD 相位相加。
- `src/romero_dickey_phase_1d.m`：真正计算 RD 解析相位公式的核心函数。
- `src/wrap_phase_2pi.m`：把相位包裹到 `[0, 2π)`，用于显示/导出。
- `src/fresnel_to_focal_fft.m`：用 Fourier lens FFT 模型检查初始相位焦平面强度。
- `src/mraf/run_mraf_refinement.m`：不生成 phase0，只接收 `phase0` 并用它作为 MRAF 初值。
- `run_one.m`：原始 RD baseline 入口。
- `run_mraf_one.m`：MRAF 入口，同时在 MRAF 前生成 RD `phase0`。
- `export_phase_for_python.m`：把 RD `phase0` 导出到 `exports/phase_rd.mat`。
- `run_initial_phase_only.m`：新增的最小初始相位检查入口。

## 参数来源表

| 参数 | 变量名 | 单位 | 默认值 | 来源文件 | 作用 |
|---|---:|---:|---:|---|---|
| 波长 | `cfg.lambda_m` | m | `532e-9` | `config/default_config.m` | RD beta 和 Fourier/focal 坐标标定 |
| 焦距 | `cfg.f_m` | m | `429e-3` | `config/default_config.m` | RD beta、FFT 焦平面坐标 |
| DOE 到透镜距离 | `cfg.doe_to_lens_m` | m | `200e-3` | `config/default_config.m` | full propagation 使用；初始 RD 解析公式本身不使用 |
| clear aperture 直径 | `cfg.aperture_diameter_m` | m | `15e-3` | `config/default_config.m` | 生成 `aperture_mask`，限制相位有效区域 |
| 入射光 1/e² 强度直径 | `cfg.input_1e2_diameter_m` | m | `5e-3` | `config/default_config.m` | 生成高斯入射振幅 |
| 入射光 1/e 振幅半径 | `cfg.input_1e_radius_m` | m | `input_1e2_radius_m/sqrt(2)` | `config/default_config.m` | RD 公式中的 `ri_m` |
| 目标 x 尺寸 | `cfg.target_size_x_m` | m | `330e-6` | `config/default_config.m` | 进入 `Ro_x_m = target_size_x_m/sqrt(pi)` |
| 目标 y 尺寸 | `cfg.target_size_y_m` | m | `120e-6` | `config/default_config.m` | 进入 `Ro_y_m = target_size_y_m/sqrt(pi)` |
| RD x 输出尺度 | `cfg.Ro_x_m` | m | `target_size_x_m/sqrt(pi)` | `config/default_config.m` | RD beta x 方向输出尺度 |
| RD y 输出尺度 | `cfg.Ro_y_m` | m | `target_size_y_m/sqrt(pi)` | `config/default_config.m` | RD beta y 方向输出尺度 |
| 网格尺寸 | `cfg.N` | samples | `2048` | `config/default_config.m` | DOE 相位矩阵尺寸 |
| 目标焦平面采样 | `cfg.requested_focal_dx_m` | m/pixel | `2.5e-6` | `config/default_config.m` | 反推 DOE 网格范围/采样 |
| DOE 网格范围 | `cfg.doe_grid_extent_m` | m | `lambda*f/requested_focal_dx` | `config/default_config.m` | `N*dx_doe_m` |
| DOE 采样间隔 | `cfg.dx_doe_m` | m/pixel | `doe_grid_extent_m/N` | `config/default_config.m` | DOE 坐标和 FFT 采样 |
| 实际焦平面采样 | `cfg.focal_dx_m` | m/pixel | `lambda*f/(N*dx_doe)` | `config/default_config.m` | 焦平面 x/y 坐标间隔 |
| 相位符号 | `cfg.phase_sign` | - | `1` | `config/default_config.m` | Fourier sign convention |
| x/y 相位倍率 | `cfg.phase_scale_x/y` | - | `1 / 1` | `config/default_config.m` | 显式相位缩放，默认不改变 RD 结果 |

## 坐标系统说明

- DOE 面坐标来自 `make_grid.m`：
  - `x_m = (-N/2:N/2-1) * cfg.dx_doe_m`
  - `y_m` 同上
  - `[X_m, Y_m] = meshgrid(x_m, y_m)`
- DOE 面单位是米，图像中常转换为 mm 或 µm。
- 15 mm 是 clear aperture / pupil 直径，由 `aperture_mask = r <= 7.5 mm` 表示。
- 5 mm 是入射高斯光 1/e² 强度直径，不是 clear aperture。
- 焦平面坐标由 Fourier lens 关系决定：
  - `focal_dx = lambda * f / (N * dx_doe)`
  - 默认约为 `2.5 µm/pixel`。
- `mraf_forward_fft.m` 使用正交 FFT 归一化用于 MRAF 内部迭代；`fresnel_to_focal_fft.m` 使用带物理尺度的 Fourier lens FFT 用于 RD baseline 诊断。

## 初始相位输出说明

- `phase0` 尺寸：默认 `2048 x 2048`。
- 单位：弧度。
- `phase_unwrapped_rad` / `phase_rd_unwrapped`：未包裹相位，当前传播和 MRAF 使用它构造 `exp(1i*phase0)`。
- `phase_wrapped_rad` / `phase_rd_wrapped`：包裹到 `[0, 2π)`，主要用于显示、导出、DOE 相位图。
- 在 `run_mraf_one.m` 中：
  - `phase_rd_unwrapped` 传给 `run_mraf_refinement(..., phase_rd_unwrapped, ...)`。
  - `common/phase_rd_wrapped.mat` 保存变量 `phase_rd_wrapped`。
  - MRAF 后保存 `common/phase_mraf_wrapped.mat`。
- 在 `run_one.m` 中：
  - `save_artifacts.m` 保存 `ideal_phase_unwrapped.mat`、`ideal_phase_wrapped.mat`、`full_phase_unwrapped.mat`、`full_phase_wrapped.mat`。
- 在 `export_phase_for_python.m` 中：
  - 输出 `exports/phase_rd.mat`。
  - 变量包括 `phase_rd_wrapped`、`phase_rd_unwrapped`、`phase_wrapped`、`phase_unwrapped`、`phase_rad`、`export_meta`。
- 在 `run_initial_phase_only.m` 中：
  - 输出 `artifacts/initial_phase_only/<timestamp>/phase0.mat`。
  - 变量包括 `phase0_unwrapped_rad`、`phase0_wrapped_rad`、`phase_info`、`x_m/y_m`、`focal_x_m/focal_y_m`。

## 初始相位和 MRAF 的关系

- `phase0` 是 MRAF 的初值，不是 MRAF 的目标。
- MRAF 的第一步用 `input_amplitude .* exp(1i * phase0)` 做 RD baseline forward FFT。
- 后续 MRAF 每轮都会修改 DOE 相位，因此最终 `result.mraf.phase` 不等于 `phase0`。
- 如果 `phase0` 本身由有限 beta 的 RD 解析式产生了较平滑的边缘，MRAF 后续可以改善均匀性或局部能量分布，但很难无代价地凭空生成无限陡边；这也是之前 shoulder/transition 调参中需要注意的物理限制。

## 当前可能存在的疑点

- `run_one.m` 使用 `fresnel_to_focal_fft.m` 的物理尺度 FFT；MRAF 内部使用 `mraf_forward_fft.m` 的正交 FFT。两者坐标一致由外部 `focal_x_m/focal_y_m` 保证，但归一化方式不同。
- `cfg.aperture_diameter_m = 15 mm` 是 clear aperture，不是入射光斑；入射光斑为 `cfg.input_1e2_diameter_m = 5 mm`。
- `cfg.target_size_x_m/y_m = 330/120 µm` 同时影响 RD 初始相位的 `Ro_x/Ro_y`，也常被 MRAF/metrics 用作目标尺寸参考。不同 target_mode 可能另有 mask 定义，但不改变 phase0 的生成。
- `cfg.requested_focal_dx_m` 通过反推 DOE 网格范围固定焦平面采样；如果修改该值，会改变 `dx_doe_m` 和相位采样。
- 当前 MATLAB 主路径中没有发现真正的迭代 point-to-point 相位求解器；项目名和注释中有 point-to-point/RD 说法，但实际 `phase0` 核心是 separable RD 解析公式。

## 如何单独运行

在 MATLAB 中：

```matlab
cd E:\program\Point2P\DOE_ROMERO_DICKEY_MATLAB
run_initial_phase_only
```

输出在：

```text
artifacts/initial_phase_only/<timestamp>/
```
