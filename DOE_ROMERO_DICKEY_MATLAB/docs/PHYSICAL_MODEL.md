# 物理模型说明

## 光路平面

本项目包含三个主要平面：

- DOE 面：加载 phase-only 相位，同时受到入射高斯光和 15 mm clear aperture 限制。
- Fourier lens / 场镜：焦距 `f = 429 mm`。
- 焦平面：目标矩形平顶所在平面。

本项目明确分两步验证，不把所有真实传播路径一开始混在一起：

1. `ideal_`：ideal Fourier-lens FFT baseline。假设 DOE 面就是有效 Fourier transform pupil plane，直接用 FFT 得到焦平面强度。这一步用于验证 Romero-Dickey separable phase 本身是否合理。
2. `full_`：full propagation model。在 baseline 正确后，再加入 DOE 到 lens 的 `200 mm` 自由空间传播、透镜面圆形口径、薄透镜二次相位、lens 到焦平面的传播。这一步用于评估真实光路偏离理想 Fourier 条件后的平顶质量、边缘 overshoot、transition width 和效率变化。

## Fourier lens 焦平面关系

ideal baseline 使用：

```text
U_focal(xf, yf) proportional to FFT{U_doe(x, y)}
```

频率坐标为：

```text
fx = (-N/2 : N/2-1) / (N*dx_doe)
xf = lambda * f * fx
```

因此焦平面采样为：

```text
focal_dx = lambda * f / (N * dx_doe)
         = lambda * f / Lgrid
```

这里 `Lgrid = N*dx_doe` 是数值 FFT 计算窗口。它不是入射光斑直径，也不是 clear aperture 直径。

## Full propagation model

full model 使用以下步骤：

```text
U_lens_in = AngularSpectrum(U_doe, z = doe_to_lens_m)
U_lens_clipped = U_lens_in * lens_aperture_mask
U_lens_out = U_lens_clipped * exp(-i*k*(x^2+y^2)/(2*f))
U_focal = FresnelSingleFFT(U_lens_out, z = f)
```

如果透镜口径设为无限大，理想标量自由传播主要引入频域相位，焦平面强度可能与 ideal FFT baseline 几乎相同。为了让 full model 能反映实际 200 mm 传播后的截断影响，默认 `apply_lens_aperture_in_full_model = true`，且 `lens_aperture_diameter_m = 15 mm`。

该模型仍是标量衍射近似、薄透镜近似，不包含像差、倾斜入射、偏心或实验误差。

## Phase-only 元件

DOE 被建模为无吸收的 phase-only 元件：

```text
U_after_DOE = U_input * exp(i * phase_rad)
```

输入能量只通过 aperture mask 被截断；DOE 本身不改变振幅。用于显示或器件加载的相位可以 wrap 到 `[0, 2*pi)`，但传播中使用 unwrapped 或 wrapped 相位在数学上等价。

## Romero-Dickey 解析相位

Romero-Dickey 设计基于 stationary phase / geometrical optics 近似，将输入 Gaussian 能量映射到 flat-top 目标。

本项目采用用户指定的一维解析式：

```text
xi = x / ri
phi_dimless(xi) = xi*sqrt(pi)/2*erf(xi) + 0.5*exp(-xi^2) - 0.5
phase_rad = beta * phi_dimless
beta = 2*pi*ri*Ro/(lambda*f)
```

矩形目标不是圆对称目标，因此不使用 radial RotationalPhaseDesign，而是采用 separable approximation：

```text
phase_2d(x, y) = phase_x(x) + phase_y(y)
```

这种近似适合作为矩形平顶的解析初始相位，但不保证直接得到完美矩形均匀光斑。
