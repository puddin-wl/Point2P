# 参数约定

## 入射高斯光斑

用户给定：

```text
input_1e2_diameter_m = 5 mm
```

这是强度下降到 `1/e^2` 的直径，不是 DOE aperture。

因此：

```text
w = input_1e2_radius_m = 2.5 mm
I(r) = exp(-2*r^2/w^2)
amplitude(r) = exp(-r^2/w^2)
```

Romero-Dickey 论文中的 `ri` 是输入强度下降到 `1/e` peak 的半径。对上述高斯强度：

```text
ri = w / sqrt(2)
   = 2.5 mm / sqrt(2)
   ≈ 1.7678 mm
```

代码中明确区分：

- `input_1e2_radius_m`：`1/e^2` 强度半径。
- `input_1e_radius_m`：论文 beta 使用的 `1/e` 强度半径 `ri`。

## Clear aperture

```text
aperture_diameter_m = 15 mm
```

这是 DOE 圆形通光孔径，不是入射光斑直径。

入射光是 `5 mm 1/e^2` 强度直径高斯光，被 `15 mm` DOE clear aperture 截断。由于 `15 mm` 远大于 `5 mm` 光斑直径，默认截断主要影响高斯远尾和数值 aperture 边界。

## 目标矩形尺寸

目标图样为：

```text
330 um x 120 um
```

这个尺寸按 50% intensity boundary 定义。代码中：

```text
target_half_x_m = 165 um
target_half_y_m = 60 um
```

## Ro 定义

默认 separable approximation 使用：

```text
Ro_x = target_size_x_m / sqrt(pi)
Ro_y = target_size_y_m / sqrt(pi)
```

也就是：

```text
Ro_x = 330 um / sqrt(pi) ≈ 186.18 um
Ro_y = 120 um / sqrt(pi) ≈ 67.70 um
```

如果需要尝试别的 `Ro` 定义，可在 `config/default_config.m` 中修改：

```matlab
cfg.Ro_definition = "...";
cfg.Ro_x_m = ...;
cfg.Ro_y_m = ...;
```

## beta 计算

Romero-Dickey beta 为：

```text
beta = 2*pi*ri*Ro/(lambda*f)
```

默认参数下：

```text
lambda = 532 nm
f = 429 mm
ri ≈ 1.7678 mm
Ro_x ≈ 186.18 um
Ro_y ≈ 67.70 um
```

得到：

```text
beta_x ≈ 9.061
beta_y ≈ 3.295
```

`beta_y` 较小，说明 y 方向目标更窄，stationary-phase 近似更吃力，边缘锐度、均匀性和 ringing 更可能成为限制。

## 数值窗口与焦平面采样

默认 `N = 2048`，希望焦平面采样约 `2.5 um/pixel`。根据：

```text
focal_dx = lambda*f/(N*dx_doe) = lambda*f/Lgrid
```

得到：

```text
Lgrid ≈ 91.2912 mm
dx_doe ≈ 44.5758 um
```

这里 `Lgrid` 只是 FFT 计算窗口，不是 5 mm 光斑，也不是 15 mm aperture。

