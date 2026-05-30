# 代码导读 — 核心算法逐段解释

> 读完这份文档, 你就能理解从初始相位到 SLM 输出的每一行关键代码在做什么。
> 所有代码片段均从本目录的 `.py` 文件中直接提取。

---

## 1. 管线数据流

```
make_phase0.py          run_rtad_mraf_gs_case.py     convert_to_slm.py
(Stage 1)               (Stage 2)                    (Stage 3)

beam=3.8mm              phase0.mat                   2048×2048 相位
    │                       │                             │
    ▼                       ▼                             ▼
Romero-Dickey  ──────►  WGS 迭代 200轮  ──────►  复振幅插值 → SLM BMP
解析相位                 (MRAF投影+权重更新)          1024×1024 @ 17μm
```

---

## 2. 目标生成 — 330×120μm 矩形平顶

**文件**: `src/rtad_target.py`

### 2.1 升余弦边缘

矩形不是硬边界, 而是通过升余弦 (raised cosine) 从 flat 区平滑过渡到零:

```python
# src/rtad_target.py, line 69-78
def raised_cosine_edge(u, u0, u1):
    """一维升余弦下降沿: u≤u0→1, u≥u1→0, 中间余弦过渡."""
    C = np.zeros_like(u, dtype=np.float32)
    C[u <= u0] = 1.0                    # flat 区 = 1
    idx = (u > u0) & (u < u1)
    t = (u[idx] - u0) / (u1 - u0)      # 归一化到 [0,1]
    C[idx] = 0.5 * (1.0 + np.cos(np.pi * t))  # 余弦曲线从 1 降到 0
    return C
```

**物理含义**: 这是一维 profile。对于 330×120μm 的目标:
- `a50 = 165μm` (X 方向半宽, 50% 强度处)
- `a0 = a50 - delta_x = 150μm` (flat 区边界, 强度=1)
- `a1 = a50 + delta_x = 180μm` (过渡区外边界, 强度=0)
- `delta_x = 15μm` 控制边缘过渡的陡峭程度

二维模板 = `I(x,y) = C_x(|x|) × C_y(|y|)`:

```python
# src/rtad_target.py, line 151-154
Ix = raised_cosine_edge(absX, a0, a1)
Iy = raised_cosine_edge(absY, b0, b1)
I_full = np.clip(Ix * Iy, 0.0, 1.0)   # 二维强度模板 (可分离)
A_full = np.sqrt(I_full)               # 振幅 = sqrt(强度)
```

### 2.2 四层 Mask 区域

MRAF 算法的核心是把焦面分成四个区, 每个区有不同约束规则:

```python
# src/rtad_target.py, line 161-171
# 从内到外四个区域:
mask_flat   = (absX <= a0) & (absY <= b0)          # 平顶核心 (150×52μm)
mask_signal = I_full >= release_level               # 约束区 (I ≥ e⁻²)
mask_signal = mask_signal | mask_flat               # 确保 flat 一定在 signal 内
mask_free   = mask_guard_window & ~mask_signal      # 自由区 (可随意演化)
mask_bg_far = ~(mask_signal | mask_free)            # 远背景 (轻微衰减)
```

**示意图** (以 X 方向为例):

```
强度
1.0 ┤  ████████████
    ┤  █  flat   █╲
0.5 ┤  █  (a0=150) █╲___ 升余弦过渡 ___
    ┤  █           █        ╲
e⁻² ┤  ██████████████████████╲_________
    ┤  │← signal →│              ╲
0   ┤  │←          guard_window         →│
    ┤  │←    free    →│←  bg_far  →│

 ←─ 区域边界: a0=150, a50=165, a1=180, a2=200μm ─→
```

---

## 3. 初始相位 — Romero-Dickey 解析法

**文件**: `make_phase0.py`

DOE 初始相位不能用随机值——那样迭代根本收不敛。Romero-Dickey 方法从高斯光束和矩形平顶的能量守恒关系出发, 解析求出可分离的初始相位:

```python
# make_phase0.py, line 44-51
def romero_dickey_phase_1d(x_m, ri_m, Ro_m, lambda_m, f_m):
    """Romero-Dickey 一维解析相位."""
    xi = x_m / ri_m                           # 归一化坐标 (ri = 1/e 振幅半径)
    # 无量纲相位函数:
    phi_dimless = (xi * np.sqrt(np.pi) / 2.0 * erf(xi)
                   + 0.5 * np.exp(-xi**2) - 0.5)
    # β: 无量纲参数, 表示几何光学近似程度
    beta = 2.0 * np.pi * ri_m * Ro_m / (lambda_m * f_m)
    # β × φ_dimless → 真实相位 (rad)
    return beta * phi_dimless
```

**β 的物理含义**:
```
β = 2π × (光束半径) × (目标半宽) / (λ × 焦距)
```
- β > 10 → 几何光学区, 稳相近似精确, WGS 只需小幅修正
- β < 10 → 衍射区, 初始相位偏差大, WGS 需要大幅修正

对于 beam=3.8mm, f=200mm: `β_x=14.8, β_y=5.4`。X 方向在几何光学区, Y 方向在衍射区边缘。

二维相位是可分离的: `φ(x,y) = φ_x(x) + φ_y(y)`

---

## 4. DOE 面的输入光场

**文件**: `src/propagation.py`

```python
# src/propagation.py, line 48-66
def make_input_gaussian(shape, dx_doe_m, gaussian_1e2_diameter_m,
                        clear_aperture_m, xp, dtype):
    """构建 DOE 平面高斯输入振幅."""
    Ny, Nx = shape
    # 坐标网格 (axis 0 = y, axis 1 = x)
    x = (xp.arange(Nx, dtype=dtype) - Nx//2) * dtype(dx_doe_m)
    y = (xp.arange(Ny, dtype=dtype) - Ny//2) * dtype(dx_doe_m)
    X, Y = xp.meshgrid(x, y)
    r2 = X*X + Y*Y

    w = dtype(gaussian_1e2_diameter_m / 2.0)  # 1/e² 强度半径
    amp = xp.exp(-r2 / (w * w))               # 高斯振幅 A = exp(-r²/w²)

    # 圆形通光孔径截断 (直径 15mm)
    aperture_radius = dtype(clear_aperture_m / 2.0)
    amp = xp.where(r2 <= aperture_radius**2, amp, 0.0)

    # 总功率归一化到 1
    return normalize_power(amp, xp, target_power=1.0)
```

**关键点**: `A(r) = exp(-r²/w²)`——注意是强度 I 的 1/e² 半径, 对应振幅的 1/e 半径。`w = 1/e²强度直径 / 2 = 1.9mm` (对 beam=3.8mm)。

---

## 5. 傅里叶透镜传播

**文件**: `src/propagation.py`

DOE 面和焦面之间通过一个傅里叶透镜耦合。前向传播 = 从 DOE 到焦面 (我们看到的), 反传 = 从焦面回到 DOE (迭代中用来更新相位)。

```python
# src/propagation.py, line 13-16
def forward_fft(field_in, xp):
    """DOE 平面 → 焦平面: FFT."""
    # ifftshift 把零频移到 (0,0) → fft2 → fftshift 把零频移回中心
    return xp.fft.fftshift(xp.fft.fft2(xp.fft.ifftshift(field_in), norm="ortho"))

def backward_fft(field_out, xp):
    """焦平面 → DOE 平面: 逆 FFT."""
    return xp.fft.fftshift(xp.fft.ifft2(xp.fft.ifftshift(field_out), norm="ortho"))
```

**关键约定**:
- `norm="ortho"`: 保证 FFT 前后总功率守恒
- 矩阵方向: `axis 0 = y (行)`, `axis 1 = x (列)`——和 MATLAB `meshgrid` 一致
- 焦面坐标: `x_focal[k] = (k - N/2) × focal_dx`, 其中 `focal_dx = 2.5μm`

**焦面像素间距**: `focal_dx = λ × f / (N × dx_doe) = 2.5μm`, 2048×2048 的焦面覆盖 ±2.56mm。

---

## 6. MRAF 焦面投影 — 三个区域, 三种规则

**文件**: `src/mraf_gs.py`

这是每次迭代里最关键的一步: 把计算出的焦面复振幅 "投影" 到物理约束上。

```python
# src/mraf_gs.py, line 287-326 (精简, 去掉已注释的 GS 分支)
def _project_farfield(farfield, target_amp_eff, masks_b, xp, method,
                       mraf_factor, bg_mode, bg_factor):
    """MRAF 焦面投影."""
    # 提取当前焦面相位
    phase_ff    = xp.angle(farfield)
    phase_factor = xp.exp(1j * phase_ff)
    signal = masks_b["mask_signal"]

    projected = farfield.copy()
    # ① signal 区: 保留相位, 替换振幅为目标振幅
    projected[signal] = target_amp_eff[signal] * phase_factor[signal]

    # ② free 区: 保留相位, 振幅 × mraf_factor (0.4)
    free = masks_b["mask_free"]
    projected[free] *= mraf_factor

    # ③ 远背景: keep / attenuate / zero
    bg = masks_b.get("mask_bg_far", ...)
    if bg_mode == "attenuate":
        projected[bg] *= bg_factor        # × 0.9, 几乎等同 keep
    elif bg_mode == "zero":
        projected[bg] = 0

    return projected
```

**三区处理总结**:

| 区域 | 振幅处理 | 相位处理 | 效果 |
|------|---------|---------|------|
| `mask_signal` | 替换为 target_amp × WGS_weight | 保留原相位 | 强制趋近目标平顶 |
| `mask_free` | × 0.4 (衰减) | 保留原相位 | 给噪声自由度, 不完全压制 |
| `mask_bg_far` | × 0.9 (几乎不变) | 保留原相位 | 不做强约束, 避免能量倒灌 |

**为什么 mraf_factor=0.4 而不是 0**: 如果 free 区完全置零, 能量无处可去, 会回流到 signal 区形成波纹。给 free 区一定自由度反而让平顶更平。

---

## 7. WGS 权重更新 — 平顶内的局部修匀

**文件**: `src/mraf_gs.py`

MRAF 只做区域级约束。WGS 在 `mask_flat` 内部做像素级权重反馈, 进一步压平不均匀。

```python
# src/mraf_gs.py, line 164-213 (精简)
def _update_flat_wgs_weights(weights, farfield_amp, update_region,
                              xp, backend, feedback_exponent,
                              clip_min, clip_max, normalize_weights):
    """WGS flat_local 权重更新."""
    eps = np.float32(1e-12)

    # 1. 提取 mask_flat 内的焦面振幅
    amp_flat = farfield_amp[update_region]

    # 2. 计算 flat 区平均振幅
    amp_mean = xp.mean(amp_flat)

    # 3. 核心公式: 每个像素的修正系数
    #    ratio = mean(|E_flat|) / |E_pixel|
    #    如果某像素比平均暗 → ratio > 1 → 权重增大 → 下次投影时要求更高
    #    如果某像素比平均亮 → ratio < 1 → 权重减小 → 下次投影时要求更低
    ratio = amp_mean / xp.maximum(amp_flat, eps)
    updated = weights[update_region] * xp.power(ratio, feedback_exponent)

    # 4. Clip: 防止权重跑飞 (默认 0.5 ~ 1.5)
    updated = xp.clip(updated, clip_min, clip_max)

    # 5. 均值归一化: 保持所有权重的平均值为 1.0
    if normalize_weights:
        updated = updated / xp.mean(updated)

    weights[update_region] = updated
    return weights
```

**直观理解**: 假设 mask_flat 内某个像素振幅是 mean 的 0.8 倍, exponent=0.8:
```
ratio = 1.0 / 0.8 = 1.25
weight_new = weight_old × 1.25^0.8 = weight_old × 1.20
```
这个像素的权重增加了 20%。下次 MRAF 投影时, 该像素被要求匹配更高目标振幅 → 迭代推动它变亮 → 趋近平坦。

---

## 8. 主迭代循环 — 完整流程

**文件**: `src/mraf_gs.py`, `run_refinement()` 函数

每轮迭代做 5 件事:

```python
# src/mraf_gs.py, line 538-609 (精简, 去掉注释和日志)
for it in range(1, total_iters + 1):          # total_iters = 200

    # ① 前向传播: DOE → 焦面
    field   = input_amp * exp(1j * phase)     # DOE 面复振幅
    farfield = forward_fft(field)             # FFT → 焦面复振幅
    farfield_amp = abs(farfield)              # 焦面振幅

    # ② WGS 权重更新 (每 5 轮一次)
    if it % update_every == 0:
        weights = _update_flat_wgs_weights(   # 只在 mask_flat 内更新
            weights, farfield_amp, update_region,
            exponent=0.8, clip=[0.5, 1.5], normalize=True)

    # ③ 加权目标振幅
    target_weighted = base_target.copy()
    target_weighted[mask_flat] *= weights     # flat 区 × WGS 权重

    # ④ MRAF 投影 + 反传 + 相位提取
    projected  = _project_farfield(           # 三区投影
        farfield, target_weighted,
        mraf_factor=0.4, bg_factor=0.9)
    nearfield  = backward_fft(projected)      # IFFT → DOE 面
    phase      = angle(nearfield) % 2π        # 提取新相位

    # ⑤ 指标记录 (每 10 轮)
    if it % metrics_interval == 0:
        记录 RMS, size50, 效率...
```

**数据流示意图**:

```
  DOE 面                        焦面                    DOE 面
┌─────────┐   ① FFT    ┌─────────────┐  ④ IFFT   ┌─────────┐
│ φ_old   │ ────────→  │ E_farfield   │ ────────→ │ φ_new   │
│ A_input │            │              │           │         │
└─────────┘            │ ② WGS权重 ↑  │           └─────────┘
                       │ ③ 加权目标  │
                       │ ④ MRAF投影  │
                       └─────────────┘
```

---

## 9. ROM → SLM — 复振幅插值

**文件**: `convert_to_slm.py`

计算用的 2048×2048 相位不能直接送到 SLM (1024×1024)。需要裁切 + 缩放到 SLM 分辨率和像素间距。

```python
# convert_to_slm.py, line 84-106 (精简)
# 1. 根据 SLM 物理尺寸 (17.4mm) / 计算网格像素间距 (20.8μm) → 裁切 837×837 像素
phys_mm  = 1024 * 17.0e-3                     # SLM 物理尺寸 = 17.4mm
crop_px  = int(phys_mm / (dx_doe_um * 1e-3))  # 裁切像素数 = 837
phase_crop = phase[cy-418:cy+419, cx-418:cx+419]  # 从 2048 中心裁出

# 2. ★ 关键: 复振幅插值 (不是直接插值相位!)
zoom_ratio = 1024 / 837                        # ≈ 1.22x 放大
cfield = np.exp(1j * phase_crop)               # 转为复振幅
real_z = zoom(cfield.real, zoom_ratio, order=3)  # 实部 cubic 插值
imag_z = zoom(cfield.imag, zoom_ratio, order=3)  # 虚部 cubic 插值
phase_slm = np.arctan2(imag_z, real_z)          # 从插值后的复振幅恢复相位
phase_slm = phase_slm % (2π)
```

**为什么不能直接插值包裹相位**: 包裹相位在 0 ↔ 2π 处有跳变, cubic spline 会把跳变当成高频信号产生巨大振铃。先转成连续平滑的 `exp(iφ)` 分量, 插值完再恢复相位, 完全避免这个问题。

---

## 10. 闪耀光栅 — 让光斑偏离中心

在 DOE 相位上叠加线性斜坡, 焦面光斑就会平移。

```python
# 光栅相位: φ_grating(x) = 2π × (Δx / focal_dx) × x / N
focal_dx_um = 2.5      # 焦面像素间距
shift_x_um  = 400.0    # 想偏移 400μm
shift_y_um  = 200.0    # 想偏移 200μm

periods_x = shift_x_um / focal_dx_um   # = 160 个周期 (跨 2048 像素)
periods_y = shift_y_um / focal_dx_um   # = 80  个周期

x = np.arange(N) - N//2
X, Y = np.meshgrid(x, x)
grating = 2π * (periods_x * X + periods_y * Y) / N

phase_with_grating = (phase_doe + grating) % (2π)
```

**物理原理**: 线性相位斜坡 = 平面波倾斜入射 = 焦面光斑侧移。1 个光栅周期 (跨整个 DOE 孔径) = 焦面偏移 1 个像素 = 2.5μm。

---

## 11. 指标计算 — 怎么判断做得好不好

**文件**: `src/metrics.py`

```python
# src/metrics.py, line 86-122 (精简)
def compute_metrics(intensity, masks, x_um, y_um, iteration):
    I = intensity
    mask_flat = masks["mask_flat"]

    # ① 归一化: I_norm = I / mean(I[mask_flat])
    #    不以全局最大值为基准, 因为那会被偶然的 hot pixel 拉偏
    I_n = I / np.mean(I[mask_flat])

    # ② RMS 不均匀度 (最核心指标)
    flat_vals = I_n[mask_flat]
    flat_rms  = np.std(flat_vals) / np.mean(flat_vals)
    # 例如: std=0.0087, mean=1.0 → RMS = 0.87%

    # ③ size50: 中心剖面半高全宽
    profile_x = I_n[N//2, :]
    size50_x  = width_at_level(x_um, profile_x, level=0.5)
    # 从中心向两侧找到 I=0.5 的 crossing 点, 线性插值

    # ④ 效率: signal 区内能量 / 总能量
    efficiency = np.sum(I[mask_signal]) / np.sum(I)

    return { "flat_rms": flat_rms, "size50_x": size50_x, ... }
```

**RMS 定义**: `std(mask_flat) / mean(mask_flat)`, 只算 flat core 内的像素。对于 beam=3.8mm, 200 轮后 RMS ≈ 0.71%。

---

## 12. 参数速查表

| 参数 | 值 | 含义 |
|------|-----|------|
| λ | 532 nm | 激光波长 |
| f | 200 mm | 透镜焦距 |
| N | 2048 | 计算网格 |
| focal_dx | 2.5 μm | 焦面像素间距 |
| beam (实测) | 3.8 mm | 1/e² 强度直径 |
| 目标 W50×H50 | 330×120 μm | 平顶 50% 尺寸 |
| δ_x / δ_y | 15 / 8 μm | 边缘过渡半宽 |
| mraf_factor | 0.4 | free 区衰减系数 |
| bg_factor | 0.9 | 背景衰减 (接近 keep) |
| feedback_exponent | 0.8 | WGS 修正力度 |
| weight_clip | [0.5, 1.5] | 权重裁剪范围 |
| update_every | 5 | 每 N 轮更新一次权重 |
| num_iters | 200 | 总迭代轮数 |

---

## 13. 文件索引

| 想了解什么 | 去看哪个文件 |
|-----------|------------|
| 初始相位怎么来的 | `make_phase0.py` (RD 公式) |
| 目标光斑怎么定义 | `src/rtad_target.py` (升余弦 + mask) |
| FFT 传播怎么写的 | `src/propagation.py` |
| 迭代循环逻辑 | `src/mraf_gs.py` → `run_refinement()` |
| 焦面投影公式 | `src/mraf_gs.py` → `_project_farfield()` |
| WGS 权重怎么更新 | `src/mraf_gs.py` → `_update_flat_wgs_weights()` |
| 怎么转成 SLM 格式 | `convert_to_slm.py` |
| 光斑质量怎么量化 | `src/metrics.py` → `compute_metrics()` |
| CLI 参数怎么传 | `run_rtad_mraf_gs_case.py` → `parse_args()` |
| GPU/CPU 怎么选 | `src/backend.py` |
| MATLAB 文件怎么读 | `src/io_mat.py` |
