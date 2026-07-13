# SLM 偏移补偿说明

这套脚本对应的是实验里常说的：

- SLM 实际安装位置有偏移
- 不改输入高斯光束位置
- 直接在原始 `2048 x 2048` 大相位上做整数像素平移补偿
- 然后再加 blaze，把焦斑偏到实验需要的位置

## 相关脚本

- `shift_sweep_slm.py`
  - 主脚本
  - 对大相位做 `X/Y` 像素位移
  - 再加 blaze
  - 输出 SLM 用文件

- `verify_shift_sweep.py`
  - 只做仿真验证
  - 用来扫不同位移量，比较上下不均匀、中心 profile 和 RMS

- `add_blaze_grating.py`
  - 单独给相位加 blaze
  - 不负责做 SLM 位置补偿扫描

## 核心逻辑

真正的 SLM 位置补偿在 `shift_sweep_slm.py` 里：

```python
ph_shifted = shift_phase(phase, sy, sx)
```

其中 `shift_phase()` 本质是：

```python
np.roll(np.roll(phase, shift_y_px, axis=0), shift_x_px, axis=1)
```

也就是把原始大相位整体平移，而不是改目标，不是改高斯光束，也不是只加 blaze。

后续再执行：

```python
ph_blazed = add_blaze(ph_shifted, dx_doe, lam, f_m, blaze_x_um, blaze_y_um)
```

最终流程是：

`原始 phase`
-> `整体位移补偿`
-> `加 blaze`
-> `导出 SLM 文件`

## 方向约定

在 `shift_sweep_slm.py` 里，约定是：

- `+Y shift`：相位在 SLM 上向下移动
- `-Y shift`：相位在 SLM 上向上移动
- `+X shift`：相位在 SLM 上向右移动
- `-X shift`：相位在 SLM 上向左移动

脚本里写的物理解释是：

- `+Y shift` 更适合补偿 “上亮下暗”
- `-Y shift` 更适合补偿 “下亮上暗”

注意：这里补偿的是“相位图相对固定光束的位置”。

## 推荐使用顺序

### 1. 先做仿真筛选

先用 `verify_shift_sweep.py` 扫一组位移，找大概范围。

示例：

```powershell
python verify_shift_sweep.py artifacts\_repro_check\phase_refined.npy --shifts-y -20,-15,-10,-5,0,5,10,15,20 --shifts-x 0 --out shift_verify_output
```

重点看：

- `shift_sweep_results.json`
- `shift_sweep_verification.png`
- `shift_sweep_profiles.png`

## 2. 再生成实验文件

用 `shift_sweep_slm.py` 生成真正给 SLM 的多组相位文件。

示例：

```powershell
python shift_sweep_slm.py artifacts\_repro_check\phase_refined.npy --out shift_sweep_output --shifts-y -20,-15,-10,-5,0,5,10,15,20 --shifts-x 0 --blaze-x 200 --blaze-y 200
```

输出里每一组通常都会包含：

- 位移后的相位记录
- 加 blaze 后的完整相位
- 给 SLM 用的输出文件
- 对应的参数摘要

## 建议基线

如果你当前就是沿用这次测试目录里的基线，建议从这里开始：

- 相位来源：`artifacts/_repro_check/phase_refined.npy`
  - 或者你确认好的最终 `phase_refined.npy`
- 先只扫 `Y` 方向
- `X` 先固定为 `0`
- blaze 先保持当前常用值

## 区分三个概念

- `shift_phase`
  - 这是 SLM 位置补偿
  - 改的是大相位阵列的位置

- `add_blaze`
  - 这是焦斑偏转
  - 改的是线性相位斜坡

- `swap_phase_xy`
  - 这是导入 MAT 时的轴交换
  - 不是实验补偿，不是 SLM 位移

## 当前目录用途

`rtad_mraf_gs_python_test_20260605` 现在可以直接用来做这类测试性工作。

如果后面要改补偿策略，优先看这两个文件：

- `shift_sweep_slm.py`
- `verify_shift_sweep.py`
