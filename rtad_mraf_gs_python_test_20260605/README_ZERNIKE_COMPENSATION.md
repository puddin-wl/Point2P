# V2 球差与离焦补偿相位导出

`export_zernike_compensated_slm.py` 在已经确认的 V2 相位上叠加标准
Noll-Zernike 球差和离焦补偿，并输出实验可直接加载的 1024×1024、8-bit BMP。

## 当前默认补偿

```text
Z40 = +0.10625 RMS waves
Z20 = +0.25000 RMS waves
pupil diameter = 15 mm
```

这组正号参数是仿真得到的负球差/负离焦主案例的反号，作为第一张实验补偿相位。
它不是最终标定值，后续应根据真实光斑指标联动微调 Z40 和 Z20。

## 处理顺序

```text
V2 phase_refined
→ V2安装补偿平移
→ 以光轴为中心叠加Z40/Z20补偿
→ 叠加X/Y闪耀光栅
→ SLM物理窗口裁切和复相位插值
→ 1024×1024、8-bit BMP
```

Zernike 补偿在安装平移之后加入，因此安装补偿只移动矩形整形相位，而球差
补偿仍以实际光束/场镜轴心为中心。

15 mm 是 Zernike 系数的归一化 pupil，并不是 6.5 mm Gaussian 光斑边界。
程序在 `rho<=1` 内严格保持原 Z40/Z20 定义；在 pupil 外默认用宽度
`0.20 rho`（物理宽度 1.5 mm）的 C2 连续径向延拓，再保持为常量活塞相位。
因此不会在 SLM 上生成“15 mm 圆外突然清零”的人工圆形相位跳变。

## 当前基准命令

```powershell
D:\software\anaconda\envs\slmrtad\python.exe export_zernike_compensated_slm.py `
  artifacts\run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm\phase_refined.npy `
  --out artifacts\run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm\SLM_LOAD_ZERNIKE_COMP_BASELINE_20260723 `
  --z40 0.10625 --z20 0.25000 `
  --extension-width-rho 0.20 `
  --install-shift-x 5 --install-shift-y 5 `
  --blaze-x 200 --blaze-y 200 `
  --verify-zero-reference artifacts\run_w50_375p020_h50_137_target_expX330_v2_beam6p5mm\SLM_LOAD_20260715\phase_20260715_targetExpX330_v2_shiftX+5_Y+5.bmp
```

`--verify-zero-reference` 会用同一程序重建零 Zernike 的历史 X+5/Y+5 相位，
要求与旧 BMP 逐像素完全一致，以确认新程序没有改变既有导出链路。

## 输出

- `phase_V2_posSphericalComp_*.bmp`：实验直接加载文件；
- 同名 `.npy`、`.mat`：SLM 相位数值记录；
- `zernike_compensation_waves_2048.npy`：补偿波前；
- `phase_compensated_pre_blaze_2048.npy`：加闪耀光栅前的完整相位；
- `phase_compensated_blazed_2048.npy`：加闪耀光栅后的完整相位；
- `zernike_compensation_preview.png`：补偿与最终 SLM 灰度预览；
- `manifest.json`：全部参数、源文件和输出 SHA256；
- `README_EXPERIMENT.txt`：现场简要说明。

第一轮不加入像散和彗差。

## Y方向整幅相位平移扫描

当需要检查补偿相位与实际光束的相对位置时，使用
`export_zernike_compensated_y_sweep.py`。该脚本先合成`V2 + Z40/Z20`，再把
完整合成相位一起平移，最后加入闪耀光栅。因此球差补偿中心会随Y平移一起移动。

当前首轮固定`X=+5`，输出：

```text
Y = -10, -5, 0, +5, +10, +15, +20 computational pixels
```

这与单张基准导出器的“只平移V2、Zernike保持光轴居中”用途不同；扫描目录的
`manifest.json`会明确记录这一点，避免混用。
