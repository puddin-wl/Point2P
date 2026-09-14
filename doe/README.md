# 15 mm 方形 DOE 转换

本目录把 2026-07-15 V2 的纯 WGS 相位转换成 15×15 mm 方形、多级刻蚀 DOE。

## 本项目保留的文件

- `prepare_v2_doe_15mm.py`：从 2048² WGS 相位生成 1024² DOE 相位、掩膜、预览、清单和 GDS；
- `kalyout_doemake_v2_15mm.m`：适配 15 mm V2 的 MATLAB 版本；
- `v2_wgs_15mm_doe_20260730/README.md`：制造参数和层定义；
- `v2_wgs_15mm_doe_20260730/metadata.json`：输出哈希与统计；
- `v2_wgs_15mm_doe_20260730/doe_v2_wgs_15mm_square_4mask.gds`：当前四掩膜版图交付文件。

用户提供的原始 `kalyout_doemake.m` 作为本机参考保留，不直接进入 Git。DOE 转换
产生的其他 `.npy`、`.mat` 和预览图均可重新生成，默认不进入 Git。唯一例外是
作为 GDS 原始输入的 V2 `phase_refined.npy`，该文件已通过 Git LFS 永久保存。

## Python 运行

默认源文件是 7 月 15 日 V2 `phase_refined.npy`，对应 `config_used.json` 也已保存。
全新克隆后先获取 LFS 文件：

```powershell
git lfs pull
```

随后可直接运行默认命令，或显式传入路径：

```powershell
python prepare_v2_doe_15mm.py `
  --source <phase_refined.npy> `
  --config <config_used.json> `
  --output v2_wgs_15mm_doe_20260730
```

具体参数以 `python prepare_v2_doe_15mm.py --help` 为准。

原始相位 SHA256 为
`cf513dd1135234c1afef33c2a347b88028171f68517ba882ad98d420f75eef61`。
相同相位和配置会重现完全相同的 1024² 相位、量化级别与 GDS 几何；MAT/GDS
文件头包含生成时间，因此重新生成的文件级 SHA256 不保证相同。

## MATLAB GDSII 工具箱

`kalyout_doemake_v2_15mm.m` 依赖第三方 GDSII Toolbox。下载源码和本机生成的 `.mexw64` 不直接纳入 Point2P；应按上游说明安装并把工具箱加入 MATLAB 路径。

当前本机测试已确认：

- MinGW-w64 编译成功；
- 23 个 MEX 组件可用；
- GDS 写入和回读通过。

## 加工前必须确认

- 衬底材料及 532 nm 折射率；
- 77.5 nm 基础刻蚀深度；
- 当前 0–14 级是否改成严格 0–15 级；
- 正胶/负胶、掩膜明暗极性和刻蚀顺序；
- Layer 100 只是 15 mm 参考边界，不是刻蚀层。
