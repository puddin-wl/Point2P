# 环境依赖说明

本目录与通用 RTAD/WGS 管线使用同一个 `slmrtad` 环境。统一依赖文件位于仓库根目录：

```powershell
python -m pip install -r ../requirements.txt
python -m pip install -r ../requirements-gpu.txt
```

运行 Zernike 导出还需要 SciPy、Pillow 和 h5py；这些已包含在基础依赖中。
