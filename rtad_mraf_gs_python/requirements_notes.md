# 环境依赖说明

推荐使用 `slmrtad` Conda 环境。仓库根目录提供：

- `requirements.txt`：CPU 和通用工具依赖；
- `requirements-gpu.txt`：CUDA 12 CuPy 与 slmsuite 增量依赖。

核心依赖包括 NumPy、SciPy、h5py、Matplotlib、OpenCV、Pillow 和 tqdm。GPU FFT 使用 `cupy-cuda12x`。`slmsuite==0.4.1` 主要用于语义参考和后续集成，当前局部 GS/MRAF/WGS 循环不依赖直接导入它。

```powershell
python -c "import cupy as cp; print(cp.cuda.runtime.getDeviceCount())"
python -c "import slmsuite; print(slmsuite.__version__)"
```
