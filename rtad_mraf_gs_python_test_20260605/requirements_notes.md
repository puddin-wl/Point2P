# Requirements Notes

This project is intended to run in the local `slmrtad` conda environment.

Expected packages:

- `slmsuite==0.4.1` for source-reference semantics and future integration
- `cupy-cuda12x` for GPU FFTs
- `numpy`
- `scipy`
- `h5py`
- `matplotlib`
- `opencv-python`
- `tqdm`

The code does not require importing slmsuite for the first implementation. It
uses a simplified local GS/MRAF/WGS loop based on the semantics read from
slmsuite 0.4.1.

Recommended checks:

```powershell
& 'D:\software\anaconda\envs\slmrtad\python.exe' -c "import cupy as cp; print(cp.cuda.runtime.getDeviceCount()); print(cp.cuda.runtime.getDeviceProperties(0)['name'])"
& 'D:\software\anaconda\envs\slmrtad\python.exe' -c "import slmsuite; print(slmsuite.__version__)"
```
