# Point2P：矩形平顶光斑与 DOE 设计项目

Point2P 用于把 532 nm 高斯光束整形成矩形平顶光斑。项目包含 Romero–Dickey 初始相位、RTAD 目标、WGS 精修、SLM 相位导出、实验图像分析、像差补偿和 15 mm 方形 DOE 制造版图。

> 当前状态以 [`docs/PROJECT_STATUS.md`](docs/PROJECT_STATUS.md) 为唯一入口。历史文档中的参数只代表当时实验，不应直接当作最新结论。

## 当前结论

| 状态 | 结论 |
|---|---|
| 仿真参考基线 | 2026-06-05，6.5 mm 入射光，目标 330×120 μm，RMS 非均匀性 1.403%，size50 为 329.56×119.82 μm |
| 当前实验基线 | 2026-07-15 V2，设计目标补偿为 W50=375.02 μm、H50=137 μm，安装补偿 X+5/Y+5 |
| 现场加载文件 | `phase_20260715_targetExpX330_v2_shiftX+5_Y+5.bmp`，包含 X/Y 闪耀光栅 |
| 纯 WGS 相位 | 同一 V2 目录中的 `phase_refined.npy`，不含安装平移和闪耀光栅 |
| 中心空洞根因 | 后续排查确认是 PBS 在极高功率下损坏；原球差/离焦拟合仅作历史排查记录 |
| Zernike 补偿 | 2026-07-23 的 Z40/Z20 空洞补偿方案已停用，不再作为后续实验依据 |
| 15 mm DOE | 2026-07-30 已生成四掩膜 GDS；脚本和制造说明纳入项目，数值中间文件按生成物管理 |

## 目录说明

```text
Point2P/
├── docs/                              # 项目状态、文档索引和数据管理约定
├── rtad_mraf_gs_python/               # f=429 mm 通用仿真管线
├── rtad_mraf_gs_python_test_20260605/ # 6.5 mm 冻结基线、V2实验和补偿工具
├── initial_phase_generation/          # MATLAB Romero–Dickey 初始相位
├── lab_test/                          # f=100/200/300 mm 与实验反馈 WGS 验证
├── real_test/                         # 实验数据分析记录
├── hollow_flattop_analysis_20260716/  # 中心空洞历史分析（已由 PBS 损伤结论取代）
├── real_world_simulation/             # 入射和光路误差容差扫描
├── fig_analysis/                      # 实验光斑图像分析
├── result_diagnostics/                # MATLAB 焦面诊断
├── truncated_beam_phase/              # 截断光束数值相位
├── doe/                               # 15 mm DOE 转换脚本与制造交付说明
├── presentation/                      # 汇报材料
└── target/                            # MATLAB 目标函数实现
```

`rtad_mraf_gs_python` 是通用算法目录；`rtad_mraf_gs_python_test_20260605` 是当前实验演进目录。为保留日志中的绝对路径，暂不改名或合并二者。

## 环境

推荐使用现有 `slmrtad` Conda 环境。CPU 基础依赖见 `requirements.txt`，CUDA 12 GPU 增量依赖见 `requirements-gpu.txt`。

```powershell
conda activate slmrtad
python -m pip install -r requirements.txt
# 使用 NVIDIA CUDA 12 时再安装：
python -m pip install -r requirements-gpu.txt
```

## 快速验证

### 运行当前 6.5 mm 冻结基线

仓库中已保存可复现的 6.5 mm 初始相位。以下命令不依赖被忽略的本地历史输出：

```powershell
cd rtad_mraf_gs_python_test_20260605
python run_rtad_mraf_gs_case.py
```

运行时间较长或没有 GPU 时，可先做 CPU 冒烟测试：

```powershell
python run_rtad_mraf_gs_case.py --iters 2 --smoke-shape 256 --no-cupy
```

### 从头生成初始相位

```powershell
cd rtad_mraf_gs_python
python make_phase0.py --beam 6.5 --out artifacts/phase0_beam6p5mm
python run_rtad_mraf_gs_case.py `
  --phase-mat artifacts/phase0_beam6p5mm/phase0.mat `
  --phase-var phase0_wrapped_rad `
  --beam-diameter 6.5 `
  --method wgs --wgs-strategy flat_local --iters 200 `
  --wgs-feedback-exponent 0.8 --wgs-weight-min 0.5 `
  --wgs-weight-max 1.5 --bg-factor 0.9 --no-swap-phase-xy
```

### 运行仓库测试

```powershell
python -m unittest discover -s rtad_mraf_gs_python_test_20260605 -p "test_*.py"
```

## 关键约定

- `phase_refined.npy` 是 WGS 精修后的连续相位，不等于最终 SLM 加载图。
- SLM 导出必须对 `exp(iφ)` 的实部和虚部分别插值，不能直接插值包裹相位。
- 设计宽度、实验测量宽度和供应商口径必须同时注明阈值定义。
- `artifacts/`、原始相机数据和大体积中间数组默认不进入 Git；稳定结果应保存清单、参数和 SHA256。
- 历史日志只追加、不覆盖；新的“当前结论”写入项目状态文档。

## 进一步阅读

- [当前状态与基线](docs/PROJECT_STATUS.md)
- [文档导航](docs/README.md)
- [数据与生成物管理](docs/DATA_POLICY.md)
- [稳定基线校验清单](baselines/README.md)
- [算法与架构说明](Point2P_Complete_Summary.md)
- [V2 实验日志索引](rtad_mraf_gs_python_test_20260605/log/README.md)
- [中心空洞根因更正](hollow_flattop_analysis_20260716/docs/ROOT_CAUSE_UPDATE_20260912.md)
