# Romero–Dickey 初始相位生成

本目录提供独立的 MATLAB 初始相位生成模块。它根据高斯入射光和矩形目标，使用 Romero–Dickey 点对点稳相方法生成 DOE 初始相位，并计算一次理想傅里叶透镜传播结果。

## 当前定位

- 默认物理配置是历史 5 mm、f=429 mm 基线；
- 当前 6.5 mm 实验基线可使用 Python `make_phase0.py` 或测试目录中已跟踪的冻结输入；
- 本目录输出位于 `artifacts/`，默认不进入 Git。

## 文件

- `default_initial_phase_config.m`：物理参数和采样配置；
- `generate_initial_phase.m`：相位生成核心函数；
- `run_initial_phase_generation.m`：单次运行入口；
- `run_initial_phase_generation_beam6p4x6p3.m` 等：椭圆光束实验入口。

## 运行

在 MATLAB 中：

```matlab
cd('E:\program\Point2P\initial_phase_generation')
run_initial_phase_generation
```

输出写入 `artifacts/<时间戳>/`，主要包括：

- `phase0.mat`：包裹/未包裹相位、坐标轴和配置；
- 初始相位图；
- 入射强度与理想焦面结果；
- 本次配置快照。

## 物理与坐标约定

初始相位按 X/Y 可分离映射计算。下游 Python 程序读取旧 MATLAB 文件时必须核对轴方向；新生成的 Python 相位通常使用 `--no-swap-phase-xy`。

修改光束直径、焦距或目标尺寸后，应生成新的命名基线，不要覆盖已有 artifact。
