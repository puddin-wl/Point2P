# 真实光路误差模拟使用说明

本目录固定 DOE 相位，扫描实验中可能出现的误差，评估平顶光斑对光束和光路参数的敏感性。它是容差分析工具，不负责重新优化相位。

## 可扫描项目

- 离焦；
- 入射波前曲率/发散角；
- DOE 面光束 X/Y 偏移；
- 入射光斑直径；
- 椭圆度；
- 指向引起的焦面偏移；
- 孔径和组合误差。

可用扫描名以 `python run_real_world_sweep.py --help` 为准。

## 命令

默认温和范围：

```powershell
python run_real_world_sweep.py --sweep all
```

单项扫描：

```powershell
python run_real_world_sweep.py --sweep defocus
python run_real_world_sweep.py --sweep beam_offset_x
python run_real_world_sweep.py --sweep divergence
```

宽范围压力测试：

```powershell
python run_real_world_sweep.py --sweep all --profile stress
```

快速 CPU 检查：

```powershell
python run_real_world_sweep.py --sweep all --smoke-size 512 --no-pdf
```

## 输出

输出位于 `artifacts/<时间戳>/`，每个扫描项包含参数表、指标表、趋势图和汇总。读取结果时优先比较：

- RMS 非均匀性；
- size50 X/Y；
- e⁻² 效率；
- 中心偏移；
- 肩峰和旁瓣；
- 误差是否造成对称或非对称形变。

该模块给出的容差只对指定固定相位和当前模型成立，不能直接当作加工公差或实验验收标准。
