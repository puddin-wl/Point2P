"""Run lightweight Python diagnostics for an existing RTAD MRAF/GS case.

诊断内容包括:
  - size50 / size13p5 / size90 (中心剖面线性插值求 crossing)
  - RMS 不均匀性 (仅用 mask_flat 内像素)
  - e⁻² 效率 (13.5% crossing 矩形内的能量比例)
  - 导数旁瓣检测 (从中心向外, dI/dr 应在 90% crossing 后单调下降)
  - 中心偏移量

调用方式: python run_diagnostics_case.py <case_dir>
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from src.diagnostics import run_case_diagnostics
from src.utils import format_metrics


def parse_args() -> argparse.Namespace:
    """解析命令行参数: case_dir 和可选的 --output-dir."""
    parser = argparse.ArgumentParser(description="Run lightweight RTAD Python diagnostics.")
    parser.add_argument("case_dir", help="Refinement case directory containing reconstruction_refined.npy and target files.")
    parser.add_argument("--output-dir", default=None, help="Diagnostics output directory. Defaults to case_dir/diagnostics_python.")
    return parser.parse_args()


def main() -> int:
    """入口: 加载 case → 运行诊断 → 打印指标."""
    args = parse_args()
    metrics, outdir = run_case_diagnostics(args.case_dir, output_dir=args.output_dir)
    print(format_metrics(metrics))
    print(f"Diagnostics output directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
