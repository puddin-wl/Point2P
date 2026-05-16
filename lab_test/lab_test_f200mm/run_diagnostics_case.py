"""Run lightweight Python diagnostics for an existing RTAD MRAF/GS case."""

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
    """Parse command-line options."""
    parser = argparse.ArgumentParser(description="Run lightweight RTAD Python diagnostics.")
    parser.add_argument("case_dir", help="Refinement case directory containing reconstruction_refined.npy and target files.")
    parser.add_argument("--output-dir", default=None, help="Diagnostics output directory. Defaults to case_dir/diagnostics_python.")
    return parser.parse_args()


def main() -> int:
    """Run diagnostics and print the metrics."""
    args = parse_args()
    metrics, outdir = run_case_diagnostics(args.case_dir, output_dir=args.output_dir)
    print(format_metrics(metrics))
    print(f"Diagnostics output directory: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
