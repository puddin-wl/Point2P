"""一键 Pipeline: Phase0 → WGS 精修 → SLM 转换 (f=429mm, SLM 1024×1024 @ 17μm).

默认使用 MATLAB baseline phase0:
  ../initial_phase_generation/artifacts/20260428-141942/phase0.mat

用法:
  python run_pipeline.py                     # 全自动: WGS 精修 → SLM 转换
  python run_pipeline.py --skip-refine       # 只做 SLM 转换 (需要已有 phase_refined.npy)
  python run_pipeline.py --phase0 my.mat     # 自定义 phase0
  python run_pipeline.py --beam 5.5          # 覆盖光束直径 (默认 5mm)
  python run_pipeline.py --outdir my_output  # 指定输出目录
"""

import argparse
import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent

PHASE0_DEFAULT = "../initial_phase_generation/artifacts/20260428-141942/phase0.mat"
PHASE0_VAR = "phase0_wrapped_rad"

RUN_REFINE = THIS_DIR / "run_rtad_mraf_gs_case.py"
CONVERT_SLM = THIS_DIR / "convert_to_slm.py"


def latest_artifacts_dir(root: Path) -> Path | None:
    """返回 artifacts/ 下最新的 (按创建时间) 子目录."""
    artifacts = root / "artifacts"
    if not artifacts.is_dir():
        return None
    dirs = [d for d in artifacts.iterdir() if d.is_dir()]
    if not dirs:
        return None
    dirs.sort(key=lambda d: d.stat().st_ctime, reverse=True)
    return dirs[0]


def run(cmd: list[str], step_name: str) -> int:
    """打印命令并执行, 失败时退出."""
    print(f"\n{'='*60}")
    print(f"  {step_name}")
    print(f"{'='*60}")
    print(f"  {' '.join(str(c) for c in cmd)}")
    print(f"{'='*60}\n")
    ret = subprocess.run(cmd, cwd=str(THIS_DIR)).returncode
    if ret != 0:
        print(f"\nERROR: {step_name} 失败 (exit code {ret})")
        sys.exit(ret)
    return ret


def main():
    parser = argparse.ArgumentParser(
        description="一键 Pipeline: Phase0 → WGS 精修 → SLM 转换 (f=429mm)"
    )
    parser.add_argument("--skip-refine", action="store_true",
                        help="跳过 WGS 精修, 只做 SLM 转换")
    parser.add_argument("--phase0", default=None,
                        help=f"phase0.mat 路径 (默认: {PHASE0_DEFAULT})")
    parser.add_argument("--beam", type=float, default=None,
                        help="光束 1/e^2 直径 (mm), 默认 5.0")
    parser.add_argument("--outdir", default=None,
                        help="指定输出目录 (默认自动生成时间戳目录)")
    args = parser.parse_args()

    phase0_path = Path(args.phase0) if args.phase0 else Path(PHASE0_DEFAULT)

    if not phase0_path.is_absolute():
        phase0_path = (THIS_DIR / phase0_path).resolve()

    if not phase0_path.exists():
        print(f"ERROR: phase0 文件不存在: {phase0_path}")
        sys.exit(1)

    # ── Stage 2: WGS 精修 ──
    if args.skip_refine:
        print("跳过 WGS 精修 (--skip-refine)")
    else:
        cmd = [
            sys.executable, str(RUN_REFINE),
            "--phase-mat", str(phase0_path),
            "--phase-var", PHASE0_VAR,
        ]
        if args.beam is not None:
            cmd.extend(["--beam-diameter", str(args.beam)])
        if args.outdir is not None:
            cmd.extend(["--outdir", args.outdir])
        run(cmd, "Stage 2: WGS 精修 (f=429mm)")

    # ── Stage 3: SLM 转换 ──
    artifacts_root = THIS_DIR
    if args.outdir:
        latest = Path(args.outdir)
    else:
        latest = latest_artifacts_dir(artifacts_root)

    if latest is None:
        print("ERROR: 找不到 artifacts 目录, 请先运行 WGS 精修")
        sys.exit(1)

    phase_npy = latest / "phase_refined.npy"
    if not phase_npy.exists():
        print(f"ERROR: 找不到 phase_refined.npy: {phase_npy}")
        sys.exit(1)

    slm_out = latest / "slm_output"
    cmd = [
        sys.executable, str(CONVERT_SLM),
        str(phase_npy),
        "--out", str(slm_out),
        "--f", "0.429",
    ]
    run(cmd, "Stage 3: SLM 转换 (1024x1024 @ 17um)")

    print(f"\n{'='*60}")
    print(f"  Pipeline 完成!")
    print(f"  SLM 相位输出: {slm_out}")
    print(f"  文件: phase_slm_1024x1024.npy / .mat / .bmp")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
