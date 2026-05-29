"""Batch-run all three beam diameters (5/6/7 mm) and convert to SLM format.

自动化三步骤, 遍历三种光束直径:
  1. 调用 make_phase0.py 生成 Romero-Dickey 初始相位
  2. 调用 run_rtad_mraf_gs_case.py 进行 WGS 精修 (method=wgs, flat_local, 200 轮)
  3. 调用 convert_to_slm.py 将精修相位转为 SLM 1024×1024 @ 17μm

Usage:
  python run_all_beam_diameters.py
"""

import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
MAKE_PHASE0 = THIS_DIR / "make_phase0.py"
RUN_REFINE = THIS_DIR / "run_rtad_mraf_gs_case.py"
CONVERT_TO_SLM = THIS_DIR / "convert_to_slm.py"

BEAM_DIAMETERS = [5, 6, 7]

ARTIFACTS_DIR = THIS_DIR / "artifacts"


def find_latest_artifacts_subdir() -> Path | None:
    """在 artifacts/ 下查找最新的时间戳子目录."""
    if not ARTIFACTS_DIR.exists():
        return None
    dirs = sorted(
        [d for d in ARTIFACTS_DIR.iterdir() if d.is_dir() and d.name != ".gitkeep"],
        reverse=True,
    )
    return dirs[0] if dirs else None


def main() -> int:
    """对 5/6/7mm 三种光束依次执行: 生成 phase0 → WGS 精修 → SLM 转换."""
    extra = [
        "--method", "wgs",
        "--wgs-strategy", "flat_local",
        "--iters", "200",
        "--wgs-feedback-exponent", "0.8",
        "--bg-factor", "0.9",
        "--no-swap-phase-xy",
    ]

    for d_mm in BEAM_DIAMETERS:
        print(f"\n{'='*60}")
        print(f"Beam diameter {d_mm} mm")
        print(f"{'='*60}")

        # Step 1: Generate phase0
        phase0_dir = THIS_DIR / f"make_phase0_output"
        phase0_dir.mkdir(parents=True, exist_ok=True)
        phase0_mat = phase0_dir / f"phase0_beam{d_mm}mm.mat"

        print(f"\n[1/3] Generating phase0 for beam={d_mm}mm...")
        result = subprocess.run([
            sys.executable, str(MAKE_PHASE0),
            "--beam", str(float(d_mm)),
            "--out", str(phase0_dir),
        ])
        if result.returncode != 0:
            print(f"ERROR: phase0 generation failed for beam {d_mm}mm", file=sys.stderr)
            return result.returncode

        # Rename phase0 to include beam diameter for clarity
        default_mat = phase0_dir / "phase0.mat"
        if default_mat.exists() and not phase0_mat.exists():
            default_mat.rename(phase0_mat)
            # Also rename config snapshot
            default_cfg = phase0_dir / "config_snapshot.mat"
            cfg_mat = phase0_dir / f"config_snapshot_beam{d_mm}mm.mat"
            if default_cfg.exists() and not cfg_mat.exists():
                default_cfg.rename(cfg_mat)

        # Step 2: WGS refinement
        print(f"\n[2/3] Running WGS refinement for beam={d_mm}mm...")
        result = subprocess.run([
            sys.executable, str(RUN_REFINE),
            "--phase-mat", str(phase0_mat),
            "--phase-var", "phase0_wrapped_rad",
            "--beam-diameter", str(float(d_mm)),
        ] + extra)
        if result.returncode != 0:
            print(f"ERROR: Refinement failed for beam {d_mm}mm", file=sys.stderr)
            return result.returncode

        # Step 3: Convert to SLM
        latest = find_latest_artifacts_subdir()
        if latest is None:
            print(f"ERROR: No artifacts found after refinement for beam {d_mm}mm", file=sys.stderr)
            return 1

        phase_refined = latest / "phase_refined.npy"
        slm_out = latest / "slm_output"

        print(f"\n[3/3] Converting to SLM format (1024x1024 @ 17um)...")
        result = subprocess.run([
            sys.executable, str(CONVERT_TO_SLM),
            str(phase_refined),
            "--out", str(slm_out),
            "--label", f"phase_slm_beam{d_mm}mm",
        ])
        if result.returncode != 0:
            print(f"ERROR: SLM conversion failed for beam {d_mm}mm", file=sys.stderr)
            return result.returncode

        print(f"\nBeam {d_mm}mm done. SLM output: {slm_out}")

    print(f"\n{'='*60}")
    print("All three beam diameters (5/6/7 mm) completed.")
    print("SLM BMP files are in each run's artifacts/<timestamp>/slm_output/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
