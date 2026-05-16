"""Batch-run all three beam diameters (5/6/7 mm) for f=200mm lab test.

Usage:
  python run_all_beam_diameters.py --phase-mat-5mm <path> --phase-mat-6mm <path> --phase-mat-7mm <path> [-- extra args]

  All extra arguments after -- are forwarded to run_rtad_mraf_gs_case.py.
  Example:
    python run_all_beam_diameters.py \
      --phase-mat-5mm matlab/artifacts/20260511-120000/phase0.mat \
      --phase-mat-6mm matlab/artifacts/20260511-120001/phase0.mat \
      --phase-mat-7mm matlab/artifacts/20260511-120002/phase0.mat \
      -- --iters 200
"""

import argparse
import subprocess
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
MAIN_SCRIPT = THIS_DIR / "run_rtad_mraf_gs_case.py"

BEAM_DIAMETERS = [5, 6, 7]


def main() -> int:
    parser = argparse.ArgumentParser(description="Batch f=200mm refinement for 5/6/7 mm beams.")
    parser.add_argument("--phase-mat-5mm", required=True, help="Path to phase0.mat for 5mm beam.")
    parser.add_argument("--phase-mat-6mm", required=True, help="Path to phase0.mat for 6mm beam.")
    parser.add_argument("--phase-mat-7mm", required=True, help="Path to phase0.mat for 7mm beam.")
    args, extra = parser.parse_known_args()

    phase_paths = {
        5: args.phase_mat_5mm,
        6: args.phase_mat_6mm,
        7: args.phase_mat_7mm,
    }

    if extra and extra[0] == "--":
        extra = extra[1:]

    for d_mm in BEAM_DIAMETERS:
        print(f"\n{'='*60}")
        print(f"Running beam diameter {d_mm} mm...")
        print(f"{'='*60}")
        cmd = [
            sys.executable, str(MAIN_SCRIPT),
            "--phase-mat", phase_paths[d_mm],
            "--beam-diameter", str(float(d_mm)),
        ] + extra
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"ERROR: Beam {d_mm}mm failed with code {result.returncode}", file=sys.stderr)
            return result.returncode
    print("\nAll three beam diameters completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
