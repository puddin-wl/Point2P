"""Refine the V2 hollow reproduction with a two-dimensional input-amplitude dip."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from simulate_v2_hollow_scan import (
    ANALYSIS_ROOT,
    DEFAULT_BMDATA,
    DEFAULT_CASE_DIR,
    DEFAULT_HOLLOW_SUMMARY,
    DEFAULT_PERCENT_SUMMARY,
    Simulator,
    evaluate_image,
    load_experimental_target,
    plot_candidate_overview,
    plot_profile_comparison,
    write_scan_csv,
)


def build_cases() -> list[dict]:
    cases: list[dict] = []
    for depth in np.linspace(0.30, 0.90, 7):
        for sigma_x_mm in (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0):
            for sigma_y_mm in (1.2, 1.5, 1.8, 2.1, 2.4, 3.0):
                params = {
                    "amplitude_dip_family": "elliptical_input_dip",
                    "amplitude_dip_depth": float(depth),
                    "amplitude_dip_sigma_x_mm": float(sigma_x_mm),
                    "amplitude_dip_sigma_y_mm": float(sigma_y_mm),
                }
                cases.append(
                    {
                        "family": "elliptical_input_dip",
                        "name": (
                            f"elliptical_dip_d{depth:.3f}_"
                            f"sx{sigma_x_mm:.3f}_sy{sigma_y_mm:.3f}mm"
                        ),
                        "params": params,
                    }
                )
    return cases


def write_summary(path: Path, experiment: dict, rows: list[dict], elapsed: float) -> None:
    best = rows[0]
    target = experiment["metrics"]
    lines = [
        "# 二维椭圆输入凹陷细化扫描",
        "",
        f"- 案例数：{len(rows)}",
        f"- 耗时：{elapsed:.1f} s",
        "- V2 相位保持不变；正向传播仍复用项目 `forward_fft`。",
        "",
        "## 实验目标",
        "",
        f"- 中心窗口/核心均值：`{target['center_window_over_core_mean']:.4f}`",
        f"- 中间三分之一/左右两侧：`{target['middle_third_over_side_thirds']:.4f}`",
        f"- 核心 RMS：`{100.0 * target['core_rms_fraction']:.2f}%`",
        "",
        "## 最佳仿真",
        "",
        f"- 名称：`{best['name']}`",
        f"- 参数：`{json.dumps(best['params'], ensure_ascii=False)}`",
        f"- 分数：`{best['score']:.4f}`",
        f"- 中心窗口/核心均值："
        f"`{best['metrics']['center_window_over_core_mean']:.4f}`",
        f"- 中间三分之一/左右两侧："
        f"`{best['metrics']['middle_third_over_side_thirds']:.4f}`",
        f"- 核心 RMS：`{100.0 * best['metrics']['core_rms_fraction']:.2f}%`",
        "",
        "## 前十名",
        "",
        "| 排名 | 参数 | 分数 | 中心比 | 中部/两侧 | RMS |",
        "|---:|---|---:|---:|---:|---:|",
    ]
    for index, row in enumerate(rows[:10], start=1):
        lines.append(
            f"| {index} | `{json.dumps(row['params'], ensure_ascii=False)}` | "
            f"{row['score']:.3f} | "
            f"{row['metrics']['center_window_over_core_mean']:.3f} | "
            f"{row['metrics']['middle_third_over_side_thirds']:.3f} | "
            f"{100.0 * row['metrics']['core_rms_fraction']:.1f}% |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    experiment = load_experimental_target(
        args.bmdata, args.percent_summary, args.hollow_summary
    )
    simulator = Simulator(args.case_dir)
    cases = build_cases()
    rows: list[dict] = []
    started = time.perf_counter()
    for index, case in enumerate(cases, start=1):
        image = simulator.reconstruct(case["params"])
        _, metrics, score_details = evaluate_image(image, simulator, experiment)
        rows.append(
            {
                **case,
                "score": score_details["score"],
                "metrics": metrics,
                "score_details": score_details,
            }
        )
        if index == 1 or index % 40 == 0 or index == len(cases):
            print(
                f"[{index:3d}/{len(cases)}] {case['name']} "
                f"score={score_details['score']:.3f}",
                flush=True,
            )
    elapsed = time.perf_counter() - started
    rows.sort(key=lambda row: row["score"])

    (args.outdir / "refinement_results.json").write_text(
        json.dumps(
            {
                "case_dir": str(args.case_dir),
                "elapsed_seconds": elapsed,
                "candidate_count": len(rows),
                "experiment_metrics": experiment["metrics"],
                "candidates": rows,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    write_scan_csv(args.outdir / "refinement_results.csv", rows)
    best = rows[0]
    best_image = simulator.reconstruct(best["params"])
    np.save(args.outdir / "best_elliptical_dip_intensity.npy", best_image)
    (args.outdir / "best_elliptical_dip.json").write_text(
        json.dumps(best, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    plot_candidate_overview(
        args.outdir / "elliptical_dip_overview.png", experiment, simulator, rows
    )
    plot_profile_comparison(
        args.outdir / "elliptical_dip_profiles.png", experiment, simulator, rows
    )
    write_summary(args.outdir / "SUMMARY_REFINEMENT.md", experiment, rows, elapsed)
    print(
        json.dumps(
            {
                "candidate_count": len(rows),
                "elapsed_seconds": elapsed,
                "best": {
                    "name": best["name"],
                    "params": best["params"],
                    "score": best["score"],
                    "metrics": {
                        key: value
                        for key, value in best["metrics"].items()
                        if not key.startswith("unit_profile")
                    },
                },
            },
            ensure_ascii=False,
            indent=2,
        )
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--bmdata", type=Path, default=DEFAULT_BMDATA)
    parser.add_argument(
        "--percent-summary", type=Path, default=DEFAULT_PERCENT_SUMMARY
    )
    parser.add_argument("--hollow-summary", type=Path, default=DEFAULT_HOLLOW_SUMMARY)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ANALYSIS_ROOT / "results" / "02_elliptical_input_counterexample",
    )
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
