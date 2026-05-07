"""Plot and PDF report generation for real-world sweeps."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from .metrics import LEVEL_E2, TargetData, normalize_by_flat
from .utils import value_label


def _extent(target: TargetData) -> list[float]:
    return [float(target.x_um[0]), float(target.x_um[-1]), float(target.y_um[0]), float(target.y_um[-1])]


def _plot_rect(ax: Any, half_x: float, half_y: float, style: str, color: str, label: str) -> None:
    x = [-half_x, half_x, half_x, -half_x, -half_x]
    y = [-half_y, -half_y, half_y, half_y, -half_y]
    ax.plot(x, y, style, color=color, linewidth=1.0, label=label)


def _zoom_target(ax: Any, target: TargetData, margin_um: float = 65.0) -> None:
    params = target.params
    ax.set_xlim(-float(params.get("a2_um", 200.0)) - margin_um, float(params.get("a2_um", 200.0)) + margin_um)
    ax.set_ylim(-float(params.get("b2_um", 100.0)) - margin_um, float(params.get("b2_um", 100.0)) + margin_um)


def _add_target_overlays(ax: Any, target: TargetData) -> None:
    params = target.params
    _plot_rect(ax, float(params.get("a0_um", 150.0)), float(params.get("b0_um", 52.0)), "--", "white", "flat")
    _plot_rect(ax, float(params.get("a50_um", 165.0)), float(params.get("b50_um", 60.0)), "-", "cyan", "50%")
    _plot_rect(ax, float(params.get("a1_um", 180.0)), float(params.get("b1_um", 68.0)), ":", "white", "edge")


def normalized_intensity(intensity: np.ndarray, target: TargetData) -> np.ndarray:
    """Return display-normalized intensity."""
    return normalize_by_flat(intensity, target.mask_flat)[0]


def plot_intensity_image(
    intensity: np.ndarray,
    target: TargetData,
    outpath: str | Path,
    title: str,
    dpi: int = 150,
) -> None:
    """Save a normalized intensity image."""
    In = normalized_intensity(intensity, target)
    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    vmax = max(2.0, float(np.nanpercentile(In, 99.5)))
    image = ax.imshow(In, extent=_extent(target), origin="lower", cmap="magma", vmin=0, vmax=vmax)
    fig.colorbar(image, ax=ax, label="I / mean(flat)")
    _add_target_overlays(ax, target)
    _zoom_target(ax, target)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def _is_negative_defocus(label: str) -> bool:
    """Detect negative defocus from common label patterns."""
    import re
    m = re.search(r"df=([+-]?\d*\.?\d+)", label)
    if m:
        return float(m.group(1)) < 0
    return False


def make_profiles_figure(profile_results: list[dict[str, Any]], target: TargetData) -> plt.Figure:
    """Create a center-profile overlay figure.

    Positive defocus → solid line; negative defocus → dashed line.
    Legend is placed in a separate row below the plots.
    """
    n_labels = len(profile_results)
    legend_rows = max(1, (n_labels + 5) // 6)

    fig = plt.figure(figsize=(13.0, 4.5 + legend_rows * 0.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[5, legend_rows], hspace=0.35, wspace=0.28,
                          left=0.06, right=0.98, top=0.94, bottom=0.06)

    ax_x = fig.add_subplot(gs[0, 0])
    ax_y = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[1, :])
    ax_legend.axis("off")

    params = target.params
    for ax, axis_name in [(ax_x, "x"), (ax_y, "y")]:
        coord_key = f"{axis_name}_um"
        prof_key = f"{axis_name}_profile"
        half0 = float(params.get("a0_um" if axis_name == "x" else "b0_um"))
        half50 = float(params.get("a50_um" if axis_name == "x" else "b50_um"))
        half1 = float(params.get("a1_um" if axis_name == "x" else "b1_um"))
        for result in profile_results:
            label = result["label"]
            ls = "--" if _is_negative_defocus(label) else "-"
            profiles = result["details"]["profiles"]
            ax.plot(
                profiles[coord_key],
                profiles[prof_key],
                linewidth=1.0, alpha=0.82, linestyle=ls,
                label=label,
            )
        ax.axhline(0.9, linestyle=":", color="0.45", linewidth=0.9)
        ax.axhline(0.5, linestyle=":", color="0.25", linewidth=0.9)
        ax.axhline(LEVEL_E2, linestyle="-.", color="tab:blue", linewidth=0.9)
        ax.axvline(-half0, linestyle="--", color="0.55", linewidth=0.8)
        ax.axvline(half0, linestyle="--", color="0.55", linewidth=0.8)
        ax.axvline(-half50, linestyle="-", color="red", linewidth=0.8)
        ax.axvline(half50, linestyle="-", color="red", linewidth=0.8)
        ax.axvline(-half1, linestyle=":", color="0.55", linewidth=0.8)
        ax.axvline(half1, linestyle=":", color="0.55", linewidth=0.8)
        ax.set_xlim(-half1 - 90.0, half1 + 90.0)
        ax.set_ylim(0.0, 2.2)
        ax.set_xlabel(f"{axis_name} / um")
        ax.set_ylabel("I / mean(flat)")
        ax.set_title(f"{axis_name} center profile")
        ax.grid(True, alpha=0.25)

    handles, labels = ax_y.get_legend_handles_labels()
    ax_legend.legend(
        handles, labels,
        fontsize=4.5, ncol=6, loc="upper center",
        frameon=True, fancybox=True, framealpha=0.9,
    )
    return fig


def plot_profiles_overlay(
    profile_results: list[dict[str, Any]],
    target: TargetData,
    outpath: str | Path,
    dpi: int = 150,
) -> None:
    """Save center-profile overlays for a sweep."""
    fig = make_profiles_figure(profile_results, target)
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


TREND_METRICS = [
    ("size50_x_um", "size50 x / um"),
    ("size50_y_um", "size50 y / um"),
    ("rms_nonuniformity_percent", "RMS nonuniformity / %"),
    ("efficiency_e2_percent", "e^-2 efficiency / %"),
    ("aperture_throughput_percent", "aperture throughput / %"),
    ("center_offset_x_um", "center offset x / um"),
]


def _trend_x(rows: list[dict[str, Any]]) -> tuple[np.ndarray, list[str], str]:
    raw = [row.get("sweep_value_raw", row.get("sweep_value")) for row in rows]
    numeric: list[float] = []
    for idx, value in enumerate(raw):
        if isinstance(value, (int, float, np.integer, np.floating)):
            numeric.append(float(value))
        else:
            numeric.append(float(idx))
    labels = [value_label(value) for value in raw]
    x_label = str(rows[0].get("sweep_parameter", "case")) if rows else "case"
    return np.asarray(numeric, dtype=np.float64), labels, x_label


def make_metric_trends_figure(rows: list[dict[str, Any]]) -> plt.Figure:
    """Create metric trend plots for a sweep."""
    x, labels, x_label = _trend_x(rows)
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2), constrained_layout=True)
    for ax, (key, label) in zip(axes.ravel(), TREND_METRICS):
        y = np.asarray([row.get(key, np.nan) for row in rows], dtype=np.float64)
        ax.plot(x, y, marker="o", linewidth=1.2)
        ax.set_title(label)
        ax.set_xlabel(x_label)
        ax.grid(True, alpha=0.25)
        if len(labels) <= 9 and len(set(labels)) == len(labels):
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right")
    return fig


def plot_metric_trends(rows: list[dict[str, Any]], outpath: str | Path, dpi: int = 150) -> None:
    """Save metric trend plots."""
    fig = make_metric_trends_figure(rows)
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def _add_text_page(pdf: PdfPages, title: str, rows: list[dict[str, Any]], warnings: list[str]) -> None:
    fig = plt.figure(figsize=(8.5, 11.0))
    ax = fig.add_subplot(111)
    ax.axis("off")
    lines = [title, "", f"cases: {len(rows)}"]
    if rows:
        worst = min(
            rows,
            key=lambda row: row.get("efficiency_e2_percent", np.inf)
            if np.isfinite(row.get("efficiency_e2_percent", np.nan))
            else np.inf,
        )
        lines.extend(
            [
                f"worst by e^-2 efficiency: {worst.get('case_label')}",
                f"worst e^-2 efficiency: {worst.get('efficiency_e2_percent', np.nan):.6g} %",
                f"worst RMS nonuniformity: {worst.get('rms_nonuniformity_percent', np.nan):.6g} %",
            ]
        )
    if warnings:
        lines.extend(["", "warnings:"])
        lines.extend(f"- {warning}" for warning in warnings[:20])
    ax.text(0.06, 0.95, "\n".join(lines), va="top", ha="left", fontsize=11, family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def _add_intensity_page(pdf: PdfPages, intensity: np.ndarray, target: TargetData, title: str) -> None:
    In = normalized_intensity(intensity, target)
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
    vmax = max(2.0, float(np.nanpercentile(In, 99.5)))
    image = ax.imshow(In, extent=_extent(target), origin="lower", cmap="magma", vmin=0, vmax=vmax)
    fig.colorbar(image, ax=ax, label="I / mean(flat)")
    _add_target_overlays(ax, target)
    _zoom_target(ax, target)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


def write_pdf_report(
    outpath: str | Path,
    sweep_name: str,
    rows: list[dict[str, Any]],
    profile_results: list[dict[str, Any]],
    target: TargetData,
    results: list[dict[str, Any]],
    warnings: list[str],
) -> None:
    """Write a compact PDF report for one sweep."""
    with PdfPages(outpath) as pdf:
        _add_text_page(pdf, f"Real-world simulation: {sweep_name}", rows, warnings)
        trend_fig = make_metric_trends_figure(rows)
        pdf.savefig(trend_fig)
        plt.close(trend_fig)
        profile_fig = make_profiles_figure(profile_results, target)
        pdf.savefig(profile_fig)
        plt.close(profile_fig)
        for result in results:
            _add_intensity_page(pdf, result["intensity"], target, result["label"])

