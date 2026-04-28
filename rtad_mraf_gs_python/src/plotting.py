"""Plot generation for RTAD MRAF/GS runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .metrics import center_profiles, normalize_intensity


def _extent(x_um: np.ndarray, y_um: np.ndarray) -> list[float]:
    """Return imshow extent for x/y axes."""
    return [float(x_um[0]), float(x_um[-1]), float(y_um[0]), float(y_um[-1])]


def _plot_rect(ax: Any, a_um: float, b_um: float, style: str, label: str, color: str) -> None:
    """Overlay a centered rectangle."""
    xx = [-a_um, a_um, a_um, -a_um, -a_um]
    yy = [-b_um, -b_um, b_um, b_um, -b_um]
    ax.plot(xx, yy, style, color=color, linewidth=1.1, label=label)


def _set_zoom(ax: Any, params: dict[str, Any], margin_um: float = 60.0) -> None:
    """Zoom plots around the RTAD support and guard band."""
    ax.set_xlim(-params["a2_um"] - margin_um, params["a2_um"] + margin_um)
    ax.set_ylim(-params["b2_um"] - margin_um, params["b2_um"] + margin_um)


def plot_target(target: Any, outdir: str | Path, dpi: int = 150) -> None:
    """Save target intensity and amplitude plots."""
    outdir = Path(outdir)
    for name, arr, title in [
        ("target_intensity.png", target.I_target, "RTAD target intensity"),
        ("target_amplitude.png", target.A_target, "RTAD target amplitude"),
    ]:
        fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
        im = ax.imshow(arr, extent=_extent(target.x_um, target.y_um), origin="lower", cmap="viridis")
        fig.colorbar(im, ax=ax, label="normalized")
        _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
        _plot_rect(ax, target.params["a50_um"], target.params["b50_um"], "-", "size50", "red")
        _plot_rect(ax, target.params["a1_um"], target.params["b1_um"], ":", "outer edge", "white")
        _set_zoom(ax, target.params)
        ax.set_xlabel("x / um")
        ax.set_ylabel("y / um")
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=8)
        fig.savefig(outdir / name, dpi=dpi)
        plt.close(fig)


def plot_masks(target: Any, outdir: str | Path, dpi: int = 150) -> None:
    """Save a discrete mask map."""
    outdir = Path(outdir)
    mask_show = np.zeros_like(target.I_target, dtype=np.uint8)
    mask_show[target.mask_bg] = 0
    mask_show[target.mask_free] = 3
    mask_show[target.mask_edge] = 2
    mask_show[target.mask_flat] = 1
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(mask_show, extent=_extent(target.x_um, target.y_um), origin="lower", cmap="tab10", vmin=0, vmax=3)
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3])
    cbar.ax.set_yticklabels(["bg", "flat", "edge", "free"])
    _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
    _plot_rect(ax, target.params["a50_um"], target.params["b50_um"], "-", "size50", "red")
    _plot_rect(ax, target.params["a1_um"], target.params["b1_um"], ":", "outer edge", "white")
    _set_zoom(ax, target.params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("RTAD masks")
    fig.savefig(outdir / "masks.png", dpi=dpi)
    plt.close(fig)


def plot_intensity(
    intensity: np.ndarray,
    masks: dict[str, np.ndarray],
    x_um: np.ndarray,
    y_um: np.ndarray,
    params: dict[str, Any],
    outpath: str | Path,
    title: str,
    dpi: int = 150,
) -> None:
    """Save a normalized reconstruction intensity image."""
    In = normalize_intensity(intensity, masks["mask_flat"], mode="flat_mean")
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(In, extent=_extent(x_um, y_um), origin="lower", cmap="magma", vmin=0, vmax=max(2.0, float(np.nanpercentile(In, 99.5))))
    fig.colorbar(im, ax=ax, label="I / mean(flat)")
    _plot_rect(ax, params["a0_um"], params["b0_um"], "--", "flat core", "white")
    _plot_rect(ax, params["a50_um"], params["b50_um"], "-", "size50", "cyan")
    _plot_rect(ax, params["a1_um"], params["b1_um"], ":", "outer edge", "white")
    _set_zoom(ax, params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def plot_phase(phase: np.ndarray, outpath: str | Path, title: str, dpi: int = 150) -> None:
    """Save a wrapped phase image."""
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    im = ax.imshow(np.mod(phase, 2 * np.pi), origin="lower", cmap="twilight", vmin=0, vmax=2 * np.pi)
    fig.colorbar(im, ax=ax, label="rad")
    ax.set_title(title)
    ax.set_xlabel("col / x")
    ax.set_ylabel("row / y")
    fig.savefig(outpath, dpi=dpi)
    plt.close(fig)


def plot_profiles_compare(
    target: Any,
    initial_intensity: np.ndarray,
    refined_intensity: np.ndarray,
    masks: dict[str, np.ndarray],
    outdir: str | Path,
    dpi: int = 150,
) -> None:
    """Save x/y center profile comparison plots."""
    outdir = Path(outdir)
    init_n = normalize_intensity(initial_intensity, masks["mask_flat"], mode="flat_mean")
    ref_n = normalize_intensity(refined_intensity, masks["mask_flat"], mode="flat_mean")
    init_p = center_profiles(init_n, target.x_um, target.y_um)
    ref_p = center_profiles(ref_n, target.x_um, target.y_um)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, coord, tprof, iprof, rprof, half50, half0, half1, label in [
        (axes[0], target.profiles["x_um"], target.profiles["x_I"], init_p["x_profile"], ref_p["x_profile"], target.params["a50_um"], target.params["a0_um"], target.params["a1_um"], "x"),
        (axes[1], target.profiles["y_um"], target.profiles["y_I"], init_p["y_profile"], ref_p["y_profile"], target.params["b50_um"], target.params["b0_um"], target.params["b1_um"], "y"),
    ]:
        ax.plot(coord, tprof, label="target I", linewidth=1.4)
        ax.plot(coord, iprof, label="initial", linewidth=1.0)
        ax.plot(coord, rprof, label="refined", linewidth=1.0)
        ax.axhline(0.5, color="0.4", linestyle=":", linewidth=1.0)
        ax.axvline(-half0, color="0.5", linestyle="--", linewidth=0.9)
        ax.axvline(half0, color="0.5", linestyle="--", linewidth=0.9)
        ax.axvline(-half50, color="red", linestyle="-", linewidth=0.9)
        ax.axvline(half50, color="red", linestyle="-", linewidth=0.9)
        ax.axvline(-half1, color="0.5", linestyle=":", linewidth=0.9)
        ax.axvline(half1, color="0.5", linestyle=":", linewidth=0.9)
        ax.set_xlim(-half1 - 80, half1 + 80)
        ax.set_ylim(0, max(2.0, float(np.nanpercentile(rprof, 99.5))))
        ax.grid(True, alpha=0.25)
        ax.set_xlabel(f"{label} / um")
        ax.set_ylabel("I / mean(flat)")
        ax.set_title(f"{label} center profile")
        ax.legend(fontsize=8)
    fig.savefig(outdir / "center_profiles_compare.png", dpi=dpi)
    plt.close(fig)


def plot_convergence(metrics_history: list[dict[str, Any]], outdir: str | Path, dpi: int = 150) -> None:
    """Save convergence curves from the metrics history."""
    outdir = Path(outdir)
    if not metrics_history:
        return
    it = np.asarray([m["iteration"] for m in metrics_history], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    panels = [
        ("flat_rms", "flat RMS std/mean"),
        ("efficiency_support", "support efficiency"),
        ("background_fraction", "background fraction"),
        ("size50_x_um", "size50 x / um"),
    ]
    for ax, (key, title) in zip(axes.ravel(), panels):
        values = np.asarray([m.get(key, np.nan) for m in metrics_history], dtype=float)
        ax.plot(it, values, marker="o", markersize=3)
        ax.set_title(title)
        ax.set_xlabel("iteration")
        ax.grid(True, alpha=0.25)
    fig.savefig(outdir / "convergence_metrics.png", dpi=dpi)
    plt.close(fig)
