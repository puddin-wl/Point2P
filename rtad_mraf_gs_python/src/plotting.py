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
    """Save target template and truncated-constraint plots."""
    outdir = Path(outdir)
    release_level = float(target.params.get("release_level", np.exp(-2.0)))

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(target.I_full, extent=_extent(target.x_um, target.y_um), origin="lower", cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, label="I_full")
    ax.contour(target.x_um, target.y_um, target.I_full, levels=[release_level], colors="cyan", linewidths=1.1)
    _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
    _plot_rect(ax, target.params["a50_um"], target.params["b50_um"], "-", "size50", "red")
    _plot_rect(ax, target.params["a1_um"], target.params["b1_um"], ":", "outer zero", "white")
    _set_zoom(ax, target.params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title(f"Full RTAD intensity template, release={release_level:.6g}")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(outdir / "target_full_intensity.png", dpi=dpi)
    fig.savefig(outdir / "target_intensity.png", dpi=dpi)
    plt.close(fig)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="0.6")
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(target.I_constraint, extent=_extent(target.x_um, target.y_um), origin="lower", cmap=cmap, vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, label="I_constraint")
    ax.contour(target.x_um, target.y_um, target.mask_signal.astype(float), levels=[0.5], colors="cyan", linewidths=1.0)
    ax.contour(target.x_um, target.y_um, target.mask_free.astype(float), levels=[0.5], colors="lime", linewidths=1.0)
    _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
    _set_zoom(ax, target.params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("Truncated RTAD constraint: cyan=signal, green=free")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(outdir / "target_constraint_intensity.png", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(target.A_full, extent=_extent(target.x_um, target.y_um), origin="lower", cmap="viridis", vmin=0, vmax=1)
    fig.colorbar(im, ax=ax, label="A_full")
    _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
    _plot_rect(ax, target.params["a50_um"], target.params["b50_um"], "-", "size50", "red")
    _plot_rect(ax, target.params["a1_um"], target.params["b1_um"], ":", "outer zero", "white")
    _set_zoom(ax, target.params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("Full RTAD amplitude template")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(outdir / "target_amplitude.png", dpi=dpi)
    plt.close(fig)

    plot_target_profiles(target, outdir, dpi=dpi)


def plot_masks(target: Any, outdir: str | Path, dpi: int = 150) -> None:
    """Save a discrete mask map."""
    outdir = Path(outdir)
    mask_show = np.zeros_like(target.I_full, dtype=np.uint8)
    mask_show[target.mask_bg_far] = 4
    mask_show[target.mask_free] = 3
    mask_show[target.mask_edge_lock] = 2
    mask_show[target.mask_flat] = 1
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    im = ax.imshow(mask_show, extent=_extent(target.x_um, target.y_um), origin="lower", cmap="tab10", vmin=1, vmax=4)
    cbar = fig.colorbar(im, ax=ax, ticks=[1, 2, 3, 4])
    cbar.ax.set_yticklabels(["flat", "edge_lock", "free", "far_bg"])
    _plot_rect(ax, target.params["a0_um"], target.params["b0_um"], "--", "flat core", "white")
    _plot_rect(ax, target.params["a50_um"], target.params["b50_um"], "-", "size50", "red")
    _plot_rect(ax, target.params["a1_um"], target.params["b1_um"], ":", "outer edge", "white")
    _set_zoom(ax, target.params)
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("Masks: 1=flat, 2=edge_lock, 3=free, 4=far_bg")
    fig.savefig(outdir / "masks.png", dpi=dpi)
    plt.close(fig)


def plot_target_profiles(target: Any, outdir: str | Path, dpi: int = 150) -> None:
    """Save full-template center profiles with release and geometry markers."""
    outdir = Path(outdir)
    release_level = float(target.params.get("release_level", np.exp(-2.0)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    panels = [
        (axes[0], target.profiles["x_um"], target.profiles["x_I"], "x", ["a0_um", "a50_um", "a1_um", "a2_um"]),
        (axes[1], target.profiles["y_um"], target.profiles["y_I"], "y", ["b0_um", "b50_um", "b1_um", "b2_um"]),
    ]
    for ax, coord, profile, label, keys in panels:
        ax.plot(coord, profile, linewidth=1.4, label=f"I_full {label} profile")
        for level, name in [(0.9, "90%"), (0.5, "50%"), (np.exp(-2.0), "13.5%"), (release_level, "release")]:
            ax.axhline(level, linestyle=":", linewidth=0.9, label=name)
        styles = ["--", "-", ":", "-."]
        colors = ["0.4", "red", "0.4", "tab:green"]
        for key, style, color in zip(keys, styles, colors):
            val = float(target.params[key])
            ax.axvline(-val, color=color, linestyle=style, linewidth=0.9)
            ax.axvline(val, color=color, linestyle=style, linewidth=0.9, label=key)
        ax.set_xlabel(f"{label} / um")
        ax.set_ylabel("I_full")
        ax.set_ylim(-0.05, 1.08)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(outdir / "target_profiles.png", dpi=dpi)
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
    center_x = float(target.params.get("center_x_um", 0.0))
    center_y = float(target.params.get("center_y_um", 0.0))
    ix0 = int(np.argmin(np.abs(target.x_um - center_x)))
    iy0 = int(np.argmin(np.abs(target.y_um - center_y)))
    constraint_x = np.asarray(target.I_constraint[iy0, :], dtype=float)
    constraint_y = np.asarray(target.I_constraint[:, ix0], dtype=float)
    release_level = float(target.params.get("release_level", np.exp(-2.0)))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, coord, full_prof, constraint_prof, iprof, rprof, half50, half0, half1, label in [
        (axes[0], target.profiles["x_um"], target.profiles["x_I"], constraint_x, init_p["x_profile"], ref_p["x_profile"], target.params["a50_um"], target.params["a0_um"], target.params["a1_um"], "x"),
        (axes[1], target.profiles["y_um"], target.profiles["y_I"], constraint_y, init_p["y_profile"], ref_p["y_profile"], target.params["b50_um"], target.params["b0_um"], target.params["b1_um"], "y"),
    ]:
        ax.plot(coord, full_prof, label="I_full template", linewidth=1.2, color="0.45", linestyle="--")
        ax.plot(coord, constraint_prof, label="I_constraint signal", linewidth=1.5, color="tab:blue")
        ax.plot(coord, iprof, label="initial", linewidth=1.0)
        ax.plot(coord, rprof, label="refined", linewidth=1.0)
        ax.axhline(0.5, color="0.4", linestyle=":", linewidth=1.0)
        ax.axhline(release_level, color="tab:blue", linestyle="-.", linewidth=0.9, label="release_level")
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
