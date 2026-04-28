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


def plot_profiles_compare_initial_mraf_wgs(
    target: Any,
    initial_intensity: np.ndarray,
    mraf_intensity: np.ndarray,
    wgs_intensity: np.ndarray,
    masks: dict[str, np.ndarray],
    outdir: str | Path,
    dpi: int = 150,
) -> None:
    """Save center profiles for initial, MRAF-end, and final WGS outputs."""
    outdir = Path(outdir)
    init_n = normalize_intensity(initial_intensity, masks["mask_flat"], mode="flat_mean")
    mraf_n = normalize_intensity(mraf_intensity, masks["mask_flat"], mode="flat_mean")
    wgs_n = normalize_intensity(wgs_intensity, masks["mask_flat"], mode="flat_mean")
    init_p = center_profiles(init_n, target.x_um, target.y_um)
    mraf_p = center_profiles(mraf_n, target.x_um, target.y_um)
    wgs_p = center_profiles(wgs_n, target.x_um, target.y_um)
    release_level = float(target.params.get("release_level", np.exp(-2.0)))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, coord, iprof, mprof, wprof, half50, half0, half1, label in [
        (
            axes[0],
            target.profiles["x_um"],
            init_p["x_profile"],
            mraf_p["x_profile"],
            wgs_p["x_profile"],
            target.params["a50_um"],
            target.params["a0_um"],
            target.params["a1_um"],
            "x",
        ),
        (
            axes[1],
            target.profiles["y_um"],
            init_p["y_profile"],
            mraf_p["y_profile"],
            wgs_p["y_profile"],
            target.params["b50_um"],
            target.params["b0_um"],
            target.params["b1_um"],
            "y",
        ),
    ]:
        ax.plot(coord, iprof, label="initial", linewidth=0.9, alpha=0.8)
        ax.plot(coord, mprof, label="after MRAF", linewidth=1.1)
        ax.plot(coord, wprof, label="after WGS", linewidth=1.1)
        ax.axhline(0.9, color="0.45", linestyle=":", linewidth=0.8, label="90%")
        ax.axhline(0.5, color="0.35", linestyle=":", linewidth=0.9, label="50%")
        ax.axhline(np.exp(-2.0), color="0.45", linestyle=":", linewidth=0.8, label="13.5%")
        ax.axhline(release_level, color="tab:blue", linestyle="-.", linewidth=0.8, label="release")
        ax.axvline(-half0, color="0.5", linestyle="--", linewidth=0.8)
        ax.axvline(half0, color="0.5", linestyle="--", linewidth=0.8)
        ax.axvline(-half50, color="red", linestyle="-", linewidth=0.8)
        ax.axvline(half50, color="red", linestyle="-", linewidth=0.8)
        ax.axvline(-half1, color="0.5", linestyle=":", linewidth=0.8)
        ax.axvline(half1, color="0.5", linestyle=":", linewidth=0.8)
        ax.set_xlim(-half1 - 80, half1 + 80)
        ax.set_ylim(0, max(2.0, float(np.nanpercentile(wprof, 99.5))))
        ax.grid(True, alpha=0.25)
        ax.set_xlabel(f"{label} / um")
        ax.set_ylabel("I / mean(flat)")
        ax.set_title(f"{label} center profile")
        ax.legend(fontsize=8)
    fig.savefig(outdir / "center_profiles_compare_initial_mraf_wgs.png", dpi=dpi)
    plt.close(fig)


def plot_wgs_weights(
    weights: np.ndarray,
    mask_flat: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
    outdir: str | Path,
    dpi: int = 150,
) -> None:
    """Save final WGS flat-core weights as an image and histogram."""
    outdir = Path(outdir)
    weights = np.asarray(weights, dtype=np.float32)
    mask_flat = np.asarray(mask_flat, dtype=bool)
    display = np.full_like(weights, np.nan, dtype=np.float32)
    display[mask_flat] = weights[mask_flat]
    vals = weights[mask_flat]

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="0.65")
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    if vals.size and np.all(np.isfinite(vals)):
        vmin = float(np.nanpercentile(vals, 1))
        vmax = float(np.nanpercentile(vals, 99))
        if abs(vmax - vmin) < 1e-6:
            vmin, vmax = float(np.nanmin(vals)), float(np.nanmax(vals) + 1e-6)
    else:
        vmin, vmax = 0.5, 2.0
    im = ax.imshow(display, extent=_extent(x_um, y_um), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, label="WGS weight")
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_title("Final WGS weights in mask_flat")
    if vals.size:
        x_vals = np.asarray(x_um)
        y_vals = np.asarray(y_um)
        yy, xx = np.where(mask_flat)
        ax.set_xlim(float(x_vals[max(0, xx.min() - 20)]), float(x_vals[min(x_vals.size - 1, xx.max() + 20)]))
        ax.set_ylim(float(y_vals[max(0, yy.min() - 20)]), float(y_vals[min(y_vals.size - 1, yy.max() + 20)]))
    fig.savefig(outdir / "wgs_weights_final.png", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.hist(vals[np.isfinite(vals)], bins=60, color="tab:blue", alpha=0.85)
    ax.set_xlabel("WGS weight in mask_flat")
    ax.set_ylabel("count")
    ax.set_title("Final WGS weight distribution")
    ax.grid(True, alpha=0.25)
    fig.savefig(outdir / "wgs_weight_hist.png", dpi=dpi)
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
