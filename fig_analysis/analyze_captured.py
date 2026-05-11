"""Analyze a captured flat-top spot image.

Supports .mat (MATLAB), .bmp, .png, .tif, .jpg formats.

Usage:
  python analyze_captured.py <path/to/image> [--variable snapshot1] [--pixel-um 3.45]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.ndimage import binary_fill_holes, binary_closing, label, binary_dilation, uniform_filter1d
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _load_image(path: str, var_name: str = "snapshot1") -> np.ndarray:
    """Load image from .mat, .bmp, .png, .tif, or .jpg."""
    ext = Path(path).suffix.lower()
    if ext == ".mat":
        data = loadmat(path)
        if var_name in data:
            return np.asarray(data[var_name], dtype=np.float64)
        for k, v in data.items():
            if not k.startswith("__"):
                arr = np.asarray(v)
                if arr.ndim == 2 and arr.shape[0] > 100:
                    return arr.astype(np.float64)
        raise KeyError(f"No 2D array found in {path}")
    else:
        return np.array(Image.open(path), dtype=np.float64)


def analyze_captured(mat_path: str, var_name: str = "snapshot1", pixel_um: float = 3.45):
    img = _load_image(mat_path, var_name)
    print(f"File: {mat_path}")
    print(f"Shape: {img.shape}, dtype: {img.dtype}")
    print(f"Range: [{img.min():.0f}, {img.max():.0f}]")

    # Background
    corners = np.concatenate([
        img[:30, :30].ravel(), img[:30, -30:].ravel(),
        img[-30:, :30].ravel(), img[-30:, -30:].ravel(),
    ])
    bg_val = float(np.median(corners))
    bg_std = float(np.std(corners))
    print(f"Background: {bg_val:.1f} +/- {bg_std:.1f}")

    # Find bright region
    threshold = bg_val + 10 * bg_std
    bright = img > threshold
    bright = binary_closing(bright, structure=np.ones((3, 3)))
    bright = binary_fill_holes(bright)
    labeled, nf = label(bright)

    if nf == 0:
        print("No bright region found.")
        return

    comp_sizes = [(labeled == i).sum() for i in range(1, nf + 1)]
    main_id = int(np.argmax(comp_sizes)) + 1
    mask = labeled == main_id

    # Center of mass (more robust than bbox center)
    total_mass = float(img[mask].sum())
    y_idx, x_idx = np.where(mask)
    cx = float((x_idx * img[mask]).sum() / total_mass) if total_mass > 0 else float(x_idx.mean())
    cy = float((y_idx * img[mask]).sum() / total_mass) if total_mass > 0 else float(y_idx.mean())

    ymin, ymax = int(y_idx.min()), int(y_idx.max())
    xmin, xmax = int(x_idx.min()), int(x_idx.max())
    bbox_w, bbox_h = xmax - xmin + 1, ymax - ymin + 1

    flat_vals = img[mask]
    flat_mean = float(flat_vals.mean())
    flat_std = float(flat_vals.std())
    flat_min = float(flat_vals.min())
    flat_max = float(flat_vals.max())
    n_sat = int((flat_vals >= 255).sum())

    rms_pct = float(flat_std / flat_mean * 100) if flat_mean > 0 else float("nan")
    pv_pct = float((flat_max - flat_min) / flat_mean * 100) if flat_mean > 0 else float("nan")
    saturation_pct = float(100 * n_sat / len(flat_vals))

    print(f"\n--- Bright region ---")
    print(f"  Bbox: {bbox_w} x {bbox_h} px, center=({cx:.0f}, {cy:.0f})")
    print(f"  Pixels: {len(flat_vals)}")
    print(f"  Mean: {flat_mean:.1f}, std: {flat_std:.1f}")
    print(f"  RMS: {rms_pct:.2f}%  (WARNING: clipped if saturated)")
    print(f"  PV: {pv_pct:.1f}%")

    # Width analysis — gradient-edge-based method
    # 1. Extract center profiles
    # 2. Find gradient peaks → these bracket the flat-top region
    # 3. Flat level = median between gradient edges (robust against hotspots)
    # 4. Find 50%/90%/13.5% crossings relative to this flat level
    cxi, cyi = int(round(cx)), int(round(cy))
    prof_x = img[cyi, :].astype(np.float64)
    prof_y = img[:, cxi].astype(np.float64)
    prof_x_s = uniform_filter1d(prof_x, 7)
    prof_y_s = uniform_filter1d(prof_y, 7)

    x_bg = float(np.median(np.concatenate([prof_x[:50], prof_x[-50:]])))
    y_bg = float(np.median(np.concatenate([prof_y[:50], prof_y[-50:]])))

    # Gradient magnitude
    gx = np.abs(np.gradient(prof_x_s))
    gy = np.abs(np.gradient(prof_y_s))

    def _find_edge_pair(grad, center, search_radius=300):
        """Find the two dominant gradient peaks bracketing the flat-top."""
        left_grad = grad[max(0, center - search_radius):center]
        right_grad = grad[center:min(len(grad), center + search_radius)]
        thr = 0.15 * float(grad.max())

        def _peaks(arr, offset):
            out = []
            for i in range(2, len(arr) - 2):
                if arr[i] > thr and arr[i] >= arr[i-1] and arr[i] >= arr[i-2] \
                   and arr[i] > arr[i+1] and arr[i] > arr[i+2]:
                    out.append((offset + i, float(arr[i])))
            return out

        left_peaks = _peaks(left_grad, max(0, center - search_radius))
        right_peaks = _peaks(right_grad, center)

        left_edge = min((p[0] for p in sorted(left_peaks, key=lambda x: -x[1])[:3]), default=center - 50)
        right_edge = max((p[0] for p in sorted(right_peaks, key=lambda x: -x[1])[:3]), default=center + 50)
        return left_edge, right_edge

    lx, rx = _find_edge_pair(gx, cxi)
    ly, ry = _find_edge_pair(gy, cyi)

    # Flat level = median between gradient edges (exclude edge transition zones)
    margin = 8
    flat_x_vals = prof_x_s[max(0, lx + margin):min(len(prof_x_s), rx - margin)]
    flat_y_vals = prof_y_s[max(0, ly + margin):min(len(prof_y_s), ry - margin)]
    flat_level_x = float(np.median(flat_x_vals)) if len(flat_x_vals) > 0 else float(np.median(flat_vals))
    flat_level_y = float(np.median(flat_y_vals)) if len(flat_y_vals) > 0 else float(np.median(flat_vals))
    hotspot_x = float(flat_x_vals.max() / flat_level_x) if flat_level_x > 0 else 1.0
    hotspot_y = float(flat_y_vals.max() / flat_level_y) if flat_level_y > 0 else 1.0
    x_peak = float(prof_x_s.max())
    y_peak = float(prof_y_s.max())

    def _find_crossings(prof, bg, flat_level, center):
        results = {}
        for lbl, frac in [("90", 0.9), ("50", 0.5), ("13.5", 0.135)]:
            thr = bg + frac * (flat_level - bg)
            left = center
            for i in range(center, 0, -1):
                if prof[i] < thr:
                    left = i + (thr - prof[i]) / (prof[i+1] - prof[i]) if i+1 < len(prof) and abs(prof[i+1] - prof[i]) > 0.1 else float(i)
                    break
            right = center
            for i in range(center, len(prof)):
                if prof[i] < thr:
                    right = i - (prof[i] - thr) / (prof[i] - prof[i-1]) if i > 0 and abs(prof[i] - prof[i-1]) > 0.1 else float(i)
                    break
            results[lbl] = (float(left), float(right), float(right - left))
        return results

    x_cross = _find_crossings(prof_x_s, x_bg, flat_level_x, cxi)
    y_cross = _find_crossings(prof_y_s, y_bg, flat_level_y, cyi)

    w50_x = x_cross["50"][2]
    w50_y = y_cross["50"][2]
    w13_x = x_cross["13.5"][2]
    w13_y = y_cross["13.5"][2]
    w90_x = x_cross["90"][2]
    w90_y = y_cross["90"][2]
    trans_x = (w13_x - w90_x) / 2.0 if not np.isnan(w13_x) and not np.isnan(w90_x) else float("nan")
    trans_y = (w13_y - w90_y) / 2.0 if not np.isnan(w13_y) and not np.isnan(w90_y) else float("nan")
    aspect = w50_x / w50_y if w50_y and w50_y > 0 else float("nan")

    print(f"\n--- Gradient edges ---")
    print(f"  X edges: [{lx}, {rx}] (width={rx-lx}px)")
    print(f"  Y edges: [{ly}, {ry}] (width={ry-ly}px)")
    print(f"  Flat level: X={flat_level_x:.1f}, Y={flat_level_y:.1f}")
    print(f"  Hotspot ratio: X={hotspot_x:.2f}, Y={hotspot_y:.2f}  (>1.1 = non-uniform)")

    print(f"\n--- Widths ---")
    print(f"  size50:    {w50_x:.0f} x {w50_y:.0f} px  =  {w50_x*pixel_um:.0f} x {w50_y*pixel_um:.0f} um")
    print(f"  size13.5:  {w13_x:.0f} x {w13_y:.0f} px  =  {w13_x*pixel_um:.0f} x {w13_y*pixel_um:.0f} um")
    print(f"  size90:    {w90_x:.0f} x {w90_y:.0f} px  =  {w90_x*pixel_um:.0f} x {w90_y*pixel_um:.0f} um")
    print(f"  transition: {trans_x:.0f} x {trans_y:.0f} px")
    print(f"  aspect (50): {aspect:.2f}  (target: 2.75)")

    # Overexposure assessment
    print(f"\n--- Overexposure ---")
    print(f"  Saturated pixels in flat region: {n_sat}/{len(flat_vals)} ({saturation_pct:.1f}%)")
    print(f"  Flat mean / 255 = {flat_mean/255:.3f}")
    if saturation_pct > 1:
        print(f"  !! SEVERELY OVEREXPOSED !!")
        print(f"  {saturation_pct:.0f}% of the flat top is clipped at 255.")
        print(f"  True flat-top uniformity CANNOT be measured.")
        print(f"  Action: reduce laser power or camera exposure/gain.")
    elif saturation_pct > 0:
        print(f"  Mildly overexposed ({saturation_pct:.1f}% saturated).")
    else:
        print(f"  Exposure OK — no saturated pixels.")

    # ---- PLOT ----
    outdir = Path(__file__).resolve().parent / "output"
    outdir.mkdir(exist_ok=True)
    stem = Path(mat_path).stem

    # Zoom window for profiles (around center of mass)
    zoom_margin = max(bbox_w, bbox_h) + 60
    px0 = max(0, cxi - zoom_margin)
    px1 = min(img.shape[1], cxi + zoom_margin)
    py0 = max(0, cyi - zoom_margin)
    py1 = min(img.shape[0], cyi + zoom_margin)

    fig = plt.figure(figsize=(20, 12))

    # 1. Full image with ROI
    ax1 = fig.add_subplot(2, 3, (1, 2))
    ax1.imshow(img, cmap="gray", origin="upper", vmin=0, vmax=255, aspect="equal")
    ax1.add_patch(plt.Rectangle((px0, py0), px1-px0, py1-py0,
                                 fill=False, ec="lime", lw=2))
    ax1.add_patch(plt.Rectangle((xmin, ymin), bbox_w, bbox_h,
                                 fill=False, ec="red", lw=1.5, ls="--"))
    ax1.plot(cx, cy, "r+", ms=18, mew=2)
    ax1.set_title(f"Full image [{img.shape[1]}x{img.shape[0]}] | green=ROI | red=flat bbox\n"
                  f"bg={bg_val:.0f} | threshold={threshold:.0f}")
    plt.colorbar(ax1.images[0], ax=ax1, shrink=0.85)

    # 2. ROI zoom
    ax2 = fig.add_subplot(2, 3, 3)
    ax2.imshow(img, cmap="gray", origin="upper", vmin=0, vmax=255, aspect="equal")
    ax2.plot(cx, cy, "r+", ms=16, mew=2)
    ax2.add_patch(plt.Rectangle((xmin, ymin), bbox_w, bbox_h,
                                 fill=False, ec="red", lw=2))
    ax2.set_xlim(px0 - 10, px1 + 10)
    ax2.set_ylim(py1 + 10, py0 - 10)
    ax2.set_title(f"ROI zoom | bright region {bbox_w}x{bbox_h} px\n"
                  f"mean={flat_mean:.0f} saturation={saturation_pct:.1f}%")
    plt.colorbar(ax2.images[0], ax=ax2, shrink=0.85)

    # 3. X profile with gradient edges
    ax3 = fig.add_subplot(2, 3, 4)
    xx_full = np.arange(len(prof_x))
    ax3.plot(xx_full, prof_x, "k-", lw=0.5, alpha=0.4)
    ax3.plot(xx_full, prof_x_s, "b-", lw=1.5, label="smoothed")
    for lbl, (xl, xr, _) in x_cross.items():
        ls = "-" if lbl == "50" else "--"
        ax3.axvline(xl, color="r", ls=ls, lw=0.8)
        ax3.axvline(xr, color="r", ls=ls, lw=0.8)
    ax3.axvline(lx, color="orange", ls=":", lw=1.5, label=f"grad edges [{lx},{rx}]")
    ax3.axvline(rx, color="orange", ls=":", lw=1.5)
    ax3.axhline(x_bg, color="gray", ls=":", label=f"bg={x_bg:.0f}")
    ax3.axhline(flat_level_x, color="g", ls="--", lw=1.5, label=f"flat={flat_level_x:.0f}")
    ax3.axhline(255, color="r", ls=":", lw=0.6)
    ax3.set_xlim(px0, px1)
    ax3.set_xlabel("x / px"); ax3.set_ylabel("Intensity")
    ax3.set_title(f"X profile (y={cyi}) | size50={w50_x:.0f}px | hotspot={hotspot_x:.2f}")
    ax3.legend(fontsize=7); ax3.grid(alpha=0.3)

    # 4. Y profile with gradient edges
    ax4 = fig.add_subplot(2, 3, 5)
    yy_full = np.arange(len(prof_y))
    ax4.plot(yy_full, prof_y, "k-", lw=0.5, alpha=0.4)
    ax4.plot(yy_full, prof_y_s, "b-", lw=1.5, label="smoothed")
    for lbl, (yl, yr, _) in y_cross.items():
        ls = "-" if lbl == "50" else "--"
        ax4.axvline(yl, color="r", ls=ls, lw=0.8)
        ax4.axvline(yr, color="r", ls=ls, lw=0.8)
    ax4.axvline(ly, color="orange", ls=":", lw=1.5, label=f"grad edges [{ly},{ry}]")
    ax4.axvline(ry, color="orange", ls=":", lw=1.5)
    ax4.axhline(y_bg, color="gray", ls=":", label=f"bg={y_bg:.0f}")
    ax4.axhline(flat_level_y, color="g", ls="--", lw=1.5, label=f"flat={flat_level_y:.0f}")
    ax4.axhline(255, color="r", ls=":", lw=0.6)
    ax4.set_xlim(py0, py1)
    ax4.set_xlabel("y / px"); ax4.set_ylabel("Intensity")
    ax4.set_title(f"Y profile (x={cxi}) | size50={w50_y:.0f}px | hotspot={hotspot_y:.2f}")
    ax4.legend(fontsize=7); ax4.grid(alpha=0.3)

    # 5. 3D surface
    ax5 = fig.add_subplot(2, 3, 6, projection="3d")
    pd = 10
    y3 = np.arange(max(0, ymin - pd), min(img.shape[0], ymax + pd + 1))
    x3 = np.arange(max(0, xmin - pd), min(img.shape[1], xmax + pd + 1))
    XX, YY = np.meshgrid(x3, y3)
    ZZ = img[y3[0]:y3[-1] + 1, x3[0]:x3[-1] + 1]
    ax5.plot_surface(XX, YY, ZZ, cmap="hot", edgecolor="none", alpha=0.9, vmin=0, vmax=255)
    ax5.set_title("3D surface (hot colormap)" if saturation_pct > 0 else "3D surface")
    ax5.view_init(35, -50)

    fig.tight_layout()
    for ext in ["png", "jpg"]:
        outpath = outdir / f"{stem}_analysis.{ext}"
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
    plt.close(fig)

    # Save JSON summary
    summary = {
        "file": str(Path(mat_path).resolve()),
        "image_shape": list(img.shape),
        "background": bg_val,
        "background_std": bg_std,
        "gradient_edges_x_px": [lx, rx],
        "gradient_edges_y_px": [ly, ry],
        "flat_level_x": flat_level_x,
        "flat_level_y": flat_level_y,
        "hotspot_ratio_x": hotspot_x,
        "hotspot_ratio_y": hotspot_y,
        "flat_center_px": [cx, cy],
        "size50_px": [w50_x, w50_y],
        "size13p5_px": [w13_x, w13_y],
        "size90_px": [w90_x, w90_y],
        "transition_px": [trans_x, trans_y],
        "aspect_ratio": aspect,
        "rms_percent": rms_pct,
        "pv_percent": pv_pct,
        "saturated_percent": saturation_pct,
        "pixel_um": pixel_um,
        "size50_um": [w50_x * pixel_um, w50_y * pixel_um],
        "size13p5_um": [w13_x * pixel_um, w13_y * pixel_um],
    }
    json_path = outdir / f"{stem}_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved: {json_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze captured flat-top spot image.")
    parser.add_argument("mat_path", help="Path to .mat file")
    parser.add_argument("--variable", default="snapshot1", help="Variable name in .mat")
    parser.add_argument("--pixel-um", type=float, default=3.45, help="Camera pixel size in um")
    args = parser.parse_args()
    analyze_captured(args.mat_path, args.variable, args.pixel_um)
