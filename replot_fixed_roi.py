"""Revised analysis: lower threshold for wider ROI, show full image, detect internal dark regions"""
import sys, json
import numpy as np
from scipy.io import loadmat
from scipy.ndimage import uniform_filter1d, binary_fill_holes, binary_closing, label, binary_dilation
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mat_path = r"E:\program\Point2P\lab_test_f300mm\captured_2.mat"
data = loadmat(mat_path)

for k, v in data.items():
    if not k.startswith("__"):
        arr = np.asarray(v)
        if arr.ndim == 2 and arr.shape[0] > 100:
            img = arr.astype(np.float64)
            name = k
            break

print(f"Variable: {name}, shape: {img.shape}, range: [{img.min():.0f}, {img.max():.0f}]")

# Estimate background from corners
corners = np.concatenate([
    img[:30, :30].ravel(), img[:30, -30:].ravel(),
    img[-30:, :30].ravel(), img[-30:, -30:].ravel()
])
bg_val = np.median(corners)
bg_std = np.std(corners)
print(f"Background (corners): {bg_val:.1f} +/- {bg_std:.1f}")

# Use a lower threshold — 8*sigma above background
threshold = bg_val + 6 * bg_std
print(f"Detection threshold: {threshold:.1f}")

bright_mask = img > threshold
# Clean up but don't overdo it
bright_mask = binary_closing(bright_mask, structure=np.ones((3,3)))
bright_mask = binary_fill_holes(bright_mask)

labeled, nf = label(bright_mask)
print(f"Connected components: {nf}")

# Show component sizes
comp_sizes = [(labeled == i).sum() for i in range(1, nf + 1)]
for i, sz in enumerate(comp_sizes, 1):
    yc, xc = np.where(labeled == i)
    print(f"  Component {i}: {sz} px, center~({xc.mean():.0f},{yc.mean():.0f})")

# Take largest
main = np.argmax(comp_sizes) + 1
bright_mask = labeled == main
print(f"Kept component {main}: {comp_sizes[main-1]} px")

# Dilate the mask slightly to capture edges
bright_mask_dilated = binary_dilation(bright_mask, structure=np.ones((5,5)), iterations=2)

y_idx, x_idx = np.where(bright_mask_dilated)
ymin, ymax = y_idx.min(), y_idx.max()
xmin, xmax = x_idx.min(), x_idx.max()
cx, cy = (xmin+xmax)/2, (ymin+ymax)/2
bbox_w, bbox_h = xmax-xmin+1, ymax-ymin+1
print(f"Flat bbox (dilated): {bbox_w:.0f}x{bbox_h:.0f} px, center=({cx:.0f},{cy:.0f})")

# ---- Flat-top uniformity using the ORIGINAL (non-dilated) mask ----
flat_vals = img[bright_mask]
flat_mean = flat_vals.mean()
flat_std = flat_vals.std()
flat_min = flat_vals.min()
flat_max = flat_vals.max()
rms_pct = flat_std / flat_mean * 100
n_sat = int((flat_vals == 255).sum())

# Internal dark regions: pixels within flat mask that are below flat_mean - 1*std
dark_threshold = flat_mean - 1.0 * flat_std
internal_dark = bright_mask & (img < dark_threshold)
n_dark = internal_dark.sum()
dark_mean = img[internal_dark].mean() if n_dark > 0 else 0

print(f"\n--- Uniformity ---")
print(f"  Flat mean: {flat_mean:.1f}, std: {flat_std:.1f}")
print(f"  Flat range: [{flat_min:.0f}, {flat_max:.0f}]")
print(f"  RMS: {rms_pct:.2f}%")
print(f"  (max-min)/mean: {(flat_max-flat_min)/flat_mean*100:.1f}%")
print(f"  Saturated: {n_sat}/{len(flat_vals)} ({100*n_sat/len(flat_vals):.1f}%)")
print(f"  Internal dark pixels (< mean-1*std): {n_dark} ({100*n_dark/len(flat_vals):.1f}% of flat)")
if n_dark > 0:
    print(f"    dark region mean: {dark_mean:.1f}, darkest: {img[internal_dark].min():.0f}")

# ---- Marginals for width ----
pad = 80
mx0, mx1 = max(0, xmin-pad), min(img.shape[1], xmax+pad)
my0, my1 = max(0, ymin-pad), min(img.shape[0], ymax+pad)

x_marg = img[ymin:ymax+1, mx0:mx1+1].mean(axis=0)
y_marg = img[my0:my1+1, xmin:xmax+1].mean(axis=1)
x_marg_s = uniform_filter1d(x_marg.astype(np.float64), 5)
y_marg_s = uniform_filter1d(y_marg.astype(np.float64), 5)

x_bg_marg = np.median(np.concatenate([x_marg[:10], x_marg[-10:]]))
y_bg_marg = np.median(np.concatenate([y_marg[:10], y_marg[-10:]]))
x_pk, y_pk = x_marg_s.max(), y_marg_s.max()

def find_w(prof, bg, pk, frac):
    thr = bg + frac*(pk-bg)
    idx = np.where(prof > thr)[0]
    if len(idx)==0: return np.nan, np.nan, np.nan, np.nan
    return float(idx[-1]-idx[0]), float(idx[0]), float(idx[-1]), float(thr)

w50_x, *_ = find_w(x_marg_s, x_bg_marg, x_pk, 0.5)
w50_y, *_ = find_w(y_marg_s, y_bg_marg, y_pk, 0.5)
w13p5_x, *_ = find_w(x_marg_s, x_bg_marg, x_pk, 0.135)
w13p5_y, *_ = find_w(y_marg_s, y_bg_marg, y_pk, 0.135)
w90_x, *_ = find_w(x_marg_s, x_bg_marg, x_pk, 0.9)
w90_y, *_ = find_w(y_marg_s, y_bg_marg, y_pk, 0.9)
trans_x = (w13p5_x-w90_x)/2 if not np.isnan(w13p5_x) and not np.isnan(w90_x) else np.nan
trans_y = (w13p5_y-w90_y)/2 if not np.isnan(w13p5_y) and not np.isnan(w90_y) else np.nan

pixel = 3.45
aspect = w50_x/w50_y if w50_y > 0 else np.nan

print(f"\n--- Widths ---")
print(f"  size50:    {w50_x:.0f} x {w50_y:.0f} px  =  {w50_x*pixel:.0f} x {w50_y*pixel:.0f} um")
print(f"  size13p5:  {w13p5_x:.0f} x {w13p5_y:.0f} px  =  {w13p5_x*pixel:.0f} x {w13p5_y*pixel:.0f} um")
print(f"  size90:    {w90_x:.0f} x {w90_y:.0f} px  =  {w90_x*pixel:.0f} x {w90_y*pixel:.0f} um")
print(f"  transition: {trans_x:.0f} x {trans_y:.0f} px  =  {trans_x*pixel:.0f} x {trans_y*pixel:.0f} um")
print(f"  aspect: {aspect:.2f}  (target: {330/120:.2f})")

# ==== PLOT ====
fig = plt.figure(figsize=(20, 11))

# ---- 1. Full image (MATLAB imshow style) ----
ax1 = fig.add_subplot(2, 3, (1, 2))
ax1.imshow(img, cmap="gray", origin="upper", vmin=0, vmax=255, aspect="auto")
# Show ROI (generous)
roi_rect = plt.Rectangle((mx0, my0), mx1-mx0, my1-my0,
                          fill=False, ec="lime", lw=2, ls="-")
ax1.add_patch(roi_rect)
# Show flat bbox (dilated)
box_rect = plt.Rectangle((xmin, ymin), bbox_w-1, bbox_h-1,
                          fill=False, ec="red", lw=1.5, ls="--")
ax1.add_patch(box_rect)
ax1.plot(cx, cy, "r+", ms=18, mew=2)
ax1.set_title(f"Full image [{img.shape[1]}x{img.shape[0]}] | green=ROI | red=flat bbox\n"
              f"Background={bg_val:.0f}, threshold={threshold:.0f}")
ax1.set_xlabel("x / px")
ax1.set_ylabel("y / px")
plt.colorbar(ax1.images[0], ax=ax1, shrink=0.85)

# ---- 2. ROI zoom with internal dark regions marked ----
ax2 = fig.add_subplot(2, 3, 3)
# Overlay: base image + dark region marker
ax2.imshow(img, cmap="gray", origin="upper", vmin=0, vmax=255, aspect="auto")
# Show dark internal regions as blue overlay
dark_overlay = np.zeros((*img.shape, 4))
dark_overlay[..., 2] = 1.0  # blue channel
dark_overlay[..., 3] = internal_dark.astype(float) * 0.5  # alpha
ax2.imshow(dark_overlay, origin="upper", aspect="auto")

ax2.plot(cx, cy, "r+", ms=16, mew=2)
ax2.add_patch(plt.Rectangle((xmin, ymin), bbox_w-1, bbox_h-1,
                             fill=False, ec="red", lw=1.5, ls="--"))
ax2.set_xlim(mx0-10, mx1+10)
ax2.set_ylim(my1+10, my0-10)
ax2.set_title(f"ROI zoom | blue=internal dark regions\n"
              f"({n_dark} px < mean-1*std, darkest={internal_dark.min() if n_dark>0 else 'N/A'})")
ax2.set_xlabel("x / px")
ax2.set_ylabel("y / px")
plt.colorbar(ax2.images[0], ax=ax2, shrink=0.85)

# ---- 3. X center profile ----
prof_x = img[int(round(cy)), :]
prof_x_s = uniform_filter1d(prof_x.astype(np.float64), 3)
ax3 = fig.add_subplot(2, 3, 4)
xx = np.arange(len(prof_x))
ax3.plot(xx, prof_x, "k-", lw=0.7, alpha=0.6)
ax3.plot(xx, prof_x_s, "b-", lw=1.5, label="Smoothed")
ax3.axhline(flat_mean, color="g", ls="--", lw=1.2, label=f"Flat mean={flat_mean:.0f}")
ax3.axhline(flat_mean - flat_std, color="orange", ls=":", lw=1, label=f"mean-1*std={flat_mean-flat_std:.0f}")
ax3.axvspan(xmin, xmax, alpha=0.1, color="green")
ax3.set_xlim(mx0-5, mx1+5)
ax3.set_ylim(0, 260)
ax3.set_xlabel("x / px")
ax3.set_ylabel("Intensity")
ax3.set_title(f"X center profile (y={cy:.0f})")
ax3.legend(fontsize=8)
ax3.grid(alpha=0.3)

# ---- 4. Y center profile ----
prof_y = img[:, int(round(cx))]
prof_y_s = uniform_filter1d(prof_y.astype(np.float64), 3)
ax4 = fig.add_subplot(2, 3, 5)
yy = np.arange(len(prof_y))
ax4.plot(yy, prof_y, "k-", lw=0.7, alpha=0.6)
ax4.plot(yy, prof_y_s, "b-", lw=1.5, label="Smoothed")
ax4.axhline(flat_mean, color="g", ls="--", lw=1.2)
ax4.axhline(flat_mean - flat_std, color="orange", ls=":", lw=1)
ax4.axvspan(ymin, ymax, alpha=0.1, color="green")
ax4.set_xlim(my0-5, my1+5)
ax4.set_ylim(0, 260)
ax4.set_xlabel("y / px")
ax4.set_ylabel("Intensity")
ax4.set_title(f"Y center profile (x={cx:.0f})")
ax4.legend(fontsize=8)
ax4.grid(alpha=0.3)

# ---- 6. 3D surface ----
ax6 = fig.add_subplot(2, 3, 6, projection="3d")
pd_surf = 10
y3 = np.arange(max(0,ymin-pd_surf), min(img.shape[0],ymax+pd_surf+1))
x3 = np.arange(max(0,xmin-pd_surf), min(img.shape[1],xmax+pd_surf+1))
XX, YY = np.meshgrid(x3, y3)
ZZ = img[y3[0]:y3[-1]+1, x3[0]:x3[-1]+1]
ax6.plot_surface(XX, YY, ZZ, cmap="gray", edgecolor="none", alpha=0.9, vmin=0, vmax=255)
ax6.set_title("3D surface (gray, [0,255])")
ax6.set_xlabel("x/px"); ax6.set_ylabel("y/px")
ax6.view_init(35, -50)

fig.tight_layout()
outdir = Path(r"E:\program\Point2P\fig_analysis")
outdir.mkdir(exist_ok=True)
for fmt, ext in [("png", "png"), ("jpg", "jpg")]:
    path = outdir / f"captured_2_fixed_roi.{ext}"
    if ext == "jpg":
        from PIL import Image as PILImage
        tmp = outdir / "_tmp.png"
        fig.savefig(tmp, dpi=150)
        PILImage.open(tmp).convert('RGB').save(path, quality=92)
        tmp.unlink()
    else:
        fig.savefig(path, dpi=150)
    print(f"Saved: {path}")
plt.close(fig)

print("\nDone.")
