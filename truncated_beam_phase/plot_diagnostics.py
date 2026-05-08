"""Quick diagnostics: focal-plane intensity + centre profiles for refined phase."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = THIS_DIR.parent
sys.path.insert(0, str(PROJECT_DIR))
sys.path.insert(0, str(PROJECT_DIR / "real_world_simulation"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.metrics import generate_target, normalize_by_flat, LEVEL_E2

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
ARTIFACT_DIR = THIS_DIR / "artifacts" / "20260508-185317_wgs_refined"

phase_np = np.load(ARTIFACT_DIR / "phase_refined.npy").astype(np.float32)
recon_np = np.load(ARTIFACT_DIR / "reconstruction_refined.npy").astype(np.float32)
amp_np = np.load(ARTIFACT_DIR / "input_amplitude.npy").astype(np.float32)

N = phase_np.shape[0]
focal_dx_um = 2.5

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------
target = generate_target(
    shape=(N, N), focal_dx_um=focal_dx_um, focal_dy_um=focal_dx_um,
    target_config={
        "W50_um": 330.0, "H50_um": 120.0,
        "delta_x_um": 15.0, "delta_y_um": 8.0,
        "guard_x_um": 20.0, "guard_y_um": 12.0,
        "constraint_mode": "truncated_rtad",
        "release_level": float(np.exp(-2.0)),
    },
)

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
I = recon_np
I_norm, flat_mean = normalize_by_flat(I, target.mask_flat)
flat_vals = I_norm[target.mask_flat]
rms = 100.0 * float(np.std(flat_vals)) / float(np.mean(flat_vals))

# Centre profiles
x_um = target.x_um
y_um = target.y_um
cy = N // 2
cx = N // 2
prof_x = I_norm[cy, :]
prof_y = I_norm[:, cx]

def crossing(axis, profile, level):
    mid = len(profile) // 2
    if profile[mid] < level:
        return np.nan
    # Left crossing (negative side)
    left = np.nan
    for i in range(mid, 0, -1):
        v0, v1 = profile[i-1], profile[i]
        if (v0 <= level <= v1) or (v1 <= level <= v0):
            t = (level - v0) / (v1 - v0) if abs(v1 - v0) > 1e-30 else 0.0
            left = axis[i-1] + t * (axis[i] - axis[i-1])
            break
    # Right crossing (positive side)
    right = np.nan
    for i in range(mid, len(profile) - 1):
        v0, v1 = profile[i], profile[i+1]
        if (v0 >= level >= v1) or (v1 >= level >= v0):
            t = (level - v0) / (v1 - v0) if abs(v1 - v0) > 1e-30 else 0.0
            right = axis[i] + t * (axis[i+1] - axis[i])
            break
    if np.isfinite(left) and np.isfinite(right):
        return float(right - left)
    return np.nan

size50_x = crossing(x_um, prof_x, 0.5)
size50_y = crossing(y_um, prof_y, 0.5)
size13x = crossing(x_um, prof_x, LEVEL_E2)
size13y = crossing(y_um, prof_y, LEVEL_E2)

# Efficiency
total = float(np.sum(I))
XX, YY = np.meshgrid(x_um, y_um)
if np.isfinite(size13x) and np.isfinite(size13y):
    x_roi = np.abs(XX) <= size13x / 2.0
    y_roi = np.abs(YY) <= size13y / 2.0
    eff = 100.0 * float(np.sum(I[y_roi & x_roi])) / total if total > 0 else np.nan
else:
    eff = np.nan

print(f"RMS nonuniformity: {rms:.4f}%")
print(f"size50_x: {size50_x:.1f} um  (target: 330)")
print(f"size50_y: {size50_y:.1f} um  (target: 120)")
print(f"size13.5_x: {size13x:.1f} um")
print(f"size13.5_y: {size13y:.1f} um")
print(f"e^-2 efficiency: {eff:.2f}%")

# ---------------------------------------------------------------------------
# Plot: intensity
# ---------------------------------------------------------------------------
a2 = float(target.params.get("a2_um", 200))
b2 = float(target.params.get("b2_um", 100))
a0 = float(target.params.get("a0_um", 150))
b0 = float(target.params.get("b0_um", 52))
a50 = float(target.params.get("a50_um", 165))
b50 = float(target.params.get("b50_um", 60))

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)

# Intensity image
ax = axes[0]
vmax = max(2.0, float(np.nanpercentile(I_norm, 99.5)))
ext = [float(x_um[0]), float(x_um[-1]), float(y_um[0]), float(y_um[-1])]
im = ax.imshow(I_norm, extent=ext, origin="lower", cmap="magma", vmin=0, vmax=vmax)
fig.colorbar(im, ax=ax, label="I / mean(flat)")

# Target overlays
for hx, hy, ls, color, lbl in [
    (a0, b0, "--", "white", "flat"),
    (a50, b50, "-", "cyan", "50%"),
]:
    xr = [-hx, hx, hx, -hx, -hx]
    yr = [-hy, -hy, hy, hy, -hy]
    ax.plot(xr, yr, ls, color=color, linewidth=1.0, label=lbl)

margin = 65
ax.set_xlim(-a2 - margin, a2 + margin)
ax.set_ylim(-b2 - margin, b2 + margin)
ax.set_xlabel("x / um")
ax.set_ylabel("y / um")
ax.set_title(f"Focal-plane intensity (RMS={rms:.2f}%)")
ax.legend(loc="upper right", fontsize=8)

# Centre profiles
ax_x = axes[1]

color = "tab:blue"
ax_x.plot(x_um, prof_x, linewidth=1.2, color="tab:blue", label="x profile")
ax_x.plot(y_um, prof_y, linewidth=1.2, color="tab:red", label="y profile")

# Target boundaries
for hx, ls, c in [(a0, "--", "0.55"), (a50, "-", "red")]:
    ax_x.axvline(-hx, linestyle=ls, color=c, linewidth=0.8)
    ax_x.axvline(hx, linestyle=ls, color=c, linewidth=0.8)
for hy, ls, c in [(b0, "--", "0.55"), (b50, "-", "red")]:
    ax_x.axvline(-hy, linestyle=ls, color="0.55", linewidth=0.6, alpha=0.5)
    ax_x.axvline(hy, linestyle=ls, color="0.55", linewidth=0.6, alpha=0.5)

# Level lines
ax_x.axhline(0.9, linestyle=":", color="0.45", linewidth=0.9)
ax_x.axhline(0.5, linestyle=":", color="0.25", linewidth=0.9)
ax_x.axhline(LEVEL_E2, linestyle="-.", color="tab:blue", linewidth=0.9)

ax_x.set_xlim(-a2 - 90, a2 + 90)
ax_x.set_ylim(0, 2.2)
ax_x.set_xlabel("um")
ax_x.set_ylabel("I / mean(flat)")
ax_x.set_title("Centre profiles (x + y)")
ax_x.legend(fontsize=8)
ax_x.grid(True, alpha=0.25)

outpath = ARTIFACT_DIR / "diagnostics_summary.png"
fig.savefig(outpath, dpi=150)
plt.close(fig)
print(f"\nSaved: {outpath}")
