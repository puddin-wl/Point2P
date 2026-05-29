"""Simulate beam mismatch: 5mm-designed phase illuminated with 3.5/4/5mm beams.
Save 2D intensity images + center profile comparison.
"""
import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, 'src')
from src.propagation import forward_fft, intensity, make_input_gaussian
from src.rtad_target import make_axis_um, make_rtad_rect_target
from src.metrics import width_at_level, center_profiles
from src.backend import get_backend

outdir = "mismatch_output"
import os; os.makedirs(outdir, exist_ok=True)

# Load 5mm designed phase
phase5 = np.load('artifacts/20260528-153911_rtad_mraf_gs_truncI0135/phase_refined.npy')
cfg = json.loads(open('artifacts/20260528-153911_rtad_mraf_gs_truncI0135/config_used.json').read())

N = phase5.shape[0]
focal_dx_um = cfg['grid']['focal_dx_um']
clear_aperture_m = cfg['physical']['clear_aperture_m']
dx_doe_m = 532e-9 * 200e-3 / (N * focal_dx_um * 1e-6)

x_um = make_axis_um(N, focal_dx_um)
y_um = make_axis_um(N, focal_dx_um)
target = make_rtad_rect_target(shape=phase5.shape, x_um=x_um, y_um=y_um,
    W50_um=330.0, H50_um=120.0,
    delta_x_um=15.0, delta_y_um=8.0,
    guard_x_um=20.0, guard_y_um=12.0,
    release_level=np.exp(-2.0), constraint_mode='truncated_rtad')
mask_flat = target.mask_flat

backend = get_backend(use_cupy=True, device_id=0, verbose=False)
xp = backend.xp
dtype = backend.float_dtype
phase_b = backend.to_backend(phase5, dtype=dtype)

results = {}
for beam_mm in [5.0, 4.0, 3.5]:
    amp = make_input_gaussian(phase5.shape, dx_doe_m,
        gaussian_1e2_diameter_m=beam_mm * 1e-3,
        clear_aperture_m=clear_aperture_m, xp=xp, dtype=dtype)
    field = amp * xp.exp(1j * phase_b).astype(backend.complex_dtype, copy=False)
    I = backend.to_numpy(intensity(forward_fft(field, xp), xp)).astype(np.float64)
    I_n = I / np.mean(I[mask_flat])
    flat_vals = I_n[mask_flat]
    rms = 100.0 * np.std(flat_vals) / np.mean(flat_vals)
    profs = center_profiles(I_n, x_um, y_um)
    s50x = width_at_level(x_um, profs['x_profile'], 0.5)
    s50y = width_at_level(y_um, profs['y_profile'], 0.5)
    results[beam_mm] = dict(I_n=I_n, I_raw=I, rms=rms, s50x=s50x, s50y=s50y,
                            x_profile=profs['x_profile'], y_profile=profs['y_profile'])

# ---- 2D intensity images (side by side) ----
a0 = target.params['a0_um']; b0 = target.params['b0_um']
a50 = target.params['a50_um']; b50 = target.params['b50_um']
span_x = a50 + 200; span_y = b50 + 120

fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
for ax, (beam_mm, label) in zip(axes, [(5.0, 'matched'), (4.0, 'mismatch'), (3.5, 'mismatch')]):
    r = results[beam_mm]
    I_n = r['I_n']
    vmax = max(2.0, float(np.nanpercentile(I_n, 99.5)))
    im = ax.imshow(I_n, extent=[x_um[0], x_um[-1], y_um[0], y_um[-1]],
                   origin='lower', cmap='magma', vmin=0, vmax=vmax)
    # target flat core rect
    xx = [-a0, a0, a0, -a0, -a0]; yy = [-b0, -b0, b0, b0, -b0]
    ax.plot(xx, yy, '--', color='white', linewidth=1.0, label='target flat')
    ax.plot([-a50, a50, a50, -a50, -a50], [-b50, -b50, b50, b50, -b50],
            '-', color='cyan', linewidth=1.0, label='target 50%')
    ax.set_xlim(-span_x, span_x); ax.set_ylim(-span_y, span_y)
    ax.set_xlabel('x / um'); ax.set_ylabel('y / um')
    ax.set_title(f'beam={beam_mm}mm ({label})\nRMS={r["rms"]:.1f}%, size50={r["s50x"]:.0f}x{r["s50y"]:.0f}um')
    ax.legend(fontsize=7, loc='upper right')
    plt.colorbar(im, ax=ax, label='I / mean(flat)', shrink=0.82)
fig.savefig(f'{outdir}/mismatch_2d_intensity.png', dpi=150)
plt.close(fig)
print(f"Saved: {outdir}/mismatch_2d_intensity.png")

# ---- Center profiles (X and Y) ----
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
colors = {5.0: 'tab:green', 4.0: 'tab:orange', 3.5: 'tab:red'}
for ax, axis_name, half50, half0, half1, prof_key in [
    (axes[0], 'x', a50, a0, target.params['a1_um'], 'x_profile'),
    (axes[1], 'y', b50, b0, target.params['b1_um'], 'y_profile')]:
    coord = x_um if axis_name == 'x' else y_um
    for beam_mm in [5.0, 4.0, 3.5]:
        r = results[beam_mm]
        prof = r[prof_key]
        ax.plot(coord, prof, color=colors[beam_mm], linewidth=1.2,
                label=f'beam={beam_mm}mm (RMS={r["rms"]:.1f}%)')
    ax.axhline(1.0, color='0.3', linestyle='--', linewidth=0.8)
    ax.axhline(0.9, color='0.5', linestyle=':', linewidth=0.7)
    ax.axhline(0.5, color='0.5', linestyle=':', linewidth=0.7)
    ax.axvline(-half0, color='0.5', linestyle='--', linewidth=0.7)
    ax.axvline(half0, color='0.5', linestyle='--', linewidth=0.7)
    ax.axvline(-half50, color='red', linestyle='-', linewidth=0.7)
    ax.axvline(half50, color='red', linestyle='-', linewidth=0.7)
    ax.set_xlim(-half1 - 80, half1 + 80)
    ax.set_ylim(0, 2.0)
    ax.set_xlabel(f'{axis_name} / um')
    ax.set_ylabel('I / mean(flat)')
    ax.set_title(f'{axis_name} center profile')
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
fig.savefig(f'{outdir}/mismatch_profiles.png', dpi=150)
plt.close(fig)
print(f"Saved: {outdir}/mismatch_profiles.png")

# Print summary
print()
print(f"{'Beam':>12s}  {'RMS%':>8s}  {'size50_x':>10s}  {'size50_y':>10s}")
print("-" * 50)
for b in [5.0, 4.0, 3.5]:
    r = results[b]
    print(f"beam={b:3.1f} {'matched' if b==5.0 else 'mismatch':>10s}  {r['rms']:7.2f}%  {r['s50x']:9.1f}  {r['s50y']:9.1f}")
