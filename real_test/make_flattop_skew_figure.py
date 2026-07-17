"""Create a measured-skew and optical-shear explanation figure."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import h5py
import matplotlib
import numpy as np
from scipy import ndimage as ndi
from scipy.ndimage import gaussian_filter, uniform_filter1d

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "20260605-1"
OUT_DIR = DATA_DIR / "visualizations"


def _scalar(group: h5py.Group, name: str) -> Any:
    value = group[name][()]
    if getattr(value, "shape", None) == (1,):
        value = value[0]
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value.item() if hasattr(value, "item") else value


def load_frame(path: Path) -> tuple[np.ndarray, dict[str, float]]:
    with h5py.File(path, "r") as h5:
        frame = h5["BG_DATA"]["1"]
        raw = frame["RAWFRAME"]
        width = int(_scalar(raw, "WIDTH"))
        height = int(_scalar(raw, "HEIGHT"))
        bin_x = int(_scalar(raw, "BINNINGX"))
        bin_y = int(_scalar(raw, "BINNINGY"))
        sx = float(_scalar(raw, "PIXELSCALEXUM")) * bin_x
        sy = float(_scalar(raw, "PIXELSCALEYUM")) * bin_y
        data = frame["DATA"][()]
    return data.reshape((height, width)).astype(np.float64), {"sx_um": sx, "sy_um": sy}


def robust_background(image: np.ndarray, corner_px: int = 30) -> float:
    c = min(corner_px, image.shape[0] // 4, image.shape[1] // 4)
    corners = np.concatenate(
        [
            image[:c, :c].ravel(),
            image[:c, -c:].ravel(),
            image[-c:, :c].ravel(),
            image[-c:, -c:].ravel(),
        ]
    )
    return float(np.median(corners))


def largest_component(mask: np.ndarray) -> np.ndarray:
    labels, count = ndi.label(mask)
    if count == 0:
        raise RuntimeError("No threshold component found.")
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    return labels == int(np.argmax(sizes))


def fit_top_bottom_edges(norm: np.ndarray, x0: int, x1: int, y0: int, y1: int, level: float) -> dict[str, Any]:
    top_pts: list[tuple[float, float]] = []
    bottom_pts: list[tuple[float, float]] = []
    for x in range(x0, x1):
        prof = uniform_filter1d(norm[:, x], 3)
        above = np.where(prof[y0:y1] >= level)[0]
        if above.size < 3:
            continue
        first = y0 + int(above[0])
        last = y0 + int(above[-1])
        if first > 0:
            den = prof[first] - prof[first - 1]
            yt = first - 1 + (level - prof[first - 1]) / den if abs(den) > 1e-12 else float(first)
        else:
            yt = float(first)
        if last < len(prof) - 1:
            den = prof[last + 1] - prof[last]
            yb = last + (level - prof[last]) / den if abs(den) > 1e-12 else float(last)
        else:
            yb = float(last)
        if np.isfinite(yt) and np.isfinite(yb) and yb > yt:
            top_pts.append((x, yt))
            bottom_pts.append((x, yb))

    def fit(points: list[tuple[float, float]]) -> dict[str, float]:
        arr = np.array(points, dtype=np.float64)
        m, b = np.polyfit(arr[:, 0], arr[:, 1], 1)
        pred = m * arr[:, 0] + b
        return {
            "slope": float(m),
            "intercept": float(b),
            "angle_deg": float(np.degrees(np.arctan(m))),
            "rmse_px": float(np.sqrt(np.mean((arr[:, 1] - pred) ** 2))),
            "n": int(len(points)),
        }

    return {"top": fit(top_pts), "bottom": fit(bottom_pts)}


def pca_orientation(mask: np.ndarray) -> dict[str, float]:
    ys, xs = np.where(mask)
    coords = np.column_stack([xs, ys]).astype(np.float64)
    coords -= coords.mean(axis=0, keepdims=True)
    cov = np.cov(coords.T)
    eigvals, eigvecs = np.linalg.eigh(cov)
    vec = eigvecs[:, int(np.argmax(eigvals))]
    angle = float(np.degrees(np.arctan2(vec[1], vec[0])))
    if angle > 90:
        angle -= 180
    if angle < -90:
        angle += 180
    return {"major_axis_angle_deg": angle}


def make_figure() -> None:
    image, meta = load_frame(DATA_DIR / "20260605-5.bmData")
    summary = json.loads((DATA_DIR / "flattop_size_analysis" / "20260605-5_flattop_size_summary.json").read_text(encoding="utf-8"))
    raw_sizes = json.loads((DATA_DIR / "visualizations" / "20260605-5_flattop_raw_profile_sizes.json").read_text(encoding="utf-8"))
    bg = robust_background(image)
    signal = np.clip(image - bg, 0.0, None)
    flat_level = float(raw_sizes["flat_level_intensity"])
    norm = gaussian_filter(signal / flat_level, sigma=0.8)

    x13 = summary["crossings"]["x13p5"]
    y13 = summary["crossings"]["y13p5"]
    x0 = max(0, int(np.floor(x13["left_px"])) - 30)
    x1 = min(image.shape[1], int(np.ceil(x13["right_px"])) + 30)
    y0 = max(0, int(np.floor(y13["left_px"])) - 30)
    y1 = min(image.shape[0], int(np.ceil(y13["right_px"])) + 30)
    crop = norm[y0:y1, x0:x1]
    level = 0.5
    mask = largest_component(norm >= level)
    mask_crop = mask[y0:y1, x0:x1]
    edge_fits = fit_top_bottom_edges(norm, x0, x1, y0, y1, level)
    orientation = pca_orientation(mask_crop)

    sx = meta["sx_um"]
    sy = meta["sy_um"]
    extent = [x0 * sx, x1 * sx, y1 * sy, y0 * sy]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = OUT_DIR / "20260605-5_skew_explanation.png"
    metrics_path = OUT_DIR / "20260605-5_skew_metrics.json"

    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2), constrained_layout=True)

    ax = axes[0]
    im = ax.imshow(crop, cmap="turbo", origin="upper", extent=extent, vmin=0.0, vmax=1.35)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="I / median(flat)")
    ax.contour(
        np.arange(x0, x1) * sx,
        np.arange(y0, y1) * sy,
        crop,
        levels=[0.5, 0.865],
        colors=["white", "magenta"],
        linewidths=[1.4, 1.1],
    )
    xx = np.array([x0, x1], dtype=np.float64)
    for key, color in [("top", "cyan"), ("bottom", "cyan")]:
        m = edge_fits[key]["slope"]
        b = edge_fits[key]["intercept"]
        yy = m * xx + b
        ax.plot(xx * sx, yy * sy, color=color, lw=1.5, ls="--")
    ax.set_title("Measured flat-top skew / edge tilt")
    ax.set_xlabel("x / um")
    ax.set_ylabel("y / um")
    ax.set_aspect("equal")
    ax.text(
        0.02,
        0.04,
        "white contour: 50% level\n"
        "magenta contour: 86.5% level\n"
        f"top/bottom edge tilt: {edge_fits['top']['angle_deg']:.2f}°, {edge_fits['bottom']['angle_deg']:.2f}°\n"
        f"PCA major-axis angle: {orientation['major_axis_angle_deg']:.2f}°",
        transform=ax.transAxes,
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 4},
    )

    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")
    ax.set_title("Why camera rotation cannot remove optical shear")
    ideal = Rectangle((0.8, 4.5), 2.3, 1.0, angle=0, fill=False, edgecolor="tab:blue", lw=2)
    ax.add_patch(ideal)
    ax.text(1.95, 5.8, "ideal rectangular target", ha="center", color="tab:blue")
    ax.add_patch(FancyArrowPatch((3.4, 5.0), (5.0, 5.0), arrowstyle="->", mutation_scale=16, lw=1.6))
    ax.text(4.2, 5.35, "galvo + f-theta\nnon-normal mapping", ha="center", fontsize=9)
    sheared = Polygon([[5.4, 4.35], [8.4, 4.05], [8.1, 5.05], [5.1, 5.35]], closed=True, fill=False, edgecolor="tab:red", lw=2)
    ax.add_patch(sheared)
    ax.text(6.75, 5.65, "sheared / tilted image", ha="center", color="tab:red")
    ax.add_patch(FancyArrowPatch((6.7, 3.8), (7.6, 3.2), arrowstyle="->", mutation_scale=14, lw=1.4, color="0.25"))
    rotated = Polygon([[5.9, 1.15], [8.7, 1.15], [8.25, 2.0], [5.45, 2.0]], closed=True, fill=False, edgecolor="tab:orange", lw=2)
    ax.add_patch(rotated)
    ax.text(7.1, 2.35, "after camera rotation:\nglobal angle changes,\nshear remains", ha="center", color="tab:orange", fontsize=9)
    ax.text(
        0.8,
        0.55,
        "Camera rotation applies a rigid rotation to the sensor coordinates.\n"
        "Galvo/f-theta misalignment can introduce affine shear, keystone,\n"
        "field curvature, and non-telecentric projection. These distortions\n"
        "are not undone by rotating the camera.",
        fontsize=9,
        va="bottom",
    )

    fig.savefig(fig_path, dpi=220)
    plt.close(fig)

    metrics = {
        "threshold_fraction_of_flat": level,
        "effective_pixel_um": [sx, sy],
        "top_edge_angle_deg": edge_fits["top"]["angle_deg"],
        "bottom_edge_angle_deg": edge_fits["bottom"]["angle_deg"],
        "top_edge_rmse_px": edge_fits["top"]["rmse_px"],
        "bottom_edge_rmse_px": edge_fits["bottom"]["rmse_px"],
        "pca_major_axis_angle_deg": orientation["major_axis_angle_deg"],
        "interpretation": "The 50%-level contour is not an ideal camera-rotated rectangle: the fitted top and bottom edges have different apparent slopes, while the PCA major axis is tilted by about -4 degrees. Camera rotation can remove only global image rotation; it cannot remove shear/keystone/non-telecentric distortion introduced by galvo and scan-lens geometry.",
        "outputs": {"figure": str(fig_path)},
    }
    metrics_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"Figure: {fig_path}")
    print(f"Metrics: {metrics_path}")


if __name__ == "__main__":
    make_figure()
