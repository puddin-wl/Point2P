"""Small utilities for run management and reports."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


def timestamp() -> str:
    """Return a filesystem-safe timestamp."""
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def ensure_unique_dir(path: str | Path) -> Path:
    """Create a new output directory without overwriting old results."""
    base = Path(path)
    candidate = base
    i = 1
    while candidate.exists():
        candidate = base.with_name(f"{base.name}_{i:02d}")
        i += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def to_jsonable(value: Any) -> Any:
    """Convert NumPy and path objects to JSON-friendly data."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def save_json(path: str | Path, data: dict[str, Any]) -> None:
    """Save indented JSON."""
    Path(path).write_text(json.dumps(to_jsonable(data), indent=2), encoding="utf-8")


def save_metrics_csv(path: str | Path, metrics_history: list[dict[str, Any]]) -> None:
    """Save metrics history as CSV."""
    if not metrics_history:
        Path(path).write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in metrics_history:
        for key in row:
            if key not in keys:
                keys.append(key)
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in metrics_history:
            writer.writerow(row)


def latest_existing_phase(root: str | Path) -> Path | None:
    """Find the newest phase0.mat under the old initial phase artifacts, if any."""
    root = Path(root)
    matches = sorted(root.glob("initial_phase_generation/artifacts/*/phase0.mat"))
    return matches[-1] if matches else None


def format_metrics(metrics: dict[str, Any]) -> str:
    """Format metrics for text reports."""
    lines = []
    for key, value in metrics.items():
        if isinstance(value, float):
            lines.append(f"{key}: {value:.8g}")
        else:
            lines.append(f"{key}: {value}")
    return "\n".join(lines)
