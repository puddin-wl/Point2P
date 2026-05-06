"""Small utilities for standalone sweep runs."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


def timestamp() -> str:
    """Return a filesystem-safe local timestamp."""
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def ensure_unique_dir(path: str | Path) -> Path:
    """Create a directory without overwriting an existing run."""
    base = Path(path)
    candidate = base
    index = 1
    while candidate.exists():
        candidate = base.with_name(f"{base.name}_{index:02d}")
        index += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def center_crop(array: np.ndarray, size: int) -> np.ndarray:
    """Return a centered square crop."""
    if size <= 0:
        raise ValueError("Crop size must be positive.")
    ny, nx = array.shape[:2]
    if size > ny or size > nx:
        raise ValueError(f"Crop size {size} exceeds array shape {(ny, nx)}.")
    y0 = ny // 2 - size // 2
    x0 = nx // 2 - size // 2
    return np.asarray(array[y0 : y0 + size, x0 : x0 + size])


def load_json(path: str | Path) -> dict[str, Any]:
    """Load JSON with UTF-8 encoding."""
    return json.loads(Path(path).read_text(encoding="utf-8"))


def to_jsonable(value: Any) -> Any:
    """Convert NumPy, Path, and scalar values into JSON-friendly data."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    return value


def save_json(path: str | Path, data: dict[str, Any]) -> None:
    """Save an indented JSON file."""
    Path(path).write_text(json.dumps(to_jsonable(data), indent=2), encoding="utf-8")


def write_csv(path: str | Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    """Write summary rows using a stable field order."""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def value_label(value: Any) -> str:
    """Return a compact label for scalar or pair sweep values."""
    if isinstance(value, (list, tuple)):
        return "(" + ", ".join(f"{float(v):g}" for v in value) + ")"
    if isinstance(value, np.ndarray):
        return value_label(value.tolist())
    if isinstance(value, (int, float, np.floating, np.integer)):
        return f"{float(value):g}"
    return str(value)


def jsonish_to_dict(value: Any) -> dict[str, Any]:
    """Decode JSON-like scalar fields saved in NPZ files."""
    if isinstance(value, np.ndarray):
        if value.dtype.kind in {"U", "S"}:
            value = "".join(str(v) for v in value.ravel())
            return json.loads(value)
        value = value.squeeze()
        if value.shape == ():
            value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        return json.loads(value)
    return dict(value)

