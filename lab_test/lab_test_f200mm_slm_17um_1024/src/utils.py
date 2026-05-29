"""Small utilities for run management and reports.
---
时间戳、目录创建、JSON/CSV 导出、指标格式化等辅助工具。
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


def timestamp() -> str:
    """生成文件名安全的时间戳 (YYYYmmdd-HHMMSS)."""
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def ensure_unique_dir(path: str | Path) -> Path:
    """创建不重名的输出目录 (已存在则追加 _01, _02...)."""
    base = Path(path)
    candidate = base
    i = 1
    while candidate.exists():
        candidate = base.with_name(f"{base.name}_{i:02d}")
        i += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def to_jsonable(value: Any) -> Any:
    """将 NumPy/Path 等对象递归转为 JSON 可序列化格式."""
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
    """保存带缩进的 JSON 文件."""
    Path(path).write_text(json.dumps(to_jsonable(data), indent=2), encoding="utf-8")


def save_metrics_csv(path: str | Path, metrics_history: list[dict[str, Any]]) -> None:
    """将指标历史保存为 CSV (每行一个迭代)."""
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
    """在 artifacts 目录中查找最新的 phase0.mat, 找不到返回 None."""
    root = Path(root)
    matches = sorted(root.glob("initial_phase_generation/artifacts/*/phase0.mat"))
    return matches[-1] if matches else None


def format_metrics(metrics: dict[str, Any]) -> str:
    """将指标字典格式化为可读的文本字符串."""
    lines = []
    for key, value in metrics.items():
        if isinstance(value, float):
            lines.append(f"{key}: {value:.8g}")
        else:
            lines.append(f"{key}: {value}")
    return "\n".join(lines)
