"""Shared iterator validation."""

from typing import Any

import numpy as np


class ConvergenceError(RuntimeError):
    """Raised when an iterative solver exhausts its iteration budget."""


def vector(values: Any, size: int = -1, name: str = "value") -> np.ndarray:
    result = np.asarray(values, dtype=float).reshape(-1)
    if size >= 0 and result.size != size:
        raise ValueError(f"{name} has {result.size} entries, expected {size}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} contains non-finite values")
    return result


def validate_options(tolerance: float, max_iterations: int) -> None:
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least one")
