"""Index normalization shared by Python translations of legacy MATLAB data."""

from typing import Any

import numpy as np


def to_zero_based(values: Any, size: int, name: str) -> np.ndarray:
    """Normalize a flat legacy index array while preserving zero-based input."""
    indices = np.asarray(values, dtype=int).reshape(-1)
    if indices.size == 0:
        return indices
    if np.any(indices < 0):
        raise ValueError(f"{name} contains a negative index")
    if np.any(indices == 0):
        if np.any(indices >= size):
            raise ValueError(f"{name} contains an index outside 0..{size - 1}")
        return indices
    if np.any(indices > size):
        raise ValueError(f"{name} contains an index outside 1..{size}")
    return indices - 1
