"""Batched LPEW netas translated from ``+fs/+lpew/+v2/netas.m``."""

from typing import Any, Dict

import numpy as np


def netas(
    FS: Dict[str, Any],
    T_all: np.ndarray,
    P_all: np.ndarray,
    O_all: np.ndarray,
    Q_corner: np.ndarray,
) -> np.ndarray:
    """Return the two geometric eta coefficients for every mesh corner."""
    csr = FS["csr"]
    t_current = np.asarray(csr["tCurrent"], dtype=int)
    t_next = np.asarray(csr["tNext"], dtype=int)

    T_all = np.asarray(T_all, dtype=float)
    P_all = np.asarray(P_all, dtype=float)
    O_all = np.asarray(O_all, dtype=float)
    Q_corner = np.asarray(Q_corner, dtype=float)

    tk = T_all[t_current]
    tk1 = T_all[t_next]
    pk = P_all[t_current]
    pk1 = P_all[t_next]

    v1 = O_all - Q_corner
    v2_col1 = pk - Q_corner
    v2_col2 = pk1 - Q_corner

    h1 = _cross_norm(v1, v2_col1) / np.linalg.norm(v2_col1, axis=1)
    h2 = _cross_norm(v1, v2_col2) / np.linalg.norm(v2_col2, axis=1)

    tk_dist = np.linalg.norm(tk - Q_corner, axis=1)
    tk1_dist = np.linalg.norm(tk1 - Q_corner, axis=1)

    return np.column_stack((tk_dist / h1, tk1_dist / h2))


def _cross_norm(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    if left.shape[1] == 2 and right.shape[1] == 2:
        return np.abs(left[:, 0] * right[:, 1] - left[:, 1] * right[:, 0])
    if left.shape[1] == 3 and right.shape[1] == 3:
        return np.linalg.norm(np.cross(left, right), axis=1)
    raise ValueError("netas expects coordinates with two or three dimensions")
