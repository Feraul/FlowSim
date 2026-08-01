"""Batched LPEW2 corner angles translated from ``+fs/+lpew/+v2/angulos.m``."""

from typing import Any, Dict, Tuple

import numpy as np


def angulos(
    FS: Dict[str, Any],
    T_all: np.ndarray,
    O_all: np.ndarray,
    Q_corner: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``ve2, ve1, theta2, theta1`` for every mesh corner."""
    csr = FS["csr"]
    t_current = np.asarray(csr["tCurrent"], dtype=int)
    t_next = np.asarray(csr["tNext"], dtype=int)

    T_all = np.asarray(T_all, dtype=float)
    O_all = np.asarray(O_all, dtype=float)
    Q_corner = np.asarray(Q_corner, dtype=float)

    tk = T_all[t_current]
    tk1 = T_all[t_next]

    v0 = O_all - Q_corner
    vth1 = tk - Q_corner
    vth2 = tk1 - Q_corner
    v1 = tk1 - tk

    n_v0 = np.linalg.norm(v0, axis=1)
    n_v1 = np.linalg.norm(v1, axis=1)
    n_vth1 = np.linalg.norm(vth1, axis=1)
    n_vth2 = np.linalg.norm(vth2, axis=1)

    d_vth1_v1 = -np.sum(vth1 * v1, axis=1)
    d_vth2_v1 = np.sum(vth2 * v1, axis=1)
    d_v0_vth1 = np.sum(v0 * vth1, axis=1)
    d_v0_vth2 = np.sum(v0 * vth2, axis=1)

    ve1 = np.arccos(d_vth1_v1 / (n_v1 * n_vth1))
    ve2 = np.arccos(d_vth2_v1 / (n_v1 * n_vth2))
    theta1 = np.arccos(d_v0_vth1 / (n_v0 * n_vth1))
    theta2 = np.arccos(d_v0_vth2 / (n_v0 * n_vth2))

    return ve2, ve1, theta2, theta1
