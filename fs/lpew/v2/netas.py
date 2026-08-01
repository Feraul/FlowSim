# -*- coding: utf-8 -*-
"""
fs.lpew.v2.netas — Batched LPEW netas (translated from +fs/+lpew/+v2/netas.m)

Inputs:
- FS: dict with FS['csr']['tCurrent'], FS['csr']['tNext']
- T_all: (nAllNeighbors, dim)
- P_all: (nAllNeighbors, dim)
- O_all: (nCorners, dim)
- Q_corner: (nCorners, dim)

Output:
- out: (nCorners, 2) array with columns [Tk_dist/h1, Tk1_dist/h2]

Notes:
- Expects 0-based indices. Pads 2D vectors to 3D for cross-product compatibility.
- Protects against division by zero with safe denominators.
"""
from typing import Dict, Any
import numpy as np


def _ensure_3d(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    if arr.shape[1] == 2:
        # pad z=0
        z = np.zeros((arr.shape[0], 1), dtype=arr.dtype)
        arr = np.concatenate([arr, z], axis=1)
    return arr


def netas(FS: Dict[str, Any], T_all, P_all, O_all, Q_corner):
    csr = FS.get('csr', {})
    tCurrent = np.asarray(csr.get('tCurrent')).astype(int)
    tNext = np.asarray(csr.get('tNext')).astype(int)

    T_all = _ensure_3d(np.asarray(T_all))
    P_all = _ensure_3d(np.asarray(P_all))
    O_all = _ensure_3d(np.asarray(O_all))
    Q_corner = _ensure_3d(np.asarray(Q_corner))

    Tk = T_all[tCurrent, :]
    Tk1 = T_all[tNext, :]
    Pk = P_all[tCurrent, :]
    Pk1 = P_all[tNext, :]

    v1 = O_all - Q_corner
    v2_col1 = Pk - Q_corner
    v2_col2 = Pk1 - Q_corner

    # cross products (nCorners x 3)
    ce1 = np.cross(v1, v2_col1)
    ce2 = np.cross(v1, v2_col2)

    n_ce1 = np.linalg.norm(ce1, axis=1)
    n_ce2 = np.linalg.norm(ce2, axis=1)
    n_v2_col1 = np.linalg.norm(v2_col1, axis=1)
    n_v2_col2 = np.linalg.norm(v2_col2, axis=1)

    # safe denominators
    denom1 = np.where(n_v2_col1 == 0, 1.0, n_v2_col1)
    denom2 = np.where(n_v2_col2 == 0, 1.0, n_v2_col2)

    h1 = n_ce1 / denom1
    h2 = n_ce2 / denom2

    Tk_dist = np.linalg.norm(Tk - Q_corner, axis=1)
    Tk1_dist = np.linalg.norm(Tk1 - Q_corner, axis=1)

    # avoid division by zero in final ratio
    h1_safe = np.where(h1 == 0, 1.0, h1)
    h2_safe = np.where(h2 == 0, 1.0, h2)

    out = np.column_stack((Tk_dist / h1_safe, Tk1_dist / h2_safe))
    return out
