# -*- coding: utf-8 -*-
"""
fs.lpew.v2.angulos — Batched LPEW2 corner angles (translated from +fs/+lpew/+v2/angulos.m)

Inputs:
- FS (dict), expects FS['csr']['tCurrent'], FS['csr']['tNext'] (arrays length nCorners)
- T_all: (nAllNeighbors, dim)
- O_all: (nCorners, dim)
- Q_corner: (nCorners, dim)

Outputs:
- ve2, ve1, theta2, theta1: each (nCorners,) numpy arrays
"""
from typing import Dict, Any, Tuple
import numpy as np


def angulos(FS: Dict[str, Any], T_all: np.ndarray, O_all: np.ndarray, Q_corner: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    csr = FS.get('csr', {})
    tCurrent = np.asarray(csr.get('tCurrent')).astype(int)
    tNext = np.asarray(csr.get('tNext')).astype(int)

    Tk = np.asarray(T_all)[tCurrent, :]
    Tk1 = np.asarray(T_all)[tNext, :]

    v0 = np.asarray(O_all) - np.asarray(Q_corner)
    vth1 = Tk - np.asarray(Q_corner)
    vth2 = Tk1 - np.asarray(Q_corner)
    v1 = Tk1 - Tk

    n_v0 = np.linalg.norm(v0, axis=1)
    n_v1 = np.linalg.norm(v1, axis=1)
    n_vth1 = np.linalg.norm(vth1, axis=1)
    n_vth2 = np.linalg.norm(vth2, axis=1)

    # Dot products
    d_vth1_v1 = -np.sum(vth1 * v1, axis=1)
    d_vth2_v1 = np.sum(vth2 * v1, axis=1)
    d_v0_vth1 = np.sum(v0 * vth1, axis=1)
    d_v0_vth2 = np.sum(v0 * vth2, axis=1)

    # Guard against division by zero
    denom_ve1 = np.where((n_v1 * n_vth1) == 0, 1.0, (n_v1 * n_vth1))
    denom_ve2 = np.where((n_v1 * n_vth2) == 0, 1.0, (n_v1 * n_vth2))
    denom_th1 = np.where((n_v0 * n_vth1) == 0, 1.0, (n_v0 * n_vth1))
    denom_th2 = np.where((n_v0 * n_vth2) == 0, 1.0, (n_v0 * n_vth2))

    ve1 = np.arccos(np.clip(d_vth1_v1 / denom_ve1, -1.0, 1.0))
    ve2 = np.arccos(np.clip(d_vth2_v1 / denom_ve2, -1.0, 1.0))
    theta1 = np.arccos(np.clip(d_v0_vth1 / denom_th1, -1.0, 1.0))
    theta2 = np.arccos(np.clip(d_v0_vth2 / denom_th2, -1.0, 1.0))

    return ve2, ve1, theta2, theta1
