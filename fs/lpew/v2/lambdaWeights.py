# -*- coding: utf-8 -*-
"""
fs.lpew.v2.lambdaWeights — Batched LPEW2 lambda and r (translated from +fs/+lpew/+v2/lambdaWeights.m)

This implements the interior fully-vectorized formula and a per-node fallback
for boundary nodes as in the MATLAB original.

Inputs:
- FS: dict with FS['csr'] and FS['mesh'] structures
- Kt1, Kt2, Kn1, Kn2: arrays per corner (Kt1: nCorners x 2, Kn1: nCorners x 2, Kt2/Kn2: nCorners)
- theta1, theta2, ve1, ve2: (nCorners,)
- netas: (nCorners, 2)
- T_all: (nAllNeighbors, dim)
- Q_corner: (nCorners, dim)

Returns:
- lambda_arr: (nCorners,) per-corner lambda
- r: (nNodes, 2) boundary Neumann correction (zeros for interior nodes)
"""
from typing import Dict, Any, Tuple
import numpy as np


def _cot(x: np.ndarray) -> np.ndarray:
    # safe cotangent: 1 / tan(x) with protection against tan==0
    x = np.asarray(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.tan(x)
        cotx = 1.0 / t
        # replace infinite values with large finite numbers
        cotx = np.where(np.isfinite(cotx), cotx, np.sign(cotx) * 1e12)
    return cotx


def lambdaWeights(FS: Dict[str, Any], Kt1: np.ndarray, Kt2: np.ndarray, Kn1: np.ndarray, Kn2: np.ndarray,
                  theta1: np.ndarray, theta2: np.ndarray, ve1: np.ndarray, ve2: np.ndarray,
                  netas: np.ndarray, T_all: np.ndarray, Q_corner: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    csr = FS.get('csr', {})
    mesh = FS.get('mesh', {})

    nCorners = int(csr.get('nCorners'))
    nNodes = int(mesh.get('nNodes'))

    Kt1 = np.asarray(Kt1)
    Kt2 = np.asarray(Kt2).reshape(-1)
    Kn1 = np.asarray(Kn1)
    Kn2 = np.asarray(Kn2).reshape(-1)
    theta1 = np.asarray(theta1).reshape(-1)
    theta2 = np.asarray(theta2).reshape(-1)
    ve1 = np.asarray(ve1).reshape(-1)
    ve2 = np.asarray(ve2).reshape(-1)
    netas = np.asarray(netas)
    T_all = np.asarray(T_all)
    Q_corner = np.asarray(Q_corner)

    lambda_arr = np.zeros(nCorners, dtype=float)
    r = np.zeros((nNodes, 2), dtype=float)

    # Interior mask per corner: nodeIsInterior mapped by cornerNode
    cornerNode = np.asarray(csr.get('cornerNode')).astype(int)
    nodeIsInterior = np.asarray(csr.get('nodeIsInterior')).astype(bool)  # per node
    isInterior = nodeIsInterior[cornerNode]

    iprev = np.asarray(csr.get('cornerPrev')).astype(int)
    inext = np.asarray(csr.get('cornerNext')).astype(int)

    # Precompute cotangents
    cot_ve1_prev = _cot(ve1[iprev])
    cot_ve2_cur = _cot(ve2)
    cot_theta2_prev = _cot(theta2[iprev])
    cot_theta1_cur = _cot(theta1)

    # zeta numerator and denominator (vectorized for all corners)
    zeta_num = Kn2[iprev] * cot_ve1_prev + Kn2 * cot_ve2_cur + Kt2[iprev] - Kt2
    zeta_den = Kn1[iprev, 1] * cot_theta2_prev + Kn1[:, 0] * cot_theta1_cur - Kt1[iprev, 1] + Kt1[:, 0]

    # safe division
    zeta_den_safe = np.where(zeta_den == 0, 1.0, zeta_den)
    zeta = zeta_num / zeta_den_safe

    # lambda for interior corners
    lambda_int = (Kn1[:, 0] * netas[:, 0] * zeta) + (Kn1[:, 1] * netas[:, 1] * zeta[inext])
    # assign only for interior corners (others will be set in boundary loop)
    lambda_arr[isInterior] = lambda_int[isInterior]

    # Boundary nodes: per-node fallback (legacy path)
    # esurn2 / nsurn2 are prefix-sum arrays; in MATLAB indexing used +1 slices, here use 0-based slices
    esurn2 = np.asarray(mesh.get('esurn2')).astype(int)
    nsurn2 = np.asarray(mesh.get('nsurn2')).astype(int)
    coord = np.asarray(mesh.get('coord'))

    # boundary node indices
    boundary_nodes = np.nonzero(~nodeIsInterior)[0]

    for No in boundary_nodes:
        start_e = int(esurn2[No])
        end_e = int(esurn2[No + 1])
        e_range = np.arange(start_e, end_e, dtype=int)

        start_p = int(nsurn2[No])
        end_p = int(nsurn2[No + 1])
        p_range = np.arange(start_p, end_p, dtype=int)

        nec = e_range.size
        if nec == 0:
            continue

        Qo = coord[No, :]

        # per-node slices
        Kt1_n = Kt1[e_range, :]
        Kt2_n = Kt2[e_range]
        Kn1_n = Kn1[e_range, :]
        Kn2_n = Kn2[e_range]
        th1_n = theta1[e_range]
        th2_n = theta2[e_range]
        ve1_n = ve1[e_range]
        ve2_n = ve2[e_range]
        netas_n = netas[e_range, :]
        T_n = T_all[p_range, :]

        # compute zeta_n length nec+1
        zeta_n = np.zeros(nec + 1, dtype=float)
        # we will compute norms using numpy.linalg.norm
        for kk in range(0, nec + 1):
            if kk == 0:
                # MATLAB kk==1 case
                zn = Kn2_n[0] * _cot(ve2_n[0]) - Kt2_n[0]
                zd = Kn1_n[0, 0] * _cot(th1_n[0]) + Kt1_n[0, 0]
                denom = zd if zd != 0 else 1.0
                zeta_n[kk] = zn / denom
                r[No, 0] = (1.0 + zeta_n[kk]) * np.linalg.norm(Qo - T_n[0, :])
            elif kk == nec:
                # MATLAB kk==nec+1 case
                idx = nec - 1
                zn = Kn2_n[idx] * _cot(ve1_n[idx]) + Kt2_n[idx]
                zd = Kn1_n[idx, 1] * _cot(th2_n[idx]) - Kt1_n[idx, 1]
                denom = zd if zd != 0 else 1.0
                zeta_n[kk] = zn / denom
                r[No, 1] = (1.0 + zeta_n[kk]) * np.linalg.norm(Qo - T_n[nec, :])
            else:
                # interior positions within node
                idxm = kk - 1
                idx = kk
                zn = (Kn2_n[idxm] * _cot(ve1_n[idxm]) + Kn2_n[idx] * _cot(ve2_n[idx]) +
                      Kt2_n[idxm] - Kt2_n[idx])
                zd = (Kn1_n[idxm, 1] * _cot(th2_n[idxm]) + Kn1_n[idx, 0] * _cot(th1_n[idx]) -
                      Kt1_n[idxm, 1] + Kt1_n[idx, 0])
                denom = zd if zd != 0 else 1.0
                zeta_n[kk] = zn / denom

        # Lambda for boundary corners with wrap
        for kk in range(0, nec):
            if kk == nec - 1:
                lambda_arr[e_range[kk]] = (Kn1_n[kk, 0] * netas_n[kk, 0] * zeta_n[kk] +
                                           Kn1_n[kk, 1] * netas_n[kk, 1] * zeta_n[0])
            else:
                lambda_arr[e_range[kk]] = (Kn1_n[kk, 0] * netas_n[kk, 0] * zeta_n[kk] +
                                           Kn1_n[kk, 1] * netas_n[kk, 1] * zeta_n[kk + 1])

    return lambda_arr, r
