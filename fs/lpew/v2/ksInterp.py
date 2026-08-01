# -*- coding: utf-8 -*-
"""
fs.lpew.v2.ksInterp — Batched permeability-tensor projections (translated from +fs/+lpew/+v2/ksInterp.m)

Inputs:
- FS: dict-like structure with keys:
    FS['mesh']['elem'] (nElems x ?), material id in column index 4 (0-based)
    FS['perm']['kmap'] (nElems x >=5) - columns 1..4 -> K11,K12,K21,K22 (0-based indices)
    FS['csr']['cornerElem'], FS['csr']['tCurrent'], FS['csr']['tNext'] (arrays of corner length)
- T_all: (nAllNeighbors x 2 or 3) midpoints from OPT
- Q_corner: (nCorners x 2 or 3) owning node coord per corner

Returns:
- Kt1: (nCorners, 2)
- Kt2: (nCorners,)
- Kn1: (nCorners, 2)
- Kn2: (nCorners,)

Notes:
- Expects 0-based indices everywhere.
- For phasekey == 2 (two-phase), raises NotImplementedError as in MATLAB.
"""
from typing import Dict, Any, Tuple
import numpy as np


def ksInterp(FS: Dict[str, Any], T_all: np.ndarray, Q_corner: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mesh = FS.get('mesh', {})
    perm = FS.get('perm', {})
    csr = FS.get('csr', {})
    cfg = FS.get('cfg', {})

    elem = np.asarray(mesh.get('elem'))
    kmap = np.asarray(perm.get('kmap'))
    cornerElem = np.asarray(csr.get('cornerElem')).astype(int)
    tCurrent = np.asarray(csr.get('tCurrent')).astype(int)
    tNext = np.asarray(csr.get('tNext')).astype(int)

    phasekey = cfg.get('phasekey', None)
    if phasekey == 2:
        raise NotImplementedError('Two-phase (phasekey=2) mobility scaling not yet vectorized.')

    # material ids per corner (MAT ids are stored in elem[:,4] in MATLAB 1-based; here assume 0-based indexing into kmap)
    # elem may be 2D array; material id column index is 4 (0-based) per repository comment
    matIds = elem[cornerElem, 4].astype(int)

    # Extract tensor components per corner from kmap (assumes kmap rows align with material ids)
    # MATLAB: K11 = kmap(matIds, 2) meaning column index 2 (1-based) -> in 0-based python it's column 1
    # However earlier code used kmap(matIds,2..5) — repository README said columns 2-5 are K11..K22. In Python kmap columns assumed same ordering.
    # Here we assume kmap columns: [idx, K11, K12, K21, K22, ...] so K11 = kmap[matIds,1]
    # Adjust if your kmap layout differs.
    K11 = kmap[matIds, 1]
    K12 = kmap[matIds, 2]
    K21 = kmap[matIds, 3]
    K22 = kmap[matIds, 4]

    # Ensure arrays
    T_all = np.asarray(T_all)
    Q_corner = np.asarray(Q_corner)

    # Column i=1: use T(k) i.e. tCurrent
    v_k = T_all[tCurrent, :] - Q_corner  # shape (nCorners, dim)
    vx = v_k[:, 0]
    vy = v_k[:, 1]
    len2_k = vx * vx + vy * vy
    # avoid division by zero
    denom_k = np.where(len2_k == 0, 1.0, len2_k)

    Kn1_c1 = (K11 * vy * vy - (K12 + K21) * vx * vy + K22 * vx * vx) / denom_k
    Kt1_c1 = ((K11 - K22) * vx * vy + K12 * vy * vy - K21 * vx * vx) / denom_k

    # Column i=2: use T(k+1) -> tNext
    v_k1 = T_all[tNext, :] - Q_corner
    vx1 = v_k1[:, 0]
    vy1 = v_k1[:, 1]
    len2_k1 = vx1 * vx1 + vy1 * vy1
    denom_k1 = np.where(len2_k1 == 0, 1.0, len2_k1)

    Kn1_c2 = (K11 * vy1 * vy1 - (K12 + K21) * vx1 * vy1 + K22 * vx1 * vx1) / denom_k1
    Kt1_c2 = ((K11 - K22) * vx1 * vy1 + K12 * vy1 * vy1 - K21 * vx1 * vx1) / denom_k1

    Kn1 = np.column_stack((Kn1_c1, Kn1_c2))
    Kt1 = np.column_stack((Kt1_c1, Kt1_c2))

    # Kn2, Kt2: use v = T(k+1) - T(k)
    v_diff = T_all[tNext, :] - T_all[tCurrent, :]
    vx_d = v_diff[:, 0]
    vy_d = v_diff[:, 1]
    len2_diff = vx_d * vx_d + vy_d * vy_d
    denom_diff = np.where(len2_diff == 0, 1.0, len2_diff)

    Kn2 = (K11 * vy_d * vy_d - (K12 + K21) * vx_d * vy_d + K22 * vx_d * vx_d) / denom_diff
    Kt2 = ((K11 - K22) * vx_d * vy_d + K12 * vy_d * vy_d - K21 * vx_d * vx_d) / denom_diff

    return Kt1, Kt2, Kn1, Kn2
