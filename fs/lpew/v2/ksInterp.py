"""Permeability projections translated from ``+fs/+lpew/+v2/ksInterp.m``."""

from typing import Any, Dict, Tuple

import numpy as np


def ksInterp(
    FS: Dict[str, Any],
    T_all: np.ndarray,
    Q_corner: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return tangential and normal permeability projections per corner."""
    if FS.get("cfg", {}).get("phasekey") == 2:
        raise NotImplementedError(
            "Two-phase mobility scaling has not been vectorized for ksInterp"
        )

    mesh = FS["mesh"]
    csr = FS["csr"]
    kmap = np.asarray(FS["perm"]["kmap"], dtype=float)
    elem = np.asarray(mesh["elem"])
    corner_elem = np.asarray(csr["cornerElem"], dtype=int)
    t_current = np.asarray(csr["tCurrent"], dtype=int)
    t_next = np.asarray(csr["tNext"], dtype=int)

    material_ids = np.asarray(elem[corner_elem, 4], dtype=int)
    material_rows = _material_rows(kmap, material_ids)

    k11 = kmap[material_rows, 1]
    k12 = kmap[material_rows, 2]
    k21 = kmap[material_rows, 3]
    k22 = kmap[material_rows, 4]

    T_all = np.asarray(T_all, dtype=float)
    Q_corner = np.asarray(Q_corner, dtype=float)

    kt1_col1, kn1_col1 = _project(
        T_all[t_current] - Q_corner, k11, k12, k21, k22
    )
    kt1_col2, kn1_col2 = _project(
        T_all[t_next] - Q_corner, k11, k12, k21, k22
    )
    kt2, kn2 = _project(
        T_all[t_next] - T_all[t_current], k11, k12, k21, k22
    )

    kt1 = np.column_stack((kt1_col1, kt1_col2))
    kn1 = np.column_stack((kn1_col1, kn1_col2))
    return kt1, kt2, kn1, kn2


def _material_rows(kmap: np.ndarray, material_ids: np.ndarray) -> np.ndarray:
    if kmap.ndim != 2 or kmap.shape[1] < 5:
        raise ValueError("FS.perm.kmap must contain material id and K11/K12/K21/K22")

    row_by_id = {int(material_id): row for row, material_id in enumerate(kmap[:, 0])}
    missing = sorted(set(material_ids) - row_by_id.keys())
    if missing:
        raise ValueError(f"Material ids missing from FS.perm.kmap: {missing}")
    return np.fromiter(
        (row_by_id[material_id] for material_id in material_ids),
        dtype=int,
        count=material_ids.size,
    )


def _project(
    vectors: np.ndarray,
    k11: np.ndarray,
    k12: np.ndarray,
    k21: np.ndarray,
    k22: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    vx = vectors[:, 0]
    vy = vectors[:, 1]
    length_squared = vx * vx + vy * vy
    if np.any(length_squared == 0):
        raise ValueError("Cannot project permeability along a zero-length vector")

    normal = (
        k11 * vy * vy
        - (k12 + k21) * vx * vy
        + k22 * vx * vx
    ) / length_squared
    tangential = (
        (k11 - k22) * vx * vy
        + k12 * vy * vy
        - k21 * vx * vx
    ) / length_squared
    return tangential, normal
