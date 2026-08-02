"""LPEW2 lambda weights translated from ``+fs/+lpew/+v2/lambdaWeights.m``."""

from typing import Any, Dict, Tuple

import numpy as np


def lambdaWeights(
    FS: Dict[str, Any],
    Kt1: np.ndarray,
    Kt2: np.ndarray,
    Kn1: np.ndarray,
    Kn2: np.ndarray,
    theta1: np.ndarray,
    theta2: np.ndarray,
    ve1: np.ndarray,
    ve2: np.ndarray,
    netas: np.ndarray,
    T_all: np.ndarray,
    Q_corner: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute per-corner lambda values and boundary Neumann corrections."""
    mesh = FS["mesh"]
    csr = FS["csr"]
    n_corners = int(csr.get("nCorners", len(csr["cornerNode"])))
    n_nodes = int(mesh["nNodes"])

    Kt1 = np.asarray(Kt1, dtype=float)
    Kt2 = np.asarray(Kt2, dtype=float).reshape(-1)
    Kn1 = np.asarray(Kn1, dtype=float)
    Kn2 = np.asarray(Kn2, dtype=float).reshape(-1)
    theta1 = np.asarray(theta1, dtype=float).reshape(-1)
    theta2 = np.asarray(theta2, dtype=float).reshape(-1)
    ve1 = np.asarray(ve1, dtype=float).reshape(-1)
    ve2 = np.asarray(ve2, dtype=float).reshape(-1)
    netas = np.asarray(netas, dtype=float)
    T_all = np.asarray(T_all, dtype=float)
    Q_corner = np.asarray(Q_corner, dtype=float)

    corner_node = np.asarray(csr["cornerNode"], dtype=int)
    corner_prev = np.asarray(csr["cornerPrev"], dtype=int)
    corner_next = np.asarray(csr["cornerNext"], dtype=int)
    node_is_interior = np.asarray(csr["nodeIsInterior"], dtype=bool)
    interior_corners = node_is_interior[corner_node]

    lambda_values = np.zeros(n_corners, dtype=float)
    r = np.zeros((n_nodes, 2), dtype=float)

    zeta_num = (
        Kn2[corner_prev] * _cot(ve1[corner_prev])
        + Kn2 * _cot(ve2)
        + Kt2[corner_prev]
        - Kt2
    )
    zeta_den = (
        Kn1[corner_prev, 1] * _cot(theta2[corner_prev])
        + Kn1[:, 0] * _cot(theta1)
        - Kt1[corner_prev, 1]
        + Kt1[:, 0]
    )
    zeta = zeta_num / zeta_den
    lambda_interior = (
        Kn1[:, 0] * netas[:, 0] * zeta
        + Kn1[:, 1] * netas[:, 1] * zeta[corner_next]
    )
    lambda_values[interior_corners] = lambda_interior[interior_corners]

    esurn2 = np.asarray(mesh["esurn2"], dtype=int)
    nsurn2 = np.asarray(mesh["nsurn2"], dtype=int)
    coord = np.asarray(mesh["coord"], dtype=float)

    for node in np.flatnonzero(~node_is_interior):
        corner_slice = slice(esurn2[node], esurn2[node + 1])
        neighbor_slice = slice(nsurn2[node], nsurn2[node + 1])
        corner_indices = np.arange(corner_slice.start, corner_slice.stop)
        n_node_corners = corner_indices.size

        kt1_node = Kt1[corner_slice]
        kt2_node = Kt2[corner_slice]
        kn1_node = Kn1[corner_slice]
        kn2_node = Kn2[corner_slice]
        theta1_node = theta1[corner_slice]
        theta2_node = theta2[corner_slice]
        ve1_node = ve1[corner_slice]
        ve2_node = ve2[corner_slice]
        netas_node = netas[corner_slice]
        t_node = T_all[neighbor_slice]

        if t_node.shape[0] != n_node_corners + 1:
            raise ValueError(
                f"Boundary node {node} must have one more neighbor than corner"
            )

        zeta_node = np.empty(n_node_corners + 1, dtype=float)
        first_num = kn2_node[0] * _cot(ve2_node[0]) - kt2_node[0]
        first_den = kn1_node[0, 0] * _cot(theta1_node[0]) + kt1_node[0, 0]
        zeta_node[0] = first_num / first_den
        r[node, 0] = (1.0 + zeta_node[0]) * np.linalg.norm(
            coord[node] - t_node[0]
        )

        for local_corner in range(1, n_node_corners):
            previous = local_corner - 1
            numerator = (
                kn2_node[previous] * _cot(ve1_node[previous])
                + kn2_node[local_corner] * _cot(ve2_node[local_corner])
                + kt2_node[previous]
                - kt2_node[local_corner]
            )
            denominator = (
                kn1_node[previous, 1] * _cot(theta2_node[previous])
                + kn1_node[local_corner, 0] * _cot(theta1_node[local_corner])
                - kt1_node[previous, 1]
                + kt1_node[local_corner, 0]
            )
            zeta_node[local_corner] = numerator / denominator

        last_num = kn2_node[-1] * _cot(ve1_node[-1]) + kt2_node[-1]
        last_den = kn1_node[-1, 1] * _cot(theta2_node[-1]) - kt1_node[-1, 1]
        zeta_node[-1] = last_num / last_den
        r[node, 1] = (1.0 + zeta_node[-1]) * np.linalg.norm(
            coord[node] - t_node[-1]
        )

        next_zeta = np.concatenate((zeta_node[1:n_node_corners], zeta_node[:1]))
        lambda_values[corner_slice] = (
            kn1_node[:, 0] * netas_node[:, 0] * zeta_node[:n_node_corners]
            + kn1_node[:, 1] * netas_node[:, 1] * next_zeta
        )

    return lambda_values, r


def _cot(values: np.ndarray) -> np.ndarray:
    return 1.0 / np.tan(values)
