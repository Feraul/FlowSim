"""Shared validation and indexing for vectorized face-flow kernels."""

from typing import Any, Dict, NamedTuple

import numpy as np

from fs.indexing import to_zero_based


class EdgeIndices(NamedTuple):
    boundary_node1: np.ndarray
    boundary_node2: np.ndarray
    boundary_left: np.ndarray
    internal_node1: np.ndarray
    internal_node2: np.ndarray
    internal_left: np.ndarray
    internal_right: np.ndarray


def edge_indices(
    bedge: np.ndarray,
    inedge: np.ndarray,
    n_nodes: int,
    n_elements: int,
) -> EdgeIndices:
    n_boundary = bedge.shape[0]
    n_internal = inedge.shape[0]
    boundary_nodes = to_zero_based(
        bedge[:, :2], n_nodes, "bedge node columns"
    ).reshape(n_boundary, 2)
    boundary_left = to_zero_based(
        bedge[:, 2], n_elements, "bedge[:, 2]"
    )
    internal_nodes = to_zero_based(
        inedge[:, :2], n_nodes, "inedge node columns"
    ).reshape(n_internal, 2)
    internal_elements = to_zero_based(
        inedge[:, 2:4], n_elements, "inedge element columns"
    ).reshape(n_internal, 2)
    return EdgeIndices(
        boundary_nodes[:, 0],
        boundary_nodes[:, 1],
        boundary_left,
        internal_nodes[:, 0],
        internal_nodes[:, 1],
        internal_elements[:, 0],
        internal_elements[:, 1],
    )


def vector(values: Any, size: int, name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float).reshape(-1)
    if result.size != size:
        raise ValueError(f"{name} has {result.size} entries, expected {size}")
    return result


def face_mobility(
    viscosity: Any, numcase: int, n_boundary: int, n_internal: int
) -> np.ndarray:
    n_faces = n_boundary + n_internal
    if not (30 < numcase < 200 or 200 < numcase < 300):
        return np.ones(n_faces)
    if viscosity is None:
        raise ValueError("config.visc is required for viscosity-dependent cases")

    viscosity = np.asarray(viscosity, dtype=float)
    if viscosity.shape[0] != n_faces:
        raise ValueError(f"config.visc has {viscosity.shape[0]} rows, expected {n_faces}")
    if 30 < numcase < 200:
        return np.sum(viscosity.reshape(n_faces, -1), axis=1)
    if numcase in (245, 246, 247, 248, 249, 251):
        if viscosity.size != n_faces:
            raise ValueError("config.visc must have one value per face for this case")
        return viscosity.reshape(-1)
    return np.ones(n_faces)


def apply_neumann(
    boundary_flow: np.ndarray,
    edge_length: np.ndarray,
    boundary_flags: np.ndarray,
    config: Dict[str, Any],
    flowrate_z: np.ndarray,
) -> None:
    neumann = boundary_flags >= 200
    if not np.any(neumann):
        return

    bcflag = np.asarray(config["bcflag"], dtype=float)
    flux_by_flag = {int(flag): value for flag, value in bcflag[:, :2]}
    missing = sorted(set(boundary_flags[neumann]) - flux_by_flag.keys())
    if missing:
        raise ValueError(f"Neumann flags missing from config.bcflag: {missing}")
    prescribed = np.fromiter(
        (flux_by_flag[flag] for flag in boundary_flags[neumann]),
        dtype=float,
        count=np.count_nonzero(neumann),
    )
    gravity_boundary = flowrate_z[: boundary_flags.size][boundary_flags > 200]
    if gravity_boundary.size != prescribed.size:
        raise ValueError("Neumann boundary flags must be greater than 200")
    boundary_flow[neumann] = (
        -edge_length[neumann] * prescribed + gravity_boundary
    )


def element_balance(
    flowrate: np.ndarray,
    boundary_left: np.ndarray,
    internal_left: np.ndarray,
    internal_right: np.ndarray,
    n_elements: int,
) -> np.ndarray:
    n_boundary = boundary_left.size
    result = np.bincount(
        boundary_left, weights=flowrate[:n_boundary], minlength=n_elements
    )
    result += np.bincount(
        internal_left, weights=flowrate[n_boundary:], minlength=n_elements
    )
    result -= np.bincount(
        internal_right, weights=flowrate[n_boundary:], minlength=n_elements
    )
    return result
