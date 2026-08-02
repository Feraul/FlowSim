"""Vectorized TPFA flow rate translated from ``+fs/+flow/tpfa.m``."""

from typing import Any, Dict, Tuple

import numpy as np

from fs.indexing import to_zero_based


def tpfa(
    p: np.ndarray, env: Dict[str, Any]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Compute face flow rates, per-element balance, and dispersive flow."""
    geometry = env["geometry"]
    config = env["config"]
    method = env["premethod"]["TPFA"]

    coord = np.asarray(geometry["coord"], dtype=float)
    bedge = np.asarray(geometry["bedge"])
    inedge = np.asarray(geometry["inedge"])
    cent_elem = np.asarray(geometry["centelem"], dtype=float)
    p = np.asarray(p, dtype=float).reshape(-1)

    n_nodes = coord.shape[0]
    n_elements = cent_elem.shape[0]
    n_boundary = bedge.shape[0]
    n_internal = inedge.shape[0]
    if p.size != n_elements:
        raise ValueError(f"p has {p.size} entries, expected {n_elements}")

    boundary_nodes = to_zero_based(
        bedge[:, :2], n_nodes, "bedge node columns"
    ).reshape(n_boundary, 2)
    boundary_node1 = boundary_nodes[:, 0]
    boundary_node2 = boundary_nodes[:, 1]
    boundary_left = to_zero_based(bedge[:, 2], n_elements, "bedge[:, 2]")
    internal_elements = to_zero_based(
        inedge[:, 2:4], n_elements, "inedge element columns"
    ).reshape(n_internal, 2)
    internal_left = internal_elements[:, 0]
    internal_right = internal_elements[:, 1]

    kn = _vector(method["Kn"], n_boundary, "TPFA.Kn")
    h_esq = _vector(method["Hesq"], n_boundary, "TPFA.Hesq")
    kde = _vector(method["Kde"], n_internal, "TPFA.Kde")
    flowrate_z = _vector(
        method["flowrateZ"], n_boundary + n_internal, "TPFA.flowrateZ"
    )
    nflag = np.asarray(config["nflag"], dtype=float)
    if nflag.shape != (n_nodes, 2):
        raise ValueError(f"config.nflag must have shape ({n_nodes}, 2)")

    flowrate = np.zeros(n_boundary + n_internal, dtype=float)
    flowratedif = np.zeros_like(flowrate)

    coord1 = coord[boundary_node1]
    coord2 = coord[boundary_node2]
    edge_vector = coord1 - coord2
    edge_length = np.linalg.norm(edge_vector, axis=1)
    if np.any(edge_length == 0):
        raise ValueError("Boundary edges must have nonzero length")

    centroids = cent_elem[boundary_left]
    c1 = nflag[boundary_node1, 1]
    c2 = nflag[boundary_node2, 1]
    coefficient = kn / (h_esq * edge_length)
    term1 = np.sum((centroids - coord2) * edge_vector, axis=1) * c1
    term2 = np.sum((centroids - coord1) * -edge_vector, axis=1) * c2
    boundary_flow = -coefficient * (
        term1 + term2 - edge_length * edge_length * p[boundary_left]
    )

    numcase = int(config.get("numcase", 0))
    mobility = _face_mobility(
        config.get("visc"), numcase, n_boundary, n_internal
    )
    boundary_flow *= mobility[:n_boundary]

    boundary_flags = bedge[:, 4].astype(int)
    neumann = boundary_flags >= 200
    if np.any(neumann):
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
        gravity_boundary = flowrate_z[:n_boundary][boundary_flags > 200]
        if gravity_boundary.size != prescribed.size:
            raise ValueError("Neumann boundary flags must be greater than 200")
        boundary_flow[neumann] = (
            -edge_length[neumann] * prescribed + gravity_boundary
        )
    flowrate[:n_boundary] = boundary_flow

    internal_flow = kde * (p[internal_right] - p[internal_left])
    flowrate[n_boundary:] = mobility[n_boundary:] * internal_flow

    if numcase in (431, 435, 437, 439):
        flowrate -= flowrate_z

    flowresult = np.zeros(n_elements, dtype=float)
    flowresult += np.bincount(
        boundary_left, weights=flowrate[:n_boundary], minlength=n_elements
    )
    flowresult += np.bincount(
        internal_left, weights=flowrate[n_boundary:], minlength=n_elements
    )
    flowresult -= np.bincount(
        internal_right, weights=flowrate[n_boundary:], minlength=n_elements
    )

    if 200 < numcase < 300 or 379 < numcase < 400:
        concentration = env["conpre"]
        con = _vector(concentration["Con"], n_elements, "conpre.Con")
        kdec = _vector(concentration["Kdec"], n_internal, "conpre.Kdec")
        flowratedif[n_boundary:] = kdec * (
            con[internal_right] - con[internal_left]
        )

    return flowrate, flowresult, flowratedif, 0


def _vector(values: Any, size: int, name: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float).reshape(-1)
    if vector.size != size:
        raise ValueError(f"{name} has {vector.size} entries, expected {size}")
    return vector


def _face_mobility(
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
