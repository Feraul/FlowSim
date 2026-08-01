"""Vectorized TPFA flow rate translated from ``+fs/+flow/tpfa.m``."""

from typing import Any, Dict, Tuple

import numpy as np

from ._common import (
    apply_neumann,
    edge_indices,
    element_balance,
    face_mobility,
    vector,
)


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

    edges = edge_indices(bedge, inedge, n_nodes, n_elements)

    kn = vector(method["Kn"], n_boundary, "TPFA.Kn")
    h_esq = vector(method["Hesq"], n_boundary, "TPFA.Hesq")
    kde = vector(method["Kde"], n_internal, "TPFA.Kde")
    flowrate_z = vector(
        method["flowrateZ"], n_boundary + n_internal, "TPFA.flowrateZ"
    )
    nflag = np.asarray(config["nflag"], dtype=float)
    if nflag.shape != (n_nodes, 2):
        raise ValueError(f"config.nflag must have shape ({n_nodes}, 2)")

    flowrate = np.zeros(n_boundary + n_internal, dtype=float)
    flowratedif = np.zeros_like(flowrate)

    coord1 = coord[edges.boundary_node1]
    coord2 = coord[edges.boundary_node2]
    edge_vector = coord1 - coord2
    edge_length = np.linalg.norm(edge_vector, axis=1)
    if np.any(edge_length == 0):
        raise ValueError("Boundary edges must have nonzero length")

    centroids = cent_elem[edges.boundary_left]
    c1 = nflag[edges.boundary_node1, 1]
    c2 = nflag[edges.boundary_node2, 1]
    coefficient = kn / (h_esq * edge_length)
    term1 = np.sum((centroids - coord2) * edge_vector, axis=1) * c1
    term2 = np.sum((centroids - coord1) * -edge_vector, axis=1) * c2
    boundary_flow = -coefficient * (
        term1 + term2 - edge_length * edge_length * p[edges.boundary_left]
    )

    numcase = int(config.get("numcase", 0))
    mobility = face_mobility(
        config.get("visc"), numcase, n_boundary, n_internal
    )
    boundary_flow *= mobility[:n_boundary]

    boundary_flags = bedge[:, 4].astype(int)
    apply_neumann(boundary_flow, edge_length, boundary_flags, config, flowrate_z)
    flowrate[:n_boundary] = boundary_flow

    internal_flow = kde * (p[edges.internal_right] - p[edges.internal_left])
    flowrate[n_boundary:] = mobility[n_boundary:] * internal_flow

    if numcase in (431, 435, 437, 439):
        flowrate -= flowrate_z

    flowresult = element_balance(
        flowrate,
        edges.boundary_left,
        edges.internal_left,
        edges.internal_right,
        n_elements,
    )

    if 200 < numcase < 300 or 379 < numcase < 400:
        concentration = env["conpre"]
        con = vector(concentration["Con"], n_elements, "conpre.Con")
        kdec = vector(concentration["Kdec"], n_internal, "conpre.Kdec")
        flowratedif[n_boundary:] = kdec * (
            con[edges.internal_right] - con[edges.internal_left]
        )

    return flowrate, flowresult, flowratedif, 0
