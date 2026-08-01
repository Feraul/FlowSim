"""Vectorized MPFA-D flow rate translated from ``+fs/+flow/mpfad.m``."""

from typing import Any, Dict, Tuple

import numpy as np

from fs.lpew.pinterp import pinterp as interpolate_nodes

from ._common import (
    apply_neumann,
    edge_indices,
    element_balance,
    face_mobility,
    vector,
)


def mpfad(
    p: np.ndarray,
    pinterp: np.ndarray,
    env: Dict[str, Any],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Compute MPFA-D face rates, element balance, and dispersive flow."""
    geometry = env["geometry"]
    config = env["config"]
    method = env["premethod"]["MPFAD"]

    coord = np.asarray(geometry["coord"], dtype=float)
    bedge = np.asarray(geometry["bedge"])
    inedge = np.asarray(geometry["inedge"])
    cent_elem = np.asarray(geometry["centelem"], dtype=float)
    p = np.asarray(p, dtype=float).reshape(-1)
    pinterp = np.asarray(pinterp, dtype=float).reshape(-1)

    n_nodes = coord.shape[0]
    n_elements = cent_elem.shape[0]
    n_boundary = bedge.shape[0]
    n_internal = inedge.shape[0]
    if p.size != n_elements:
        raise ValueError(f"p has {p.size} entries, expected {n_elements}")
    if pinterp.size != n_nodes:
        raise ValueError(f"pinterp has {pinterp.size} entries, expected {n_nodes}")

    edges = edge_indices(bedge, inedge, n_nodes, n_elements)
    kn = vector(method["Kn"], n_boundary, "MPFAD.Kn")
    kt = vector(method["Kt"], n_boundary, "MPFAD.Kt")
    h_esq = vector(method["Hesq"], n_boundary, "MPFAD.Hesq")
    kde = vector(method["Kde"], n_internal, "MPFAD.Kde")
    ded = vector(method["Ded"], n_internal, "MPFAD.Ded")
    flowrate_z = vector(
        method["flowrateZ"], n_boundary + n_internal, "MPFAD.flowrateZ"
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
    boundary_flow = (
        -coefficient
        * (term1 + term2 - edge_length * edge_length * p[edges.boundary_left])
        - (c2 - c1) * kt
    )

    numcase = int(config.get("numcase", 0))
    mobility = face_mobility(
        config.get("visc"), numcase, n_boundary, n_internal
    )
    boundary_flow *= mobility[:n_boundary]
    boundary_flags = bedge[:, 4].astype(int)
    apply_neumann(boundary_flow, edge_length, boundary_flags, config, flowrate_z)
    flowrate[:n_boundary] = boundary_flow

    pressure_delta = (
        p[edges.internal_right]
        - p[edges.internal_left]
        - ded
        * (
            pinterp[edges.internal_node2]
            - pinterp[edges.internal_node1]
        )
    )
    flowrate[n_boundary:] = mobility[n_boundary:] * kde * pressure_delta

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
        _, concentration_interp = interpolate_nodes(p, env)
        concentration_interp = vector(
            concentration_interp, n_nodes, "interpolated concentration"
        )
        concentration = env["conpre"]
        con = vector(concentration["Con"], n_elements, "conpre.Con")
        kdec = vector(concentration["Kdec"], n_internal, "conpre.Kdec")
        dedc = vector(concentration["Dedc"], n_internal, "conpre.Dedc")
        concentration_delta = (
            con[edges.internal_right]
            - con[edges.internal_left]
            - dedc
            * (
                concentration_interp[edges.internal_node2]
                - concentration_interp[edges.internal_node1]
            )
        )
        flowratedif[n_boundary:] = kdec * concentration_delta

    return flowrate, flowresult, flowratedif, 0
