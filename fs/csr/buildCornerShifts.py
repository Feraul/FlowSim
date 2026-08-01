"""Corner shifts translated from ``+fs/+csr/buildCornerShifts.m``."""

from typing import Any, Dict

import numpy as np


def buildCornerShifts(FS: Dict[str, Any]) -> Dict[str, Any]:
    """Precompute zero-based neighbor and corner shifts for each corner."""
    csr = FS["csr"]
    if "cornerNode" not in csr or "cornerElem" not in csr:
        raise ValueError("Run fs.csr.buildCorners before buildCornerShifts")

    mesh = FS["mesh"]
    n_nodes = int(mesh["nNodes"])
    esurn2 = np.asarray(mesh["esurn2"], dtype=int)
    nsurn2 = np.asarray(mesh["nsurn2"], dtype=int)
    corner_node = np.asarray(csr["cornerNode"], dtype=int)
    n_corners = int(csr["nCorners"])

    nec = np.diff(esurn2)
    nns = np.diff(nsurn2)
    if nec.size != n_nodes or nns.size != n_nodes:
        raise ValueError("esurn2 and nsurn2 must contain one row per node")
    if np.any((nns != nec) & (nns != nec + 1)):
        raise ValueError("Each node must have nec neighbors, or nec + 1 on a boundary")

    corner_local = np.arange(n_corners) - np.repeat(esurn2[:-1], nec)
    corner_is_first = corner_local == 0
    corner_is_last = corner_local == nec[corner_node] - 1
    node_is_interior = nns == nec
    interior_corner = node_is_interior[corner_node]

    neighbor_base = nsurn2[corner_node]
    t_current = neighbor_base + corner_local
    t_next = t_current + 1
    wrap_mask = corner_is_last & interior_corner
    t_next[wrap_mask] = neighbor_base[wrap_mask]

    corner_index = np.arange(n_corners)
    corner_base = esurn2[corner_node]

    corner_prev = corner_index - 1
    interior_first = corner_is_first & interior_corner
    corner_prev[interior_first] = (
        corner_base[interior_first] + nec[corner_node[interior_first]] - 1
    )
    boundary_first = corner_is_first & ~interior_corner
    corner_prev[boundary_first] = corner_index[boundary_first]

    corner_next = corner_index + 1
    corner_next[wrap_mask] = corner_base[wrap_mask]
    boundary_last = corner_is_last & ~interior_corner
    corner_next[boundary_last] = corner_index[boundary_last]

    csr["cornerLocal"] = corner_local
    csr["cornerIsFirst"] = corner_is_first
    csr["cornerIsLast"] = corner_is_last
    csr["nodeIsInterior"] = node_is_interior
    csr["tCurrent"] = t_current
    csr["tNext"] = t_next
    csr["cornerPrev"] = corner_prev
    csr["cornerNext"] = corner_next
    csr["boundaryFirstMask"] = boundary_first
    csr["boundaryLastMask"] = boundary_last
    return FS
