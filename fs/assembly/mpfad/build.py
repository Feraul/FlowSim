"""Vectorized MPFA-D assembly translated from ``+fs/+assembly/+mpfad/build.m``."""

from typing import Any, Callable, Dict, Tuple

import numpy as np
from scipy import sparse

from fs.flow._common import edge_indices, vector
from fs.indexing import to_zero_based


def build(
    env: Dict[str, Any], parms: Dict[str, Any]
) -> Tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Assemble the MPFA-D sparse matrix and right-hand-side vector."""
    geometry = env["geometry"]
    config = env["config"]
    method = env["premethod"]["MPFAD"]

    coord = np.asarray(geometry["coord"], dtype=float)
    elem = np.asarray(geometry["elem"])
    bedge = np.asarray(geometry["bedge"])
    inedge = np.asarray(geometry["inedge"])
    cent_elem = np.asarray(geometry["centelem"], dtype=float)
    esurn1 = to_zero_based(geometry["esurn1"], elem.shape[0], "esurn1")
    esurn2 = np.asarray(geometry["esurn2"], dtype=int).reshape(-1)

    n_nodes = coord.shape[0]
    n_elements = elem.shape[0]
    n_boundary = bedge.shape[0]
    n_internal = inedge.shape[0]
    if cent_elem.shape[0] != n_elements:
        raise ValueError("geometry.centelem must contain one row per element")
    if (
        esurn2.size != n_nodes + 1
        or esurn2[0] != 0
        or esurn2[-1] != esurn1.size
        or np.any(np.diff(esurn2) < 0)
    ):
        raise ValueError("geometry.esurn2 is not a valid node-to-element pointer")

    edges = edge_indices(bedge, inedge, n_nodes, n_elements)
    ded = vector(method["Ded"], n_internal, "MPFAD.Ded")
    kde = vector(method["Kde"], n_internal, "MPFAD.Kde")
    kn = vector(method["Kn"], n_boundary, "MPFAD.Kn")
    kt = vector(method["Kt"], n_boundary, "MPFAD.Kt")
    h_esq = vector(method["Hesq"], n_boundary, "MPFAD.Hesq")
    flowrate_z = vector(
        method["flowrateZ"], n_boundary + n_internal, "MPFAD.flowrateZ"
    )
    flowresult_z = vector(
        method["flowresultZ"], n_elements, "MPFAD.flowresultZ"
    )
    nodal_source = vector(method["s"], n_nodes, "MPFAD.s")
    weight = vector(method["weight"], esurn1.size, "MPFAD.weight")

    nflag = np.asarray(config["nflag"], dtype=float)
    if nflag.shape != (n_nodes, 2):
        raise ValueError(f"config.nflag must have shape ({n_nodes}, 2)")

    coord1 = coord[edges.boundary_node1]
    coord2 = coord[edges.boundary_node2]
    v0 = coord2 - coord1
    v1 = cent_elem[edges.boundary_left] - coord1
    v2 = cent_elem[edges.boundary_left] - coord2
    edge_length = np.linalg.norm(v0, axis=1)
    if np.any(edge_length == 0):
        raise ValueError("Boundary edges must have nonzero length")

    flags = bedge[:, 4].astype(int)
    is_dirichlet = flags < 200
    is_neumann = ~is_dirichlet
    dir_left = edges.boundary_left[is_dirichlet]

    v0_dir = v0[is_dirichlet]
    coefficient = -kn[is_dirichlet] / (
        h_esq[is_dirichlet] * edge_length[is_dirichlet]
    )
    matrix_dir = -coefficient * np.sum(v0_dir * v0_dir, axis=1)
    c1 = nflag[edges.boundary_node1[is_dirichlet], 1]
    c2 = nflag[edges.boundary_node2[is_dirichlet], 1]
    dot_v2_v0 = np.sum(v2[is_dirichlet] * -v0_dir, axis=1)
    dot_v1_v0 = np.sum(v1[is_dirichlet] * v0_dir, axis=1)
    rhs_dir = (
        -coefficient * (dot_v2_v0 * c1 + dot_v1_v0 * c2)
        + (c2 - c1) * kt[is_dirichlet]
    )

    neumann_callback = _benchmark_callback(env, "calcularNeumannBoundary")
    rhs_neumann = vector(
        neumann_callback(
            is_neumann,
            bedge,
            config["bcflag"],
            config["nflagface"],
            flowrate_z,
            edge_length,
            geometry["normals"],
            env,
        ),
        np.count_nonzero(is_neumann),
        "Neumann callback result",
    )
    neumann_left = edges.boundary_left[is_neumann]

    rhs_left = np.zeros(n_elements)
    rhs_right = np.zeros(n_elements)
    node1_flags = nflag[edges.internal_node1, 0]
    node2_flags = nflag[edges.internal_node2, 0]

    dirichlet1 = node1_flags < 200
    dirichlet2 = node2_flags < 200
    correction1 = (
        kde[dirichlet1]
        * ded[dirichlet1]
        * nflag[edges.internal_node1[dirichlet1], 1]
    )
    correction2 = (
        kde[dirichlet2]
        * ded[dirichlet2]
        * nflag[edges.internal_node2[dirichlet2], 1]
    )
    rhs_left += np.bincount(
        edges.internal_left[dirichlet1],
        weights=-correction1,
        minlength=n_elements,
    )
    rhs_right += np.bincount(
        edges.internal_right[dirichlet1],
        weights=correction1,
        minlength=n_elements,
    )
    rhs_left += np.bincount(
        edges.internal_left[dirichlet2],
        weights=correction2,
        minlength=n_elements,
    )
    rhs_right += np.bincount(
        edges.internal_right[dirichlet2],
        weights=-correction2,
        minlength=n_elements,
    )

    neumann1 = np.isin(node1_flags, (201, 202))
    neumann2 = np.isin(node2_flags, (201, 202))
    correction1 = (
        kde[neumann1]
        * ded[neumann1]
        * nodal_source[edges.internal_node1[neumann1]]
    )
    correction2 = (
        kde[neumann2]
        * ded[neumann2]
        * nodal_source[edges.internal_node2[neumann2]]
    )
    rhs_left += np.bincount(
        edges.internal_left[neumann1],
        weights=-correction1,
        minlength=n_elements,
    )
    rhs_right += np.bincount(
        edges.internal_right[neumann1],
        weights=correction1,
        minlength=n_elements,
    )
    rhs_left += np.bincount(
        edges.internal_left[neumann2],
        weights=correction2,
        minlength=n_elements,
    )
    rhs_right += np.bincount(
        edges.internal_right[neumann2],
        weights=-correction2,
        minlength=n_elements,
    )

    scatter1 = _weight_scatter(
        np.flatnonzero(node1_flags > 200),
        edges.internal_node1,
        edges.internal_left,
        edges.internal_right,
        esurn1,
        esurn2,
        weight,
        kde,
        ded,
        1.0,
    )
    scatter2 = _weight_scatter(
        np.flatnonzero(node2_flags > 200),
        edges.internal_node2,
        edges.internal_left,
        edges.internal_right,
        esurn1,
        esurn2,
        weight,
        kde,
        ded,
        -1.0,
    )

    rows = np.concatenate(
        (
            dir_left,
            edges.internal_left,
            edges.internal_left,
            edges.internal_right,
            edges.internal_right,
            scatter1[0],
            scatter2[0],
        )
    )
    columns = np.concatenate(
        (
            dir_left,
            edges.internal_left,
            edges.internal_right,
            edges.internal_right,
            edges.internal_left,
            scatter1[1],
            scatter2[1],
        )
    )
    values = np.concatenate(
        (matrix_dir, -kde, kde, -kde, kde, scatter1[2], scatter2[2])
    )
    matrix = sparse.coo_matrix(
        (values, (rows, columns)), shape=(n_elements, n_elements)
    ).tocsr()

    rhs = np.bincount(dir_left, weights=rhs_dir, minlength=n_elements)
    rhs += np.bincount(
        neumann_left, weights=rhs_neumann, minlength=n_elements
    )
    rhs += rhs_left + rhs_right

    temporal_callback = _benchmark_callback(env, "adicionarTermoTemporal")
    matrix, rhs = temporal_callback(
        matrix, rhs, parms, flowresult_z, env
    )
    matrix = sparse.csr_matrix(matrix)
    rhs = vector(rhs, n_elements, "temporal right-hand side")

    elembedge = np.empty((0, 2), dtype=float)
    if str(config.get("modflowcase", "n")).lower() == "y":
        nflagface = np.asarray(config["nflagface"], dtype=float)
        if nflagface.shape != (n_boundary, 2):
            raise ValueError(
                f"config.nflagface must have shape ({n_boundary}, 2)"
            )
        elembedge = np.column_stack(
            (dir_left, nflagface[is_dirichlet, 1])
        )
        editable = matrix.tolil()
        for element, value in elembedge:
            element = int(element)
            editable.rows[element] = [element]
            editable.data[element] = [1.0]
            rhs[element] = value
        matrix = editable.tocsr()

    return matrix, rhs, elembedge


def _weight_scatter(
    selected_edges: np.ndarray,
    edge_nodes: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    esurn1: np.ndarray,
    esurn2: np.ndarray,
    weight: np.ndarray,
    kde: np.ndarray,
    ded: np.ndarray,
    sign: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if selected_edges.size == 0:
        empty = np.empty(0, dtype=float)
        return empty.astype(int), empty.astype(int), empty

    nodes = edge_nodes[selected_edges]
    counts = np.diff(esurn2)[nodes]
    total = int(np.sum(counts))
    if total == 0:
        empty = np.empty(0, dtype=float)
        return empty.astype(int), empty.astype(int), empty

    starts = esurn2[nodes]
    repeated_starts = np.repeat(starts, counts)
    group_offsets = np.repeat(
        np.cumsum(np.concatenate(([0], counts[:-1]))), counts
    )
    corner_indices = (
        repeated_starts + np.arange(total) - group_offsets
    )
    columns = esurn1[corner_indices]
    coefficients = np.repeat(kde[selected_edges] * ded[selected_edges], counts)
    weighted = sign * coefficients * weight[corner_indices]
    rows_left = np.repeat(left[selected_edges], counts)
    rows_right = np.repeat(right[selected_edges], counts)
    return (
        np.concatenate((rows_left, rows_right)),
        np.concatenate((columns, columns)),
        np.concatenate((weighted, -weighted)),
    )


def _benchmark_callback(
    env: Dict[str, Any], name: str
) -> Callable[..., Any]:
    benchmark = env.get("benchmark")
    if isinstance(benchmark, dict):
        callback = benchmark.get(name)
    else:
        callback = getattr(benchmark, name, None)
    if not callable(callback):
        raise TypeError(f"env.benchmark.{name} must be callable")
    return callback
