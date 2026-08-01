"""Vectorized TPFA assembly translated from ``+fs/+assembly/+tpfa/build.m``."""

from typing import Any, Callable, Dict, Tuple

import numpy as np
from scipy import sparse

from fs.flow._common import edge_indices, vector


def build(
    env: Dict[str, Any], parms: Dict[str, Any]
) -> Tuple[sparse.csr_matrix, np.ndarray, np.ndarray]:
    """Assemble the TPFA sparse matrix and right-hand-side vector."""
    geometry = env["geometry"]
    config = env["config"]
    method = env["premethod"]["TPFA"]

    coord = np.asarray(geometry["coord"], dtype=float)
    elem = np.asarray(geometry["elem"])
    bedge = np.asarray(geometry["bedge"])
    inedge = np.asarray(geometry["inedge"])
    cent_elem = np.asarray(geometry["centelem"], dtype=float)

    n_nodes = coord.shape[0]
    n_elements = elem.shape[0]
    n_boundary = bedge.shape[0]
    n_internal = inedge.shape[0]
    if cent_elem.shape[0] != n_elements:
        raise ValueError("geometry.centelem must contain one row per element")

    edges = edge_indices(bedge, inedge, n_nodes, n_elements)
    kn = vector(method["Kn"], n_boundary, "TPFA.Kn")
    h_esq = vector(method["Hesq"], n_boundary, "TPFA.Hesq")
    kde = vector(method["Kde"], n_internal, "TPFA.Kde")
    flowrate_z = vector(
        method["flowrateZ"], n_boundary + n_internal, "TPFA.flowrateZ"
    )
    flowresult_z = vector(
        method["flowresultZ"], n_elements, "TPFA.flowresultZ"
    )

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
    rhs_dir = -coefficient * (dot_v2_v0 * c1 + dot_v1_v0 * c2)

    rhs_neumann = _neumann_values(
        env, is_neumann, flags, edge_length, flowrate_z
    )
    neumann_left = edges.boundary_left[is_neumann]

    rows = np.concatenate(
        (
            dir_left,
            edges.internal_left,
            edges.internal_left,
            edges.internal_right,
            edges.internal_right,
        )
    )
    columns = np.concatenate(
        (
            dir_left,
            edges.internal_left,
            edges.internal_right,
            edges.internal_right,
            edges.internal_left,
        )
    )
    values = np.concatenate((matrix_dir, -kde, kde, -kde, kde))
    matrix = sparse.coo_matrix(
        (values, (rows, columns)), shape=(n_elements, n_elements)
    ).tocsr()

    rhs = np.bincount(dir_left, weights=rhs_dir, minlength=n_elements)
    rhs += np.bincount(
        neumann_left, weights=rhs_neumann, minlength=n_elements
    )

    numcase = config.get("numcase", 0)
    if 330 <= numcase < 400 or 400 < numcase < 500:
        callback = _benchmark_callback(env, "adicionarTermoTemporal")
        matrix, rhs = callback(matrix, rhs, parms, flowresult_z, env)
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


def _neumann_values(
    env: Dict[str, Any],
    is_neumann: np.ndarray,
    flags: np.ndarray,
    edge_length: np.ndarray,
    flowrate_z: np.ndarray,
) -> np.ndarray:
    count = np.count_nonzero(is_neumann)
    if count == 0:
        return np.empty(0)

    config = env["config"]
    numcase = config.get("numcase", 0)
    if numcase in (341, 341.1):
        callback = _benchmark_callback(env, "calcularNeumannBoundary")
        values = callback(
            is_neumann,
            env["geometry"]["bedge"],
            config["bcflag"],
            config["nflagface"],
            flowrate_z,
            edge_length,
            env["geometry"]["normals"],
            env,
        )
        return vector(values, count, "Neumann callback result")

    bcflag = np.asarray(config["bcflag"], dtype=float)
    value_by_flag = {int(flag): value for flag, value in bcflag[:, :2]}
    neumann_flags = flags[is_neumann]
    missing = sorted(set(neumann_flags) - value_by_flag.keys())
    if missing:
        raise ValueError(f"Neumann flags missing from config.bcflag: {missing}")
    prescribed = np.fromiter(
        (value_by_flag[flag] for flag in neumann_flags),
        dtype=float,
        count=count,
    )
    gravity = flowrate_z[: flags.size][flags > 200]
    if gravity.size != count:
        raise ValueError("Neumann boundary flags must be greater than 200")
    return edge_length[is_neumann] * prescribed + gravity


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
