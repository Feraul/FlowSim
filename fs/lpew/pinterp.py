"""Nodal interpolation translated from ``+fs/+lpew/pinterp.m`` and its legacy kernel."""

from typing import Any, Dict, Tuple, Union

import numpy as np

from fs.indexing import to_zero_based


def pinterp(
    p: np.ndarray, env: Dict[str, Any]
) -> Tuple[np.ndarray, Union[np.ndarray, float]]:
    """Interpolate element pressure, and concentration when configured, to nodes."""
    geometry = env["geometry"]
    config = env["config"]
    n_nodes = np.asarray(geometry["coord"]).shape[0]
    n_elements = np.asarray(geometry["elem"]).shape[0]
    esurn1 = to_zero_based(geometry["esurn1"], n_elements, "esurn1")
    esurn2 = np.asarray(geometry["esurn2"], dtype=int).reshape(-1)

    if esurn2.size != n_nodes + 1 or esurn2[0] != 0 or esurn2[-1] != esurn1.size:
        raise ValueError("geometry.esurn2 is not a valid node-to-element pointer")

    node_ids = np.repeat(np.arange(n_nodes), np.diff(esurn2))
    mpfad = env["premethod"]["MPFAD"]
    pressure = _interpolate(
        values=p,
        weights=mpfad["weight"],
        source=mpfad["s"],
        flags=config["nflag"],
        element_ids=esurn1,
        node_ids=node_ids,
        n_nodes=n_nodes,
        field_name="pressure",
    )

    numcase = int(config["numcase"])
    if not (200 < numcase < 300 or 379 < numcase < 400):
        return pressure, 0.0

    concentration_data = env["conpre"]
    concentration = _interpolate(
        values=concentration_data["Con"],
        weights=concentration_data["wightc"],
        source=concentration_data["sc"],
        flags=concentration_data["nflagnoc"],
        element_ids=esurn1,
        node_ids=node_ids,
        n_nodes=n_nodes,
        field_name="concentration",
    )
    return pressure, concentration


def _interpolate(
    values: np.ndarray,
    weights: np.ndarray,
    source: np.ndarray,
    flags: np.ndarray,
    element_ids: np.ndarray,
    node_ids: np.ndarray,
    n_nodes: int,
    field_name: str,
) -> np.ndarray:
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.asarray(weights, dtype=float).reshape(-1)
    source = np.asarray(source, dtype=float).reshape(-1)
    flags = np.asarray(flags, dtype=float)

    if weights.size != element_ids.size:
        raise ValueError(
            f"{field_name} weights have {weights.size} entries, "
            f"expected {element_ids.size}"
        )
    if source.size != n_nodes:
        raise ValueError(
            f"{field_name} source has {source.size} entries, expected {n_nodes}"
        )
    if flags.shape != (n_nodes, 2):
        raise ValueError(f"{field_name} flags must have shape ({n_nodes}, 2)")
    if element_ids.size and element_ids.max() >= values.size:
        raise ValueError(f"{field_name} values do not cover every referenced element")

    interpolated = np.bincount(
        node_ids,
        weights=weights * values[element_ids],
        minlength=n_nodes,
    )
    neumann = flags[:, 0] == 202
    interpolated[neumann] += source[neumann]
    dirichlet = flags[:, 0] <= 200
    interpolated[dirichlet] = flags[dirichlet, 1]
    return interpolated
