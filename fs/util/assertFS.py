"""FS invariant checker translated from ``+fs/+util/assertFS.m``."""

from typing import Any, Dict

import numpy as np


def assertFS(FS: Dict[str, Any]) -> None:
    """Raise ``ValueError`` when a core FS structural invariant is broken."""
    top_level = ("mesh", "geom", "csr", "perm", "bc", "cfg", "state", "workspace")
    missing = [name for name in top_level if name not in FS]
    if missing:
        raise ValueError(f"Missing FS top-level fields: {', '.join(missing)}")

    mesh = FS["mesh"]
    csr = FS["csr"]
    geom = FS["geom"]
    perm = FS["perm"]

    if "nNodes" in mesh and "esurn2" in mesh:
        pointer_size = np.asarray(mesh["esurn2"]).size
        expected = int(mesh["nNodes"]) + 1
        if pointer_size != expected:
            raise ValueError(
                f"FS.mesh.esurn2 has {pointer_size} entries, expected {expected}"
            )

    if "nCorners" in csr and "nodePtr" in csr:
        node_pointer = np.asarray(csr["nodePtr"], dtype=int).reshape(-1)
        if node_pointer.size:
            n_corners = int(csr["nCorners"])
            if node_pointer[0] != 0 or node_pointer[-1] != n_corners:
                raise ValueError(
                    "FS.csr.nodePtr must start at 0 and end at FS.csr.nCorners"
                )

    if "centElem" in geom and "nElems" in mesh:
        centroids = np.asarray(geom["centElem"])
        if centroids.size and centroids.shape[0] != int(mesh["nElems"]):
            raise ValueError(
                "FS.geom.centElem row count must equal FS.mesh.nElems"
            )

    if "tensor" in perm and "nElems" in mesh:
        tensor = np.asarray(perm["tensor"])
        if tensor.size and tensor.shape[0] != int(mesh["nElems"]):
            raise ValueError(
                "FS.perm.tensor row count must equal FS.mesh.nElems"
            )
