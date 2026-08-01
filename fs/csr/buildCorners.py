"""CSR corner layout translated from ``+fs/+csr/buildCorners.m``."""

from typing import Any, Dict

import numpy as np


def buildCorners(FS: Dict[str, Any]) -> Dict[str, Any]:
    """Populate the flat node-to-corner layout in ``FS['csr']``."""
    mesh = FS["mesh"]
    if "esurn1" not in mesh or "esurn2" not in mesh:
        raise ValueError("FS.mesh.esurn1 and FS.mesh.esurn2 are required")

    esurn1 = np.asarray(mesh["esurn1"], dtype=int).reshape(-1)
    esurn2 = np.asarray(mesh["esurn2"], dtype=int).reshape(-1)
    n_nodes = int(mesh["nNodes"])
    node_nec = np.diff(esurn2)

    if node_nec.size != n_nodes:
        raise ValueError(
            f"diff(esurn2) has {node_nec.size} entries, expected {n_nodes}"
        )
    if esurn2[0] != 0 or esurn2[-1] != esurn1.size or np.any(node_nec < 0):
        raise ValueError("esurn2 is not a valid CSR row pointer")

    csr = FS.setdefault("csr", {})
    csr["nodePtr"] = esurn2
    csr["nCorners"] = esurn1.size
    csr["cornerElem"] = esurn1
    csr["nodeNec"] = node_nec
    csr["cornerNode"] = np.repeat(np.arange(n_nodes), node_nec)
    csr["maxNec"] = int(node_nec.max()) if node_nec.size else 0
    return FS
