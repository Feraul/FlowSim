# -*- coding: utf-8 -*-
"""
fs.lpew.OPT — Batched geometry gather for LPEW/interpolation (translated from +fs/+lpew/OPT.m)

Outputs (numpy arrays):
- P         : (nAllNeighbors, 2 or 3) neighbor-node coords (one row per neighbor entry)
- T         : (nAllNeighbors, 2 or 3) midpoints (P + Qneighbor)/2
- O         : (nCorners, 2 or 3)      element centroids at each corner
- Qcorner   : (nCorners, 2 or 3)      owning node's coord per corner
- Qneighbor : (nAllNeighbors, 2 or 3) owning node's coord per neighbor entry

Preconditions: FS is a dict-like structure with keys:
- FS['mesh']['coord'], FS['mesh']['nsurn1'], FS['mesh']['nsurn2']
- FS['geom']['centElem']
- FS['csr']['cornerNode'], FS['csr']['cornerElem']
- FS['mesh']['nNodes']

Indices are expected to be 0-based.
"""
from typing import Tuple, Dict, Any
import numpy as np


def OPT(FS: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mesh = FS.get('mesh', {})
    geom = FS.get('geom', {})
    csr = FS.get('csr', {})

    coord = np.asarray(mesh.get('coord'))
    nsurn1 = np.asarray(mesh.get('nsurn1'))
    nsurn2 = np.asarray(mesh.get('nsurn2'))
    centElem = np.asarray(geom.get('centElem'))
    cornerNode = np.asarray(csr.get('cornerNode')).astype(int)
    cornerElem = np.asarray(csr.get('cornerElem')).astype(int)
    nNodes = int(mesh.get('nNodes', np.max(cornerNode) + 1))

    # Neighbors (P, T): one row per (node, neighbor) pair
    # P corresponds to coord(nsurn1, :)
    P = coord[nsurn1, :]

    # number of neighbors per node
    nns = np.diff(nsurn2)

    # nodeOfNeighbor: repelem((1:nNodes).', nns) in MATLAB — 0-based here
    node_indices = np.arange(nNodes, dtype=int)
    nodeOfNeighbor = np.repeat(node_indices, nns)

    Qneighbor = coord[nodeOfNeighbor, :]
    T = 0.5 * (P + Qneighbor)

    # Corners (O): one row per (node, element) pair
    O = centElem[cornerElem, :]

    # Owning-node coord per corner
    Qcorner = coord[cornerNode, :]

    return P, T, O, Qcorner, Qneighbor
