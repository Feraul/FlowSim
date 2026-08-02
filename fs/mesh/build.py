"""Legacy environment adapter translated from ``+fs/+mesh/build.m``."""

from typing import Any, Dict

import numpy as np


def build(env: Dict[str, Any]) -> Dict[str, Any]:
    """Create the vectorization-friendly FS structure without mutating ``env``."""
    geometry = env["geometry"]
    config = env.get("config", {})

    coord = np.asarray(geometry["coord"])
    elem = np.asarray(geometry["elem"])
    bedge = np.asarray(geometry["bedge"])
    inedge = np.asarray(geometry["inedge"])

    mesh = {
        "nNodes": coord.shape[0],
        "nElems": elem.shape[0],
        "nBFaces": bedge.shape[0],
        "nIFaces": inedge.shape[0],
        "coord": coord,
        "elem": elem,
        "bedge": bedge,
        "inedge": inedge,
        "nsurn1": _zero_based(geometry["nsurn1"], coord.shape[0], "nsurn1"),
        "nsurn2": np.asarray(geometry["nsurn2"], dtype=int).reshape(-1),
        "esurn1": _zero_based(geometry["esurn1"], elem.shape[0], "esurn1"),
        "esurn2": np.asarray(geometry["esurn2"], dtype=int).reshape(-1),
    }

    geom = {
        "centElem": np.asarray(geometry["centelem"]),
        "elemArea": np.asarray(geometry["elemarea"]),
    }
    if "normals" in geometry:
        geom["normalInt"] = np.asarray(geometry["normals"])

    perm = {}
    if "perm" in config:
        perm["tensor"] = np.asarray(config["perm"])
    if "kmap" in config:
        perm["kmap"] = np.asarray(config["kmap"])

    bc = {}
    if "nflag" in config:
        bc["nflag"] = np.asarray(config["nflag"])
    if "nflagface" in config:
        bc["nflagFace"] = np.asarray(config["nflagface"])
    if "bcflag" in config:
        bc["bcflag"] = np.asarray(config["bcflag"])

    return {
        "mesh": mesh,
        "geom": geom,
        "csr": {},
        "perm": perm,
        "bc": bc,
        "cfg": dict(config),
        "state": {},
        "workspace": {},
    }


def _zero_based(values: Any, size: int, name: str) -> np.ndarray:
    indices = np.asarray(values, dtype=int).reshape(-1)
    if indices.size == 0:
        return indices
    if np.any(indices < 0):
        raise ValueError(f"{name} contains a negative index")
    if np.any(indices == 0):
        if np.any(indices >= size):
            raise ValueError(f"{name} contains an index outside 0..{size - 1}")
        return indices
    if np.any(indices > size):
        raise ValueError(f"{name} contains an index outside 1..{size}")
    return indices - 1
