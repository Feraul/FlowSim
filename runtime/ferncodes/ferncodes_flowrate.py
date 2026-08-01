# -*- coding: utf-8 -*-
"""
runtime.ferncodes.ferncodes_flowrate — tradução de legacy/ferncodes/mpfad/ferncodes_flowrate.m

Assunções:
- env é um dict com chaves geometry, config, premethod, benchmark.
- Todas as matrizes de entrada já estão em 0-based indices (se vierem do MATLAB, subtrair 1).
- p: array de pressões por elemento (nelem,)
- pinterp: array de pressões interpoladas por nó (nNodes,)

Retorna:
- flowrate: vetor (nb+ni,)
- flowresult: vetor (nelem,)
- flowratedif: vetor (nb+ni,)
- faceaux: int (compatibilidade)
"""
from typing import Tuple, Any
import numpy as np


def ferncodes_flowrate(p: np.ndarray, pinterp: np.ndarray, env: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Any]:
    coord = np.asarray(env["geometry"]["coord"])            # (nNodes, 2)
    bedge = np.asarray(env["geometry"]["bedge"])            # (nb, >=5)
    inedge = np.asarray(env["geometry"]["inedge"])          # (ni, 4)
    centelem = np.asarray(env["geometry"]["centelem"])      # (nelem, 2)
    bcflag = np.asarray(env["config"]["bcflag"])
    numcase = env["config"].get("numcase", 0)
    viscosity = env["config"].get("visc", None)
    flowrateZ = env["premethod"]["MPFAD"].get("flowrateZ", None)

    Kde = np.asarray(env["premethod"]["MPFAD"].get("Kde"))
    Ded = np.asarray(env["premethod"]["MPFAD"].get("Ded"))
    Kn = np.asarray(env["premethod"]["MPFAD"].get("Kn"))
    Kt = np.asarray(env["premethod"]["MPFAD"].get("Kt"))
    Hesq = np.asarray(env["premethod"]["MPFAD"].get("Hesq"))
    nflag = np.asarray(env["config"].get("nflag"))

    bedgesize = bedge.shape[0]
    inedgesize = inedge.shape[0]

    flowrate = np.zeros(bedgesize + inedgesize, dtype=float)
    flowratedif = np.zeros(bedgesize + inedgesize, dtype=float)
    flowresult = np.zeros(centelem.shape[0], dtype=float)

    # === BOUNDARY EDGES ===
    B1 = bedge[:, 0].astype(int)
    B2 = bedge[:, 1].astype(int)
    lef = bedge[:, 2].astype(int)

    coordB1 = coord[B1, :]
    coordB2 = coord[B2, :]

    edgevec = coordB1 - coordB2
    nor = np.linalg.norm(edgevec, axis=1)

    O = centelem[lef, :]

    dirmask = bedge[:, 4] < 200
    neumask = ~dirmask

    c1 = nflag[B1, 1]
    c2 = nflag[B2, 1]

    A = Kn / (Hesq * nor)

    term1 = np.sum((O - coordB2) * (coordB1 - coordB2), axis=1) * c1
    term2 = np.sum((O - coordB1) * (coordB2 - coordB1), axis=1) * c2

    p_lef = p[lef]
    flowrate_b = -A * (term1 + term2 - (nor ** 2) * p_lef) - (c2 - c1) * Kt

    # viscosity handling (mirrors MATLAB conditions)
    visonface = np.ones(bedgesize, dtype=float)
    if 30 < numcase < 200 and viscosity is not None:
        # viscosity assumed shaped (bedgesize, ...)
        visonface = np.sum(viscosity[0:bedgesize, :], axis=1)
    elif 200 < numcase < 300 and viscosity is not None:
        if numcase in (246, 245, 247, 248, 249, 251):
            visonface = viscosity[0:bedgesize, :].reshape(bedgesize)

    flowrate_b = visonface * flowrate_b

    # Neumann boundary
    if np.any(neumask):
        bcvals = np.zeros(bedgesize, dtype=float)
        # find location of flags
        flags = bedge[neumask, 4]
        # match flags to bcflag[:,0]
        # build mapping from flag id to value
        flag_ids = bcflag[:, 0]
        # for flags that match, get value
        # use searchsorted-like approach: build dict
        flag_to_val = {int(fid): float(val) for fid, val in bcflag[:, :2]}
        for i, flag in enumerate(flags):
            idx_global = np.nonzero(neumask)[0][i]
            bcvals[idx_global] = flag_to_val.get(int(flag), 0.0)
        mask222 = bedge[:, 4] > 200
        # collect flowrateZ entries where mask222 true
        if flowrateZ is not None:
            z_vals = np.asarray(flowrateZ).reshape(-1)
            z_selected = z_vals[mask222.astype(bool)]
            # assign to neumask positions — assume same length
            # if lengths mismatch, pad or trim accordingly
            n_neu = np.count_nonzero(neumask)
            if z_selected.size == n_neu:
                flowrate_b[neumask] = -nor[neumask] * bcvals[neumask] + z_selected
            else:
                # fallback: use first n_neu values or zeros
                z_use = z_selected[:n_neu] if z_selected.size >= n_neu else np.pad(z_selected, (0, n_neu - z_selected.size))
                flowrate_b[neumask] = -nor[neumask] * bcvals[neumask] + z_use
        else:
            flowrate_b[neumask] = -nor[neumask] * bcvals[neumask]

    flowrate[:bedgesize] = flowrate_b

    # === INTERNAL EDGES ===
    node1 = inedge[:, 0].astype(int)
    node2 = inedge[:, 1].astype(int)
    lef_i = inedge[:, 2].astype(int)
    rel = inedge[:, 3].astype(int)

    p1 = pinterp[node1]
    p2 = pinterp[node2]

    visonface_i = np.ones(inedgesize, dtype=float)
    if 30 < numcase < 200 and viscosity is not None:
        visonface_i = np.sum(viscosity[bedgesize:bedgesize + inedgesize, :], axis=1)
    elif 200 < numcase < 300 and viscosity is not None:
        if numcase in (246, 245, 247, 248, 249, 251):
            visonface_i = viscosity[bedgesize:bedgesize + inedgesize, :].reshape(inedgesize)

    flowrate[bedgesize:] = visonface_i * Kde * (p[rel] - p[lef_i] - Ded * (p2 - p1))

    if numcase in (435, 431, 437, 439):
        # subtract gravitic contribution
        if flowrateZ is not None:
            flowrate = flowrate - np.asarray(flowrateZ).reshape(-1)

    # === FLOWRESULT ===
    # accumarray(auxlef, flowrate[0:bedgesize])
    auxlef = bedge[:, 2].astype(int)
    flowresult += np.bincount(auxlef, weights=flowrate[:bedgesize], minlength=flowresult.size)

    idx = np.arange(bedgesize, bedgesize + inedgesize)
    flowresult += np.bincount(lef_i, weights=flowrate[idx], minlength=flowresult.size)
    flowresult -= np.bincount(rel, weights=flowrate[idx], minlength=flowresult.size)

    # === DISPERSIVE FLOW (placeholder) ===
    if (200 < numcase < 300) or (379 < numcase < 400):
        # The MATLAB code references variables (cinterp, Kdec, Con, Dedc) that
        # are part of concentration-coupled logic. Implement when translating
        # concentration modules. For now leave zeros.
        pass

    faceaux = 0
    return flowrate, flowresult, flowratedif, faceaux
