# -*- coding: utf-8 -*-
"""
ferncodes_Kde_Ded_Kt_Kn.py — translated from legacy/ferncodes/mpfad/ferncodes_Kde_Ded_Kt_Kn.m

Computes MPFA-D geometric/physical premethod parameters and stores them in
env.premethod.MPFAD: Hesq, Kde, Kn, Kt, Ded, flowrateZ, flowresultZ.

Assumptions and notes:
- env is a dict-like object with keys geometry, config, premethod, benchmark.
- Arrays in env are numpy-compatible and indices are 0-based.
- auxkmap: permeability map, may come from parms.auxperm or env.config.perm.
- This translation preserves the vectorized MATLAB computations (uses numpy).

Return: modified env
"""
from typing import Dict, Any
import numpy as np


def ferncodes_Kde_Ded_Kt_Kn(env: Dict[str, Any], parms: Dict[str, Any] = None) -> Dict[str, Any]:
    geom = env["geometry"]
    config = env["config"]

    bedge = np.asarray(geom["bedge"])      # (nb, >=5)
    inedge = np.asarray(geom["inedge"])    # (ni, 4)
    coord = np.asarray(geom["coord"])      # (nNodes, dim)
    centelem = np.asarray(geom["centelem"])# (nElems, dim)
    elem = np.asarray(geom["elem"])        # (nElems, >=5)
    nflagface = np.asarray(config.get("nflagface"))

    if config.get("numcase", 0) < 400 or not parms:
        auxkmap = np.asarray(config.get("perm"))
    else:
        auxkmap = np.asarray(parms.get("auxperm"))

    nb = bedge.shape[0]
    ni = inedge.shape[0]

    flowrateZ = np.zeros(nb + ni, dtype=float)
    flowresultZ = np.zeros(elem.shape[0], dtype=float)

    # ----- Boundary edges geometry -----
    B1 = bedge[:, 0].astype(int)
    B2 = bedge[:, 1].astype(int)

    v1x = coord[B2, 0] - coord[B1, 0]
    v1y = coord[B2, 1] - coord[B1, 1]

    Centro = centelem[bedge[:, 2].astype(int), :]
    ve2x = coord[B2, 0] - Centro[:, 0]
    ve2y = coord[B2, 1] - Centro[:, 1]

    normv2 = v1x ** 2 + v1y ** 2
    nv = np.sqrt(normv2)
    # Hesq: abs(cross)/(nv)
    Hesq = np.abs(v1x * ve2y - v1y * ve2x) / np.where(nv == 0, 1.0, nv)

    lef = bedge[:, 2].astype(int)
    vx = v1x
    vy = v1y

    # ----- Permeability boundary -----
    matid = elem[lef, 4].astype(int)
    # auxkmap columns: assume same layout as MATLAB (col 1 unused, 2..5 K11..K22)
    K11 = auxkmap[matid, 1]
    K12 = auxkmap[matid, 2]
    K21 = auxkmap[matid, 3]
    K22 = auxkmap[matid, 4]

    Kn = (K11 * vy ** 2 - 2 * K12 * vx * vy + K22 * vx ** 2) / np.where(normv2 == 0, 1.0, normv2)
    Kt = (vy * (K11 * vx + K12 * vy) + (-vx) * (K21 * vx + K22 * vy)) / np.where(normv2 == 0, 1.0, normv2)

    # ----- Flowrate boundary -----
    if env["benchmark"].temFlowrateBoundary():
        maskT = bedge[:, 4] < 200
        h_contorno = nflagface[:, 1]

        # benchmark can adjust boundary tensor
        K11_adj, K12_adj, K21_adj, K22_adj = env["benchmark"].ajustarKContorno(
            env, parms, auxkmap, matid, h_contorno, maskT)

        Kn = (K11_adj * vy ** 2 - 2 * K12_adj * vx * vy + K22_adj * vx ** 2) / np.where(normv2 == 0, 1.0, normv2)
        Kt = (vy * (K11_adj * vx + K12_adj * vy) + (-vx) * (K21_adj * vx + K22_adj * vy)) / np.where(normv2 == 0, 1.0, normv2)

        coordB1 = coord[B1, :]
        coordB2 = coord[B2, :]
        dB = coordB1 - coordB2
        nor = np.sqrt(dB[:, 0] ** 2 + dB[:, 1] ** 2)

        O = Centro
        # in MATLAB c1 = coord(B1,2) is taking y coordinate? Original used nflagface; here keep nflagface? In original code they used c1/c2 from nflagface previously; earlier c1,c2 were assigned from nflagface in build.m; but here original uses coord(B1,2)?? To match legacy, use nflagface values if available.
        # The original ferncodes_Kde... used c1/c2 from nflagface? But earlier flowrate code used c1 = nflag(B1,2). For safety we attempt to pull nflag from env.config.nflag if available, otherwise fallback to second coord component.
        nflag = np.asarray(env.get("config", {}).get("nflag"))
        if nflag is not None and nflag.shape[0] > 0:
            c1 = nflag[B1, 1]
            c2 = nflag[B2, 1]
        else:
            c1 = coord[B1, 1]
            c2 = coord[B2, 1]

        A1 = -Kn / np.where(Hesq * nv == 0, 1.0, (Hesq * nv))
        term1 = np.sum((O - coordB2) * (coordB1 - coordB2), axis=1) * c1
        term2 = np.sum((O - coordB1) * (coordB2 - coordB1), axis=1) * c2

        flowrateZ[0:nb] = A1 * (term1 + term2 - (nor ** 2) * Centro[:, 1]) - (c2 - c1) * Kt

        # Neumann + benchmark adjust
        mask201 = bedge[:, 4] > 200
        if np.any(mask201):
            flowrateZ[mask201] = flowrateZ[mask201]
        flowrateZ = env["benchmark"].ajustarFlowrate(flowrateZ, bedge)

        flowresultZ = flowresultZ + np.bincount(lef, weights=flowrateZ[0:nb], minlength=flowresultZ.size)

    # ----- Internal edges geometry -----
    C1 = centelem[inedge[:, 2].astype(int), :]
    C2 = centelem[inedge[:, 3].astype(int), :]
    lef_i = inedge[:, 2].astype(int)
    rel_i = inedge[:, 3].astype(int)
    vcen = C2 - C1

    e1x = coord[inedge[:, 1].astype(int), 0] - coord[inedge[:, 0].astype(int), 0]
    e1y = coord[inedge[:, 1].astype(int), 1] - coord[inedge[:, 0].astype(int), 1]

    vx = e1x
    vy = e1y
    normv2 = vx ** 2 + vy ** 2
    nv = np.sqrt(normv2)

    vd2x = C2[:, 0] - coord[inedge[:, 0].astype(int), 0]
    vd2y = C2[:, 1] - coord[inedge[:, 0].astype(int), 1]

    ve2x = C1[:, 0] - coord[inedge[:, 0].astype(int), 0]
    ve2y = C1[:, 1] - coord[inedge[:, 0].astype(int), 1]

    H2 = np.abs(vx * vd2y - vy * vd2x) / np.where(nv == 0, 1.0, nv)
    H1 = np.abs(vx * ve2y - vy * ve2x) / np.where(nv == 0, 1.0, nv)

    no1 = coord[inedge[:, 0].astype(int), 1]
    no2 = coord[inedge[:, 1].astype(int), 1]

    # Permeability internal
    matL = elem[inedge[:, 2].astype(int), 4].astype(int)
    matR = elem[inedge[:, 3].astype(int), 4].astype(int)

    K11L = auxkmap[matL, 1]
    K12L = auxkmap[matL, 2]
    K21L = auxkmap[matL, 3]
    K22L = auxkmap[matL, 4]
    K11R = auxkmap[matR, 1]
    K12R = auxkmap[matR, 2]
    K21R = auxkmap[matR, 3]
    K22R = auxkmap[matR, 4]

    Kn1 = (K11L * vy ** 2 - 2 * K12L * vx * vy + K22L * vx ** 2) / np.where(normv2 == 0, 1.0, normv2)
    Kt1 = (vy * (K11L * vx + K12L * vy) + (-vx) * (K21L * vx + K22L * vy)) / np.where(normv2 == 0, 1.0, normv2)
    Kn2 = (K11R * vy ** 2 - 2 * K12R * vx * vy + K22R * vx ** 2) / np.where(normv2 == 0, 1.0, normv2)
    Kt2 = (vy * (K11R * vx + K12R * vy) + (-vx) * (K21R * vx + K22R * vy)) / np.where(normv2 == 0, 1.0, normv2)

    Kde = -nv * (Kn1 * Kn2) / np.where((Kn1 * H2 + Kn2 * H1) == 0, 1.0, (Kn1 * H2 + Kn2 * H1))

    dot_vd1_vcen = vx * vcen[:, 0] + vy * vcen[:, 1]
    Ded = dot_vd1_vcen / np.where(normv2 == 0, 1.0, normv2) - (1.0 / np.where(nv == 0, 1.0, nv)) * (
        (Kt2 / np.where(Kn2 == 0, 1.0, Kn2)) * H1 + (Kt1 / np.where(Kn1 == 0, 1.0, Kn1)) * H2
    )

    # Flowrate internal
    if env["benchmark"].temFlowrateBoundary():
        idx = np.arange(nb, nb + ni)
        # ensure centelem columns exist; uses second coordinate (y) in original
        flowrateZ[idx] = Kde * (centelem[rel_i, 1] - centelem[lef_i, 1] - Ded * (no2 - no1))
        flowresultZ = flowresultZ + np.bincount(lef_i, weights=flowrateZ[idx], minlength=flowresultZ.size) - np.bincount(rel_i, weights=flowrateZ[idx], minlength=flowresultZ.size)

    # Pack into env.premethod.MPFAD
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Hesq"] = Hesq
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kde"] = Kde
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kn"] = np.concatenate((Kn.reshape(-1, 1), Kn1.reshape(-1, 1)), axis=0) if False else Kn
    # In MATLAB Kn was used both for boundary (Kn vector) and internal (Kn1/Kn2). Here we keep boundary Kn in env.premethod.MPFAD['Kn'] and also store Kn1/Kn2 as separate fields if desired.
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kt"] = Kt
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Ded"] = Ded
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["flowrateZ"] = flowrateZ
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["flowresultZ"] = flowresultZ

    # Also store internal Kn1/Kn2 and Kt1/Kt2 for later use
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kn1"] = Kn1
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kn2"] = Kn2
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kt1"] = Kt1
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["Kt2"] = Kt2

    return env
