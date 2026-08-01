# -*- coding: utf-8 -*-
"""
fs.lpew.v2.preLPEW2 (Python) — tradução de +fs/+lpew/+v2/preLPEW2.m

Pipeline LPEW2 vetorizado:
- constrói estruturas FS + CSR + shifts
- executa kernels batched: OPT, ksInterp, angulos, netas, lambdaWeights
- normaliza lambda por nó para obter weight (1 x nCorners)
- calcula termo de Neumann s via benchmark.calcularTermoNeumannVet
- persiste env.premethod.MPFAD.weight e env.premethod.MPFAD.s

Assunções:
- env e parms são dicts com estruturas numpy-friendly
- índices já convertidos para 0-based quando necessário
- os kernels fs.mesh.build, fs.csr.buildCorners, fs.csr.buildCornerShifts,
  fs.lpew.OPT, fs.lpew.v2.ksInterp, fs.lpew.v2.angulos, fs.lpew.v2.netas,
  fs.lpew.v2.lambdaWeights existem (serão traduzidos separadamente)

Retorna: (env, weight, s)
- weight: 1D numpy array (row) com tamanho nCorners (legacy: 1 x nCorners row-vector)
- s: 1D numpy array shape (nNodes,)
"""
from typing import Tuple, Dict, Any
import numpy as np


def preLPEW2(env: Dict[str, Any], parms: Dict[str, Any]) -> Tuple[Dict[str, Any], np.ndarray, np.ndarray]:
    # Ensure kmap as MATLAB code does
    if env.get("config", {}).get("numcase", 0) > 400 and parms.get("auxperm") is not None and parms.get("auxperm") != []:
        env.setdefault("config", {})["kmap"] = parms.get("auxperm")
    elif env.get("config", {}).get("kmap") is None or env.get("config", {}).get("kmap") == []:
        env.setdefault("config", {})["kmap"] = env.get("config", {}).get("perm")

    # Build FS scaffolding
    # These functions/modules must be implemented separately in Python.
    from fs.mesh import build as mesh_build
    from fs.csr import buildCorners, buildCornerShifts

    FS = mesh_build(env)
    FS = buildCorners(FS)
    FS = buildCornerShifts(FS)
    # Keep compatibility with MATLAB: FS.cfg.phasekey
    FS.setdefault("cfg", {})["phasekey"] = env.get("config", {}).get("phasekey")

    # Batched kernels
    # Each kernel should operate on FS and return numpy arrays shaped like in MATLAB
    # The names mirror the MATLAB +fs pipeline; implement them in fs.lpew and fs.lpew.v2
    from fs.lpew import OPT
    from fs.lpew.v2 import ksInterp, angulos, netas, lambdaWeights

    P_all, T_all, O_all, Qc_all, _ = OPT(FS)
    Kt1, Kt2, Kn1, Kn2 = ksInterp(FS, T_all, Qc_all)
    ve2, ve1, theta2, theta1 = angulos(FS, T_all, O_all, Qc_all)
    netas_all = netas(FS, T_all, P_all, O_all, Qc_all)
    lambda_arr, r = lambdaWeights(FS, Kt1, Kt2, Kn1, Kn2,
                                   theta1, theta2, ve1, ve2, netas_all, T_all, Qc_all)

    # Segmented normalization: weight = lambda / sum(lambda per node)
    # FS.csr.cornerNode: array length nCorners mapping each corner -> node index (0-based)
    csr = FS.get("csr")
    if csr is None or "cornerNode" not in csr:
        raise RuntimeError("FS.csr.cornerNode not found — ensure fs.csr.buildCorners produced it")

    cornerNode = np.asarray(csr["cornerNode"]).astype(int)
    nNodes = int(FS.get("mesh", {}).get("nNodes", np.max(cornerNode) + 1))

    # sum_lambda per node (equivalent a accumarray in MATLAB)
    sum_lambda = np.bincount(cornerNode, weights=lambda_arr, minlength=nNodes)

    # denom per corner: sum_lambda mapped by cornerNode
    denom = sum_lambda[cornerNode]
    # Guard against division by zero
    denom_safe = np.where(denom == 0, 1.0, denom)

    weight_col = np.asarray(lambda_arr) / denom_safe

    # Legacy layout: row-vector (1 x nCorners). We'll return a 1D numpy array but
    # document it as row layout (transpose if needed by callers)
    weight = weight_col.reshape(-1)

    # Neumann source term (benchmark-provided) — N is env.premethod.MPFAD.N
    N = env.get("premethod", {}).get("MPFAD", {}).get("N")
    if N is None:
        # If N is missing, we cannot compute s; default to zeros
        nNodes_local = nNodes
        s = np.zeros(nNodes_local, dtype=float)
    else:
        # Benchmark callback expects (r, sum_lambda, N, env)
        # sum_lambda must be shape (nNodes, ) as in MATLAB
        s = env["benchmark"].calcularTermoNeumannVet(r, sum_lambda, N, env)
        s = np.asarray(s).reshape(-1)

    # Persist in env.premethod
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["weight"] = weight
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["s"] = s

    return env, weight, s
