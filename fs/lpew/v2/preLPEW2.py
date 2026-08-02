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
- fs.mesh.build normaliza os índices CSR legados para 0-based

Retorna: (env, weight, s)
- weight: 1D numpy array (row) com tamanho nCorners (legacy: 1 x nCorners row-vector)
- s: 1D numpy array shape (nNodes,)
"""
from typing import Tuple, Dict, Any
import numpy as np


def preLPEW2(env: Dict[str, Any], parms: Dict[str, Any]) -> Tuple[Dict[str, Any], np.ndarray, np.ndarray]:
    # Ensure kmap as MATLAB code does
    config = env.setdefault("config", {})
    auxperm = parms.get("auxperm")
    if config.get("numcase", 0) > 400 and not _is_empty(auxperm):
        config["kmap"] = auxperm
    elif _is_empty(config.get("kmap")):
        config["kmap"] = _default_kmap(config.get("perm"))

    # Build FS scaffolding
    from fs.mesh import build as mesh_build
    from fs.csr import buildCorners, buildCornerShifts

    FS = mesh_build(env)
    FS = buildCorners(FS)
    FS = buildCornerShifts(FS)
    # Keep compatibility with MATLAB: FS.cfg.phasekey
    FS.setdefault("cfg", {})["phasekey"] = env.get("config", {}).get("phasekey")

    # Batched kernels
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
    nNodes = int(FS["mesh"]["nNodes"])

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
    N = env["premethod"]["MPFAD"]["N"]
    benchmark = env["benchmark"]
    if isinstance(benchmark, dict):
        neumann_callback = benchmark.get("calcularTermoNeumannVet")
    else:
        neumann_callback = getattr(benchmark, "calcularTermoNeumannVet", None)
    if not callable(neumann_callback):
        raise TypeError("env.benchmark.calcularTermoNeumannVet must be callable")
    s = np.asarray(neumann_callback(r, sum_lambda, N, env)).reshape(-1)

    # Persist in env.premethod
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["weight"] = weight
    env.setdefault("premethod", {}).setdefault("MPFAD", {})["s"] = s

    return env, weight, s


def _is_empty(value: Any) -> bool:
    return value is None or np.asarray(value).size == 0


def _default_kmap(perm: Any) -> np.ndarray:
    if _is_empty(perm):
        raise ValueError("env.config.perm is required when kmap is empty")
    tensor = np.asarray(perm)
    if tensor.ndim == 1:
        tensor = tensor.reshape(1, -1)
    if tensor.ndim != 2:
        raise ValueError("env.config.perm must be a one- or two-dimensional array")
    if tensor.shape[1] >= 5:
        return tensor
    if tensor.shape[1] != 4:
        raise ValueError("env.config.perm must contain K11/K12/K21/K22")
    material_ids = np.arange(1, tensor.shape[0] + 1).reshape(-1, 1)
    return np.column_stack((material_ids, tensor))
