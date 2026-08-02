"""L-scheme iteration translated from ``+fs/+iter/lscheme.m``."""

from typing import Callable, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as sparse_linalg

from .common import ConvergenceError, validate_options, vector


def lscheme(
    matrix,
    rhs: np.ndarray,
    x0: np.ndarray,
    *,
    assemble: Optional[Callable[[np.ndarray], Tuple[object, np.ndarray]]] = None,
    tolerance: float = 1e-7,
    jump_tolerance: float = 1e-3,
    max_iterations: int = 100,
    stabilization: float = 5.0,
    relaxation: float = 1.0,
    max_jump: float = 0.2,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """Solve a nonlinear system with an optional state-dependent reassembly."""
    validate_options(tolerance, max_iterations)
    if stabilization <= 0:
        raise ValueError("stabilization must be positive")
    if not 0 < relaxation <= 1:
        raise ValueError("relaxation must be in (0, 1]")
    if max_jump <= 0 or jump_tolerance < 0:
        raise ValueError("jump limits must be positive")
    if assemble is not None and not callable(assemble):
        raise TypeError("assemble must be callable")

    x = vector(x0, name="x0")
    matrix = _matrix(matrix, x.size)
    rhs = vector(rhs, x.size, "rhs")
    initial_residual = np.linalg.norm(matrix @ x - rhs)
    if initial_residual == 0:
        return x, 0, np.empty((0, 2))

    history = []
    for iteration in range(1, max_iterations + 1):
        shifted = matrix + stabilization * sparse.eye(x.size, format="csr")
        shifted_rhs = rhs + stabilization * x
        candidate, info = sparse_linalg.cg(
            shifted, shifted_rhs, rtol=1e-8, atol=0.0, maxiter=500
        )
        if info != 0:
            candidate = sparse_linalg.spsolve(shifted, shifted_rhs)
        candidate = vector(candidate, x.size, "linear solve result")

        step = relaxation * (candidate - x)
        largest_jump = np.max(np.abs(step)) if step.size else 0.0
        if largest_jump > max_jump:
            step *= max_jump / largest_jump
        x_new = x + step
        largest_jump = np.max(np.abs(step)) if step.size else 0.0

        if assemble is not None:
            matrix, rhs = assemble(x_new.copy())
            matrix = _matrix(matrix, x.size)
            rhs = vector(rhs, x.size, "assembled rhs")

        relative_residual = np.linalg.norm(matrix @ x_new - rhs) / initial_residual
        history.append((iteration, relative_residual))
        x = x_new
        if relative_residual <= tolerance and largest_jump <= jump_tolerance:
            return x, iteration, np.asarray(history)

    raise ConvergenceError(
        f"L-scheme did not converge after {max_iterations} iterations"
    )


def _matrix(value, size: int) -> sparse.csr_matrix:
    result = sparse.csr_matrix(value, dtype=float)
    if result.shape != (size, size):
        raise ValueError(f"matrix has shape {result.shape}, expected ({size}, {size})")
    if not np.all(np.isfinite(result.data)):
        raise ValueError("matrix contains non-finite values")
    return result
