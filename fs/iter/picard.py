"""Relaxed Picard iteration translated from ``+fs/+iter/picard.m``."""

from typing import Callable, Tuple

import numpy as np

from .common import ConvergenceError, validate_options, vector


def picard(
    x0: np.ndarray,
    fixed_point: Callable[[np.ndarray], np.ndarray],
    *,
    tolerance: float = 1e-7,
    max_iterations: int = 100,
    relaxation: float = 0.5,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """Solve ``x = fixed_point(x)`` with relaxed fixed-point iteration."""
    validate_options(tolerance, max_iterations)
    if not callable(fixed_point):
        raise TypeError("fixed_point must be callable")
    if not 0 < relaxation <= 1:
        raise ValueError("relaxation must be in (0, 1]")

    x = vector(x0, name="x0")
    history = []
    for iteration in range(1, max_iterations + 1):
        candidate = vector(fixed_point(x.copy()), x.size, "fixed_point result")
        x_new = x + relaxation * (candidate - x)
        residual = np.linalg.norm(x_new - x)
        history.append((iteration, residual))
        x = x_new
        if residual <= tolerance * (1.0 + np.linalg.norm(x)):
            return x, iteration, np.asarray(history)

    raise ConvergenceError(
        f"Picard iteration did not converge after {max_iterations} iterations"
    )
