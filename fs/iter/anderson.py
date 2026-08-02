"""Anderson acceleration translated from ``+fs/+iter/anderson.m``."""

from typing import Callable, Tuple

import numpy as np

from .common import ConvergenceError, validate_options, vector


def anderson(
    x0: np.ndarray,
    fixed_point: Callable[[np.ndarray], np.ndarray],
    *,
    tolerance: float = 1e-7,
    max_iterations: int = 100,
    memory: int = 10,
    start: int = 1,
    beta: float = 1.0,
) -> Tuple[np.ndarray, int, np.ndarray]:
    """Solve ``x = fixed_point(x)`` using damped Anderson acceleration."""
    validate_options(tolerance, max_iterations)
    if not callable(fixed_point):
        raise TypeError("fixed_point must be callable")
    if memory < 0:
        raise ValueError("memory must be non-negative")
    if start < 0:
        raise ValueError("start must be non-negative")
    if not 0 < beta <= 1:
        raise ValueError("beta must be in (0, 1]")

    x = vector(x0, name="x0")
    history = []
    delta_f = []
    delta_g = []
    previous_f = None
    previous_g = None

    for iteration in range(max_iterations + 1):
        g_value = vector(fixed_point(x.copy()), x.size, "fixed_point result")
        residual_vector = g_value - x
        residual = np.linalg.norm(residual_vector)
        history.append((iteration, residual))
        if residual <= tolerance * (1.0 + np.linalg.norm(g_value)):
            return g_value, iteration, np.asarray(history)

        if previous_f is not None and memory:
            delta_f.append(residual_vector - previous_f)
            delta_g.append(g_value - previous_g)
            if len(delta_f) > memory:
                delta_f.pop(0)
                delta_g.pop(0)

        if iteration < start or not delta_f:
            x_new = g_value
        else:
            residual_differences = np.column_stack(delta_f)
            value_differences = np.column_stack(delta_g)
            coefficients = np.linalg.lstsq(
                residual_differences, residual_vector, rcond=None
            )[0]
            accelerated = g_value - value_differences @ coefficients
            x_new = beta * accelerated + (1.0 - beta) * g_value

        previous_f = residual_vector
        previous_g = g_value
        x = vector(x_new, x.size, "Anderson update")

    raise ConvergenceError(
        f"Anderson iteration did not converge after {max_iterations} iterations"
    )
