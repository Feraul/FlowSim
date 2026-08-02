import unittest

import numpy as np

from fs.iter import ConvergenceError, anderson, lscheme, picard


class IteratorTests(unittest.TestCase):
    def test_picard_converges_with_relaxation(self):
        solution, iterations, history = picard(
            np.array([0.0]),
            lambda x: np.cos(x),
            tolerance=1e-10,
            max_iterations=200,
            relaxation=0.8,
        )

        np.testing.assert_allclose(solution, [0.7390851332], atol=1e-9)
        self.assertEqual(history.shape, (iterations, 2))

    def test_anderson_converges_for_linear_fixed_point(self):
        solution, iterations, history = anderson(
            np.array([0.0]),
            lambda x: 0.5 * x + 1.0,
            tolerance=1e-12,
            max_iterations=20,
            memory=3,
        )

        np.testing.assert_allclose(solution, [2.0], atol=1e-11)
        self.assertEqual(history[-1, 0], iterations)

    def test_lscheme_solves_stabilized_linear_problem(self):
        solution, iterations, history = lscheme(
            np.array([[2.0]]),
            np.array([4.0]),
            np.array([0.0]),
            tolerance=1e-10,
            jump_tolerance=1e-9,
            max_iterations=100,
            stabilization=1.0,
            max_jump=10.0,
        )

        np.testing.assert_allclose(solution, [2.0], atol=1e-8)
        self.assertEqual(history.shape, (iterations, 2))

    def test_lscheme_reassembles_state_dependent_system(self):
        calls = []

        def assemble(x):
            calls.append(x.copy())
            return np.array([[1.0 + 0.1 * x[0]]]), np.array([1.0])

        solution, _, _ = lscheme(
            np.array([[1.0]]),
            np.array([1.0]),
            np.array([0.0]),
            assemble=assemble,
            tolerance=1e-8,
            jump_tolerance=1e-7,
            max_iterations=100,
            stabilization=1.0,
            max_jump=10.0,
        )

        self.assertTrue(calls)
        self.assertAlmostEqual((1.0 + 0.1 * solution[0]) * solution[0], 1.0, places=7)

    def test_reports_nonconvergence(self):
        with self.assertRaises(ConvergenceError):
            picard(
                np.array([0.0]),
                lambda x: x + 1.0,
                max_iterations=2,
            )


if __name__ == "__main__":
    unittest.main()
