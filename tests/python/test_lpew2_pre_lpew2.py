import unittest

import numpy as np

from fs.lpew.v2.preLPEW2 import preLPEW2


class _Benchmark:
    def __init__(self):
        self.call = None

    def calcularTermoNeumannVet(self, r, sum_lambda, N, env):
        self.call = (r.copy(), sum_lambda.copy(), N, env)
        return np.sum(r, axis=1) + sum_lambda + N


class PreLPEW2Tests(unittest.TestCase):
    def test_runs_complete_pipeline_and_persists_results(self):
        benchmark = _Benchmark()
        auxperm = np.array([[1, 2.0, 0.0, 0.0, 2.0]])
        env = self._triangle_env(benchmark)

        result_env, weight, s = preLPEW2(env, {"auxperm": auxperm})

        self.assertIs(result_env, env)
        np.testing.assert_allclose(weight, np.ones(3))
        np.testing.assert_allclose(env["config"]["kmap"], auxperm)
        np.testing.assert_allclose(env["premethod"]["MPFAD"]["weight"], weight)
        np.testing.assert_allclose(env["premethod"]["MPFAD"]["s"], s)
        self.assertIsNotNone(benchmark.call)
        self.assertEqual(benchmark.call[0].shape, (3, 2))
        self.assertEqual(benchmark.call[1].shape, (3,))
        self.assertEqual(s.shape, (3,))
        self.assertTrue(np.all(np.isfinite(s)))

    def test_requires_neumann_callback(self):
        env = self._triangle_env({})

        with self.assertRaisesRegex(TypeError, "calcularTermoNeumannVet"):
            preLPEW2(env, {})

    @staticmethod
    def _triangle_env(benchmark):
        return {
            "geometry": {
                "coord": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                "elem": np.array([[1, 2, 3, 0, 1]]),
                "bedge": np.zeros((3, 5)),
                "inedge": np.zeros((0, 6)),
                "nsurn1": np.array([2, 3, 3, 1, 1, 2]),
                "nsurn2": np.array([0, 2, 4, 6]),
                "esurn1": np.array([1, 1, 1]),
                "esurn2": np.array([0, 1, 2, 3]),
                "centelem": np.array([[1.0 / 3, 1.0 / 3]]),
                "elemarea": np.array([0.5]),
            },
            "config": {
                "numcase": 401,
                "phasekey": 1,
                "perm": np.array([[2.0, 0.0, 0.0, 2.0]]),
            },
            "premethod": {"MPFAD": {"N": np.ones(3)}},
            "benchmark": benchmark,
        }


if __name__ == "__main__":
    unittest.main()
