import unittest

import numpy as np
from scipy import sparse

from fs.assembly.tpfa import build


class _TemporalBenchmark:
    def __init__(self):
        self.called = False

    def adicionarTermoTemporal(self, matrix, rhs, parms, flowresult_z, env):
        self.called = True
        return matrix + sparse.eye(matrix.shape[0]), rhs + flowresult_z


class TpfaAssemblyTests(unittest.TestCase):
    def test_assembles_sparse_matrix_and_boundary_rhs(self):
        matrix, rhs, elembedge = build(self._env(), {})

        expected_matrix = np.array([[0.0, 2.0], [2.0, -2.0]])
        expected_rhs = np.array([22.0, 3.0 * np.sqrt(2.0) + 0.5])
        np.testing.assert_allclose(matrix.toarray(), expected_matrix)
        np.testing.assert_allclose(rhs, expected_rhs)
        self.assertEqual(elembedge.shape, (0, 2))

    def test_invokes_temporal_benchmark_hook(self):
        env = self._env()
        benchmark = _TemporalBenchmark()
        env["config"]["numcase"] = 439
        env["benchmark"] = benchmark
        env["premethod"]["TPFA"]["flowresultZ"] = np.array([0.25, 0.75])

        matrix, rhs, _ = build(env, {"dt": 1.0})

        self.assertTrue(benchmark.called)
        np.testing.assert_allclose(
            matrix.toarray(), np.array([[1.0, 2.0], [2.0, -1.0]])
        )
        np.testing.assert_allclose(
            rhs, [22.25, 3.0 * np.sqrt(2.0) + 1.25]
        )

    def test_applies_modflow_element_constraints(self):
        env = self._env()
        env["config"]["modflowcase"] = "y"

        matrix, rhs, elembedge = build(env, {})

        np.testing.assert_allclose(matrix.toarray()[0], [1.0, 0.0])
        self.assertEqual(rhs[0], 7.0)
        np.testing.assert_allclose(elembedge, [[0.0, 7.0]])

    @staticmethod
    def _env():
        return {
            "geometry": {
                "coord": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                "elem": np.zeros((2, 5)),
                "bedge": np.array(
                    [[1, 2, 1, 0, 100], [2, 3, 2, 0, 201]], dtype=float
                ),
                "inedge": np.array([[1, 2, 1, 2]], dtype=float),
                "centelem": np.array([[0.5, 0.5], [0.25, 0.5]]),
                "normals": np.zeros((2, 2)),
            },
            "config": {
                "numcase": 1,
                "modflowcase": "n",
                "nflag": np.array([[100, 10.0], [100, 12.0], [100, 10.0]]),
                "nflagface": np.array([[100, 7.0], [201, 0.0]]),
                "bcflag": np.array([[201, 3.0]]),
            },
            "premethod": {
                "TPFA": {
                    "Kde": np.array([2.0]),
                    "Kn": np.array([2.0, 1.0]),
                    "Hesq": np.ones(2),
                    "flowrateZ": np.array([0.0, 0.5, 0.0]),
                    "flowresultZ": np.zeros(2),
                }
            },
        }


if __name__ == "__main__":
    unittest.main()
