import unittest

import numpy as np

from fs.assembly.mpfad import build


class _Benchmark:
    def __init__(self):
        self.temporal_called = False

    def calcularNeumannBoundary(
        self,
        is_neumann,
        bedge,
        bcflag,
        nflagface,
        flowrate_z,
        edge_length,
        normals,
        env,
    ):
        value_by_flag = {int(flag): value for flag, value in bcflag[:, :2]}
        flags = bedge[is_neumann, 4].astype(int)
        prescribed = np.array([value_by_flag[flag] for flag in flags])
        gravity = flowrate_z[: bedge.shape[0]][bedge[:, 4] > 200]
        return edge_length[is_neumann] * prescribed + gravity

    def adicionarTermoTemporal(self, matrix, rhs, parms, flowresult_z, env):
        self.temporal_called = True
        return matrix, rhs


class MpfadAssemblyTests(unittest.TestCase):
    def test_assembles_boundary_internal_and_nodal_rhs_terms(self):
        env = self._env()

        matrix, rhs, elembedge = build(env, {})

        neumann = 3.0 * np.sqrt(2.0) + 0.5
        np.testing.assert_allclose(
            matrix.toarray(), [[0.0, 2.0], [2.0, -2.0]]
        )
        np.testing.assert_allclose(rhs, [26.0, neumann - 2.0])
        self.assertTrue(env["benchmark"].temporal_called)
        self.assertEqual(elembedge.shape, (0, 2))

    def test_scatters_lpew_weights_for_non_dirichlet_vertices(self):
        env = self._env()
        env["config"]["nflag"][:2] = np.array([[203, 0.0], [203, 0.0]])

        matrix, _, _ = build(env, {})

        np.testing.assert_allclose(
            matrix.toarray(), [[0.75, 1.25], [1.25, -1.25]]
        )

    def test_applies_modflow_constraints(self):
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
                "esurn1": np.array([1, 1, 2, 2]),
                "esurn2": np.array([0, 1, 3, 4]),
            },
            "config": {
                "numcase": 1,
                "modflowcase": "n",
                "nflag": np.array([[100, 10.0], [100, 12.0], [100, 10.0]]),
                "nflagface": np.array([[100, 7.0], [201, 0.0]]),
                "bcflag": np.array([[201, 3.0]]),
            },
            "premethod": {
                "MPFAD": {
                    "Kde": np.array([2.0]),
                    "Ded": np.array([0.5]),
                    "Kn": np.array([2.0, 1.0]),
                    "Kt": np.array([1.0, 0.0]),
                    "Hesq": np.ones(2),
                    "flowrateZ": np.array([0.0, 0.5, 0.0]),
                    "flowresultZ": np.zeros(2),
                    "s": np.zeros(3),
                    "weight": np.array([1.0, 0.25, 0.75, 1.0]),
                }
            },
            "benchmark": _Benchmark(),
        }


if __name__ == "__main__":
    unittest.main()
