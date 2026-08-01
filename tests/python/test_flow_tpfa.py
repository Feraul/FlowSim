import unittest

import numpy as np

from fs.flow.tpfa import tpfa


class TpfaFlowTests(unittest.TestCase):
    def test_computes_boundary_internal_and_element_flows(self):
        env = self._env()

        flowrate, flowresult, flowratedif, faceaux = tpfa(
            np.array([5.0, 8.0]), env
        )

        neumann = -3.0 * np.sqrt(2.0) + 0.5
        np.testing.assert_allclose(flowrate, [-10.0, neumann, 6.0])
        np.testing.assert_allclose(flowresult, [-4.0, neumann - 6.0])
        np.testing.assert_allclose(flowratedif, np.zeros(3))
        self.assertEqual(faceaux, 0)

    def test_applies_face_viscosity(self):
        env = self._env()
        env["config"]["numcase"] = 50
        env["config"]["visc"] = np.array(
            [[1.0, 1.0], [0.5, 0.5], [2.0, 1.0]]
        )

        flowrate, _, _, _ = tpfa(np.array([5.0, 8.0]), env)

        np.testing.assert_allclose(flowrate[[0, 2]], [-20.0, 18.0])

    def test_computes_transport_dispersive_flow(self):
        env = self._env()
        env["config"]["numcase"] = 250
        env["config"]["visc"] = np.ones((3, 1))
        env["conpre"] = {"Con": np.array([0.2, 0.8]), "Kdec": np.array([4.0])}

        _, _, flowratedif, _ = tpfa(np.array([5.0, 8.0]), env)

        np.testing.assert_allclose(flowratedif, [0.0, 0.0, 2.4])

    def test_accepts_zero_based_edge_indices(self):
        env = self._env()
        env["geometry"]["bedge"][:, :3] -= 1
        env["geometry"]["inedge"][:, :4] -= 1

        flowrate, _, _, _ = tpfa(np.array([5.0, 8.0]), env)

        self.assertTrue(np.all(np.isfinite(flowrate)))

    @staticmethod
    def _env():
        return {
            "geometry": {
                "coord": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                "bedge": np.array(
                    [[1, 2, 1, 0, 100], [2, 3, 2, 0, 201]], dtype=float
                ),
                "inedge": np.array([[1, 2, 1, 2]], dtype=float),
                "centelem": np.array([[0.5, 0.5], [0.25, 0.5]]),
            },
            "config": {
                "numcase": 1,
                "nflag": np.array([[100, 10.0], [100, 10.0], [100, 10.0]]),
                "bcflag": np.array([[201, 3.0]]),
            },
            "premethod": {
                "TPFA": {
                    "Kde": np.array([2.0]),
                    "Kn": np.array([2.0, 1.0]),
                    "Hesq": np.ones(2),
                    "flowrateZ": np.array([0.0, 0.5, 0.0]),
                }
            },
        }


if __name__ == "__main__":
    unittest.main()
