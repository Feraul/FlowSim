import unittest

import numpy as np

from fs.flow.mpfad import mpfad


class MpfadFlowTests(unittest.TestCase):
    def test_computes_tangential_boundary_and_interpolated_internal_flow(self):
        env = self._env()

        flowrate, flowresult, flowratedif, faceaux = mpfad(
            np.array([5.0, 8.0]), np.array([7.0, 9.0, 8.0]), env
        )

        neumann = -3.0 * np.sqrt(2.0) + 0.5
        np.testing.assert_allclose(flowrate, [-14.0, neumann, 4.0])
        np.testing.assert_allclose(flowresult, [-10.0, neumann - 4.0])
        np.testing.assert_allclose(flowratedif, np.zeros(3))
        self.assertEqual(faceaux, 0)

    def test_computes_transport_dispersive_flow_with_nodal_interpolation(self):
        env = self._env()
        env["config"]["numcase"] = 250
        env["config"]["visc"] = np.ones((3, 1))
        env["geometry"]["elem"] = np.zeros((2, 5))
        env["geometry"]["esurn1"] = np.array([1, 1, 2, 2])
        env["geometry"]["esurn2"] = np.array([0, 1, 3, 4])
        env["premethod"]["MPFAD"]["weight"] = np.array([1.0, 0.5, 0.5, 1.0])
        env["premethod"]["MPFAD"]["s"] = np.zeros(3)
        env["conpre"] = {
            "Con": np.array([0.2, 0.8]),
            "wightc": np.array([1.0, 0.5, 0.5, 1.0]),
            "sc": np.zeros(3),
            "nflagnoc": np.array([[201, 0.0], [201, 0.0], [201, 0.0]]),
            "Kdec": np.array([4.0]),
            "Dedc": np.array([0.25]),
        }

        _, _, flowratedif, _ = mpfad(
            np.array([5.0, 8.0]), np.array([5.0, 6.5, 8.0]), env
        )

        np.testing.assert_allclose(flowratedif, [0.0, 0.0, 2.1])

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
                "nflag": np.array([[100, 10.0], [100, 12.0], [100, 10.0]]),
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
                }
            },
        }


if __name__ == "__main__":
    unittest.main()
