import unittest

import numpy as np

from fs.lpew.pinterp import pinterp


class PinterpTests(unittest.TestCase):
    def test_interpolates_pressure_and_applies_boundary_conditions(self):
        env = self._env(numcase=100)

        pressure, concentration = pinterp(np.array([10.0, 20.0]), env)

        np.testing.assert_allclose(pressure, [7.0, 19.5, 20.0])
        self.assertEqual(concentration, 0.0)

    def test_interpolates_transport_concentration_from_conpre(self):
        env = self._env(numcase=250)
        env["conpre"] = {
            "Con": np.array([0.2, 0.8]),
            "wightc": np.array([1.0, 0.5, 0.5, 1.0]),
            "sc": np.array([0.0, 0.1, 0.0]),
            "nflagnoc": np.array([[100, 0.9], [202, 0.0], [201, 0.0]]),
        }

        _, concentration = pinterp(np.array([10.0, 20.0]), env)

        np.testing.assert_allclose(concentration, [0.9, 0.6, 0.8])

    def test_accepts_zero_based_element_indices(self):
        env = self._env(numcase=100)
        env["geometry"]["esurn1"] = np.array([0, 0, 1, 1])

        pressure, _ = pinterp(np.array([10.0, 20.0]), env)

        np.testing.assert_allclose(pressure, [7.0, 19.5, 20.0])

    def test_rejects_weight_shape_mismatch(self):
        env = self._env(numcase=100)
        env["premethod"]["MPFAD"]["weight"] = np.ones(3)

        with self.assertRaisesRegex(ValueError, "expected 4"):
            pinterp(np.array([10.0, 20.0]), env)

    @staticmethod
    def _env(numcase):
        return {
            "geometry": {
                "coord": np.zeros((3, 2)),
                "elem": np.zeros((2, 5)),
                "esurn1": np.array([1, 1, 2, 2]),
                "esurn2": np.array([0, 1, 3, 4]),
            },
            "config": {
                "numcase": numcase,
                "nflag": np.array([[100, 7.0], [202, 0.0], [201, 0.0]]),
            },
            "premethod": {
                "MPFAD": {
                    "weight": np.array([1.0, 0.25, 0.75, 1.0]),
                    "s": np.array([0.0, 2.0, 0.0]),
                }
            },
        }


if __name__ == "__main__":
    unittest.main()
