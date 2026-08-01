import unittest

import numpy as np

from fs.lpew.v2.netas import netas


class NetasTests(unittest.TestCase):
    def test_computes_eta_coefficients_for_planar_corner(self):
        FS = {"csr": {"tCurrent": np.array([0]), "tNext": np.array([1])}}
        T_all = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        P_all = np.array([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        O_all = np.array([[1.0, 1.0, 0.0]])
        Q_corner = np.zeros((1, 3))

        result = netas(FS, T_all, P_all, O_all, Q_corner)

        np.testing.assert_allclose(result, np.ones((1, 2)))

    def test_batches_multiple_corners_using_csr_shifts(self):
        FS = {
            "csr": {
                "tCurrent": np.array([0, 2]),
                "tNext": np.array([1, 3]),
            }
        }
        T_all = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 1.0, 0.0],
                [1.0, 2.0, 0.0],
            ]
        )
        P_all = np.array(
            [
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [3.0, 1.0, 0.0],
                [1.0, 3.0, 0.0],
            ]
        )
        O_all = np.array([[1.0, 1.0, 0.0], [2.0, 2.0, 0.0]])
        Q_corner = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]])

        result = netas(FS, T_all, P_all, O_all, Q_corner)

        self.assertEqual(result.shape, (2, 2))
        self.assertTrue(np.all(np.isfinite(result)))


if __name__ == "__main__":
    unittest.main()
