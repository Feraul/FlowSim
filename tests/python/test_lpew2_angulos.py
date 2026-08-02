import unittest

import numpy as np

from fs.lpew.v2.angulos import angulos


class AngulosTests(unittest.TestCase):
    def test_computes_corner_angles_in_matlab_output_order(self):
        FS = {"csr": {"tCurrent": np.array([0]), "tNext": np.array([1])}}
        T_all = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        O_all = np.array([[1.0, 1.0, 0.0]])
        Q_corner = np.zeros((1, 3))

        ve2, ve1, theta2, theta1 = angulos(FS, T_all, O_all, Q_corner)

        expected = np.array([np.pi / 4])
        np.testing.assert_allclose(ve2, expected)
        np.testing.assert_allclose(ve1, expected)
        np.testing.assert_allclose(theta2, expected)
        np.testing.assert_allclose(theta1, expected)

    def test_batches_multiple_corners_using_csr_shifts(self):
        FS = {
            "csr": {
                "tCurrent": np.array([0, 2]),
                "tNext": np.array([1, 3]),
            }
        }
        T_all = np.array(
            [
                [2.0, 0.0],
                [0.0, 2.0],
                [2.0, 1.0],
                [1.0, 3.0],
            ]
        )
        O_all = np.array([[1.0, 1.0], [2.0, 2.0]])
        Q_corner = np.array([[0.0, 0.0], [1.0, 1.0]])

        outputs = angulos(FS, T_all, O_all, Q_corner)

        self.assertTrue(all(output.shape == (2,) for output in outputs))
        self.assertTrue(all(np.all(np.isfinite(output)) for output in outputs))


if __name__ == "__main__":
    unittest.main()
