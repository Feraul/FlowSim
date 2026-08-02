import unittest

import numpy as np

from fs.lpew.v2.lambdaWeights import lambdaWeights


class LambdaWeightsTests(unittest.TestCase):
    def test_computes_vectorized_interior_weights(self):
        FS = self._fs(node_is_interior=True, corners=3, neighbors=3)
        inputs = self._uniform_inputs(corners=3, neighbors=3)

        lambda_values, r = lambdaWeights(FS, *inputs)

        np.testing.assert_allclose(lambda_values, np.full(3, 2.0))
        np.testing.assert_allclose(r, np.zeros((1, 2)))

    def test_computes_boundary_weights_and_neumann_distances(self):
        FS = self._fs(node_is_interior=False, corners=2, neighbors=3)
        inputs = self._uniform_inputs(corners=2, neighbors=3)

        lambda_values, r = lambdaWeights(FS, *inputs)

        np.testing.assert_allclose(lambda_values, np.full(2, 2.0))
        np.testing.assert_allclose(r, np.array([[2.0, 6.0]]))

    def test_rejects_invalid_boundary_layout(self):
        FS = self._fs(node_is_interior=False, corners=2, neighbors=2)
        inputs = self._uniform_inputs(corners=2, neighbors=2)

        with self.assertRaisesRegex(ValueError, "one more neighbor"):
            lambdaWeights(FS, *inputs)

    @staticmethod
    def _fs(node_is_interior, corners, neighbors):
        return {
            "mesh": {
                "nNodes": 1,
                "esurn2": np.array([0, corners]),
                "nsurn2": np.array([0, neighbors]),
                "coord": np.zeros((1, 2)),
            },
            "csr": {
                "nCorners": corners,
                "cornerNode": np.zeros(corners, dtype=int),
                "cornerPrev": np.roll(np.arange(corners), 1),
                "cornerNext": np.roll(np.arange(corners), -1),
                "nodeIsInterior": np.array([node_is_interior]),
            },
        }

    @staticmethod
    def _uniform_inputs(corners, neighbors):
        kt1 = np.zeros((corners, 2))
        kt2 = np.zeros(corners)
        kn1 = np.ones((corners, 2))
        kn2 = np.ones(corners)
        angle = np.full(corners, np.pi / 4)
        netas = np.ones((corners, 2))
        t_all = np.column_stack(
            (np.arange(1, neighbors + 1, dtype=float), np.zeros(neighbors))
        )
        q_corner = np.zeros((corners, 2))
        return (
            kt1,
            kt2,
            kn1,
            kn2,
            angle,
            angle,
            angle,
            angle,
            netas,
            t_all,
            q_corner,
        )


if __name__ == "__main__":
    unittest.main()
