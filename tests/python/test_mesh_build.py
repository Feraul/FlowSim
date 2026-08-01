import unittest

import numpy as np

from fs.mesh import build


class MeshBuildTests(unittest.TestCase):
    def test_adapts_legacy_environment_and_normalizes_flat_indices(self):
        env = {
            "geometry": {
                "coord": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                "elem": np.array([[1, 2, 3, 0, 7]]),
                "bedge": np.zeros((3, 5)),
                "inedge": np.zeros((0, 6)),
                "nsurn1": np.array([2, 3, 1, 3, 1, 2]),
                "nsurn2": np.array([0, 2, 4, 6]),
                "esurn1": np.array([1, 1, 1]),
                "esurn2": np.array([0, 1, 2, 3]),
                "centelem": np.array([[1.0 / 3, 1.0 / 3]]),
                "elemarea": np.array([0.5]),
                "normals": np.zeros((3, 2)),
            },
            "config": {
                "phasekey": 1,
                "perm": np.array([[2.0, 0.0, 0.0, 2.0]]),
                "kmap": np.array([[7, 2.0, 0.0, 0.0, 2.0]]),
                "nflag": np.zeros((3, 2)),
            },
        }

        FS = build(env)

        self.assertEqual(FS["mesh"]["nNodes"], 3)
        self.assertEqual(FS["mesh"]["nElems"], 1)
        self.assertEqual(FS["mesh"]["nBFaces"], 3)
        self.assertEqual(FS["mesh"]["nIFaces"], 0)
        np.testing.assert_array_equal(FS["mesh"]["nsurn1"], [1, 2, 0, 2, 0, 1])
        np.testing.assert_array_equal(FS["mesh"]["esurn1"], [0, 0, 0])
        np.testing.assert_array_equal(env["geometry"]["esurn1"], [1, 1, 1])
        self.assertEqual(FS["cfg"]["phasekey"], 1)
        self.assertIn("normalInt", FS["geom"])
        self.assertEqual(FS["csr"], {})

    def test_preserves_indices_that_are_already_zero_based(self):
        env = {
            "geometry": {
                "coord": np.zeros((2, 2)),
                "elem": np.zeros((1, 5)),
                "bedge": np.zeros((0, 5)),
                "inedge": np.zeros((0, 6)),
                "nsurn1": np.array([1, 0]),
                "nsurn2": np.array([0, 1, 2]),
                "esurn1": np.array([0, 0]),
                "esurn2": np.array([0, 1, 2]),
                "centelem": np.zeros((1, 2)),
                "elemarea": np.ones(1),
            }
        }

        FS = build(env)

        np.testing.assert_array_equal(FS["mesh"]["nsurn1"], [1, 0])
        np.testing.assert_array_equal(FS["mesh"]["esurn1"], [0, 0])


if __name__ == "__main__":
    unittest.main()
