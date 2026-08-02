import unittest

import numpy as np

from fs.csr import buildCornerShifts, buildCorners


class CsrBuilderTests(unittest.TestCase):
    def test_builds_flat_corner_layout(self):
        FS = {
            "mesh": {
                "nNodes": 3,
                "esurn1": np.array([0, 0, 1, 1]),
                "esurn2": np.array([0, 1, 3, 4]),
            },
            "csr": {},
        }

        result = buildCorners(FS)

        self.assertIs(result, FS)
        self.assertEqual(FS["csr"]["nCorners"], 4)
        np.testing.assert_array_equal(FS["csr"]["nodeNec"], [1, 2, 1])
        np.testing.assert_array_equal(FS["csr"]["cornerNode"], [0, 1, 1, 2])
        np.testing.assert_array_equal(FS["csr"]["cornerElem"], [0, 0, 1, 1])
        self.assertEqual(FS["csr"]["maxNec"], 2)

    def test_builds_boundary_and_interior_shifts(self):
        FS = {
            "mesh": {
                "nNodes": 2,
                "esurn1": np.arange(5),
                "esurn2": np.array([0, 2, 5]),
                "nsurn1": np.arange(6),
                "nsurn2": np.array([0, 3, 6]),
            },
            "csr": {},
        }
        buildCorners(FS)

        buildCornerShifts(FS)

        csr = FS["csr"]
        np.testing.assert_array_equal(csr["cornerLocal"], [0, 1, 0, 1, 2])
        np.testing.assert_array_equal(csr["nodeIsInterior"], [False, True])
        np.testing.assert_array_equal(csr["tCurrent"], [0, 1, 3, 4, 5])
        np.testing.assert_array_equal(csr["tNext"], [1, 2, 4, 5, 3])
        np.testing.assert_array_equal(csr["cornerPrev"], [0, 0, 4, 2, 3])
        np.testing.assert_array_equal(csr["cornerNext"], [1, 1, 3, 4, 2])
        np.testing.assert_array_equal(
            csr["boundaryFirstMask"], [True, False, False, False, False]
        )
        np.testing.assert_array_equal(
            csr["boundaryLastMask"], [False, True, False, False, False]
        )

    def test_rejects_invalid_corner_pointer(self):
        FS = {
            "mesh": {
                "nNodes": 1,
                "esurn1": np.array([0]),
                "esurn2": np.array([1, 1]),
            },
            "csr": {},
        }

        with self.assertRaisesRegex(ValueError, "valid CSR"):
            buildCorners(FS)


if __name__ == "__main__":
    unittest.main()
