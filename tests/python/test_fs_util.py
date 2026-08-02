import unittest

import numpy as np

from fs.util import assertFS


class AssertFSTests(unittest.TestCase):
    def test_accepts_valid_populated_structure(self):
        assertFS(self._valid_fs())

    def test_rejects_missing_top_level_field(self):
        FS = self._valid_fs()
        del FS["perm"]

        with self.assertRaisesRegex(ValueError, "perm"):
            assertFS(FS)

    def test_rejects_invalid_mesh_pointer_length(self):
        FS = self._valid_fs()
        FS["mesh"]["esurn2"] = np.array([0, 1])

        with self.assertRaisesRegex(ValueError, "expected 3"):
            assertFS(FS)

    def test_rejects_invalid_corner_pointer(self):
        FS = self._valid_fs()
        FS["csr"]["nodePtr"][-1] = 1

        with self.assertRaisesRegex(ValueError, "nCorners"):
            assertFS(FS)

    def test_rejects_geometry_or_permeability_row_drift(self):
        FS = self._valid_fs()
        FS["geom"]["centElem"] = np.zeros((1, 2))
        with self.assertRaisesRegex(ValueError, "centElem"):
            assertFS(FS)

        FS = self._valid_fs()
        FS["perm"]["tensor"] = np.zeros((1, 4))
        with self.assertRaisesRegex(ValueError, "tensor"):
            assertFS(FS)

    @staticmethod
    def _valid_fs():
        return {
            "mesh": {
                "nNodes": 2,
                "nElems": 2,
                "esurn2": np.array([0, 1, 3]),
            },
            "geom": {"centElem": np.zeros((2, 2))},
            "csr": {"nCorners": 3, "nodePtr": np.array([0, 1, 3])},
            "perm": {"tensor": np.zeros((2, 4))},
            "bc": {},
            "cfg": {},
            "state": {},
            "workspace": {},
        }


if __name__ == "__main__":
    unittest.main()
