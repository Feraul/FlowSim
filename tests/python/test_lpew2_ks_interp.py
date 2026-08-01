import unittest

import numpy as np

from fs.lpew.v2.ksInterp import ksInterp


class KsInterpTests(unittest.TestCase):
    def test_projects_isotropic_permeability(self):
        FS = self._fs([[7, 2.0, 0.0, 0.0, 2.0]], [7])
        T_all = np.array([[1.0, 0.0], [0.0, 1.0]])
        Q_corner = np.zeros((1, 2))

        kt1, kt2, kn1, kn2 = ksInterp(FS, T_all, Q_corner)

        np.testing.assert_allclose(kt1, np.zeros((1, 2)))
        np.testing.assert_allclose(kt2, np.zeros(1))
        np.testing.assert_allclose(kn1, np.full((1, 2), 2.0))
        np.testing.assert_allclose(kn2, np.full(1, 2.0))

    def test_uses_material_id_column_instead_of_row_position(self):
        FS = self._fs(
            [
                [20, 4.0, 1.0, 2.0, 3.0],
                [10, 8.0, 0.0, 0.0, 6.0],
            ],
            [10, 20],
        )
        FS["csr"]["tCurrent"] = np.array([0, 2])
        FS["csr"]["tNext"] = np.array([1, 3])
        T_all = np.array(
            [[1.0, 0.0], [0.0, 1.0], [2.0, 0.0], [0.0, 2.0]]
        )
        Q_corner = np.zeros((2, 2))

        kt1, kt2, kn1, kn2 = ksInterp(FS, T_all, Q_corner)

        np.testing.assert_allclose(kt1[0], [0.0, 0.0])
        np.testing.assert_allclose(kn1[0], [6.0, 8.0])
        np.testing.assert_allclose(kt1[1], [-2.0, 1.0])
        np.testing.assert_allclose(kn1[1], [3.0, 4.0])
        self.assertEqual(kt2.shape, (2,))
        self.assertEqual(kn2.shape, (2,))

    def test_rejects_two_phase_path(self):
        FS = self._fs([[1, 1.0, 0.0, 0.0, 1.0]], [1])
        FS["cfg"]["phasekey"] = 2

        with self.assertRaises(NotImplementedError):
            ksInterp(FS, np.ones((2, 2)), np.zeros((1, 2)))

    @staticmethod
    def _fs(kmap, material_ids):
        material_ids = np.asarray(material_ids)
        return {
            "cfg": {"phasekey": 1},
            "mesh": {
                "elem": np.column_stack(
                    (
                        np.zeros((material_ids.size, 4), dtype=int),
                        material_ids,
                    )
                )
            },
            "perm": {"kmap": np.asarray(kmap, dtype=float)},
            "csr": {
                "cornerElem": np.arange(material_ids.size),
                "tCurrent": np.zeros(material_ids.size, dtype=int),
                "tNext": np.ones(material_ids.size, dtype=int),
            },
        }


if __name__ == "__main__":
    unittest.main()
