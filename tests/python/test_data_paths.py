import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fs.data import paths


class DataPathsTests(unittest.TestCase):
    def test_prefers_repository_root(self):
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as data:
            root_file = Path(root) / "spe10.mat"
            root_file.touch()
            (Path(data) / "spe10.mat").touch()

            with patch.dict(os.environ, {"FS_DATA_DIR": data}):
                result = paths("spe10", root=root)

            self.assertEqual(result, str(root_file))

    def test_uses_configured_data_directory(self):
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as data:
            data_file = Path(data) / "spe_perm.dat"
            data_file.touch()

            with patch.dict(os.environ, {"FS_DATA_DIR": data}):
                result = paths("spe_perm", root=root)

            self.assertEqual(result, str(data_file))

    def test_returns_empty_string_when_known_file_is_missing(self):
        with tempfile.TemporaryDirectory() as root:
            with patch.dict(os.environ, {"FS_DATA_DIR": ""}):
                self.assertEqual(paths("gmsh", root=root), "")

    def test_rejects_unknown_key(self):
        with self.assertRaisesRegex(ValueError, "valid keys"):
            paths("missing")


if __name__ == "__main__":
    unittest.main()
