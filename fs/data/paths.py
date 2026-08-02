"""Known data path resolver translated from ``+fs/+data/paths.m``."""

import os
from pathlib import Path
from typing import Optional, Union


_KNOWN_FILES = {
    "spe10": "spe10.mat",
    "spe_perm": "spe_perm.dat",
    "gmsh": "gmsh.exe",
}


def paths(key: str, root: Optional[Union[str, Path]] = None) -> str:
    """Return the first existing path for ``key``, or an empty string."""
    if key not in _KNOWN_FILES:
        valid = ", ".join(_KNOWN_FILES)
        raise ValueError(f"Unknown data key {key!r}; valid keys: {valid}")

    filename = _KNOWN_FILES[key]
    if root is None:
        root = os.environ.get("FS_TEST_ROOT", os.getcwd())
    candidates = [Path(root) / filename]

    home = Path.home()
    candidates.append(
        home
        / "projects"
        / "new-axon"
        / "axon"
        / "my-axon"
        / "dev-projects"
        / "flowsim-vectorize"
        / "data"
        / "legacy-binaries"
        / filename
    )

    data_directory = os.environ.get("FS_DATA_DIR")
    if data_directory:
        candidates.append(Path(data_directory) / filename)

    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return ""
