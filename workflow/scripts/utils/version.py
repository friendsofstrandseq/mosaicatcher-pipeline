"""
Centralized version management for MosaiCatcher pipeline.

This module provides utilities to read the pipeline version from the centralized VERSION file.
"""

import os
from pathlib import Path


def get_version():
    """
    Read the pipeline version from the VERSION file at the repository root.

    Returns:
        str: The current pipeline version (e.g., "2.3.5")

    Raises:
        FileNotFoundError: If VERSION file cannot be found
    """
    # Try to find VERSION file starting from this script's directory and going up
    current_dir = Path(__file__).resolve().parent

    # Go up to repository root (workflow/scripts/utils -> workflow -> root)
    repo_root = current_dir.parent.parent.parent
    version_file = repo_root / "VERSION"

    if not version_file.exists():
        raise FileNotFoundError(
            f"VERSION file not found at {version_file}. "
            "Ensure VERSION file exists at repository root."
        )

    with open(version_file, "r") as f:
        version = f.read().strip()

    return version


__version__ = get_version()


if __name__ == "__main__":
    print(f"MosaiCatcher Pipeline v{__version__}")
