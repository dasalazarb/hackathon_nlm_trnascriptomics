"""Resolve pipeline outputs without ever consulting INPUT_DIR."""
from __future__ import annotations

import os
from pathlib import Path


class OutputDirectoryError(RuntimeError):
    pass


def resolve_output_dir() -> Path:
    explicit = os.environ.get("OUTPUT_DIR")
    if explicit:
        path = Path(explicit).expanduser()
        if path.is_dir():
            return path.resolve()
        raise OutputDirectoryError(f"OUTPUT_DIR does not exist or is not a directory: {path}")
    root = os.environ.get("PROJECT_ROOT")
    if root and (Path(root).expanduser() / "output").is_dir():
        return (Path(root).expanduser() / "output").resolve()
    raise OutputDirectoryError(
        "Pipeline outputs were not found. Set OUTPUT_DIR to the completed R pipeline output directory."
    )


def run_label(path: Path) -> str:
    """Return only a non-sensitive directory name, never a full server path."""
    return path.name or "output"
