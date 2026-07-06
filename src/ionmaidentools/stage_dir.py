"""Shared directory-staging helper: build a fresh directory of symlink-free copies.

Used by the `necromerge2-stage-dir` CLI, and importable directly by other
wrapper tools (e.g. `necromerge2-stage-sage-input`) that need to stage a few
named inputs into one directory before invoking a downstream tool.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def stage(output_dir: Path, links: dict[str, str]) -> None:
    """Replace `output_dir` with a fresh directory containing one entry per `links`.

    Each `name: target` pair is materialized as `output_dir/name`, a
    reflink-if-possible copy of `target` (resolved to its real path first, so
    this also works when `target` is itself a symlink).
    """
    out = Path(output_dir)
    if out.is_symlink() or out.is_file():
        out.unlink()
    elif out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    for name, target in links.items():
        resolved = Path(target).resolve()
        subprocess.run(
            ["cp", "-a", "--reflink=auto", str(resolved), str(out / name)],
            check=True,
        )
