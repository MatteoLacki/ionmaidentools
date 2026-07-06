#!/usr/bin/env python3
"""CLI: necromerge2-stage-dir — stage a set of named inputs into a fresh directory.

Generic replacement for necromerge2's old Snakemake `_stage_dir`/`_stage_mkpmsms_input`
helpers: wipes/creates `output_dir`, then reflink-copies each `--link NAME=SRC`
target into `output_dir/NAME`.
"""
import argparse
from pathlib import Path

from ionmaidentools.stage_dir import stage


def _parse_link(spec: str) -> tuple[str, str]:
    if "=" not in spec:
        raise argparse.ArgumentTypeError(f"expected NAME=SRC_PATH, got {spec!r}")
    name, _, src = spec.partition("=")
    return name, src


def main():
    p = argparse.ArgumentParser(
        description=(
            "Stage a set of named inputs into a fresh directory, as reflink-if-possible "
            "copies. Replaces output_dir entirely (unlinks/removes it first)."
        )
    )
    p.add_argument("output_dir", type=Path, help="Directory to (re)create and stage into.")
    p.add_argument(
        "--link",
        action="append",
        dest="links",
        required=True,
        metavar="NAME=SRC_PATH",
        type=_parse_link,
        help="One staged entry: output_dir/NAME will be a copy of SRC_PATH. Repeatable.",
    )
    args = p.parse_args()
    stage(args.output_dir, dict(args.links))


if __name__ == "__main__":
    main()
