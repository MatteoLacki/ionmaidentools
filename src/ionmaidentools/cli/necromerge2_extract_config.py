#!/usr/bin/env python3
"""CLI: necromerge2-extract-config — pull one section out of the unified pipeline
config and materialize it as a standalone TOML or JSON file.

Format is chosen by out_path's suffix (.toml or .json). Keys are recursively
sorted at every nesting level before writing (array/list order is left
untouched) so a pure key-reorder edit of the unified config can never change a
downstream necroflow node's content hash. This performs no validation of
values — it is pure section-selection, format conversion, and canonicalization.
"""
import argparse
import json
import tomllib
from pathlib import Path

import tomlkit


def _get_section(data: dict, dotted_path: str):
    node = data
    for part in dotted_path.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"section {dotted_path!r} not found in config (missing {part!r})")
        node = node[part]
    return node


def _sorted(value):
    if isinstance(value, dict):
        return {k: _sorted(value[k]) for k in sorted(value)}
    if isinstance(value, list):
        return [_sorted(v) for v in value]
    return value


def extract_config(config_path: Path, section: str, out_path: Path) -> None:
    with open(config_path, "rb") as f:
        data = tomllib.load(f)
    section_data = _sorted(_get_section(data, section))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.suffix == ".json":
        out_path.write_text(json.dumps(section_data, indent=2) + "\n")
    else:
        out_path.write_text(tomlkit.dumps(section_data))


def main():
    p = argparse.ArgumentParser(
        description=(
            "Extract a dotted section from the unified pipeline config and write it "
            "as a standalone TOML or JSON file (format chosen by out_path's suffix)."
        )
    )
    p.add_argument("config_path", type=Path, help="Path to the unified config TOML.")
    p.add_argument("section", type=str, help="Dotted section path, e.g. precursor_filters.mkpmsms")
    p.add_argument("out_path", type=Path, help="Output file path; .toml or .json.")
    args = p.parse_args()
    extract_config(args.config_path, args.section, args.out_path)


if __name__ == "__main__":
    main()
