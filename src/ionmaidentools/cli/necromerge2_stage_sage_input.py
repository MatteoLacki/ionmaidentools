#!/usr/bin/env python3
"""CLI: necromerge2-stage-sage-input — assemble a Sage-ready spectra directory.

Absorbs the body of the old Snakemake `make_plain_pipeline_sage_input` (and its
tof-filtered `use rule ... with` variant): stages pmsms + tof2mz alongside a
precursors parquet (converted from the precursors mmappet) into one output
directory, which Sage reads as its `spectra` argument.
"""
import argparse
import subprocess
from pathlib import Path

from ionmaidentools.stage_dir import stage


def stage_sage_input(pmsms_dir: Path, tof2mz_dir: Path, precursors_dir: Path, output_dir: Path, python_bin: Path) -> None:
    out = Path(output_dir)
    stage(out, {"pmsms.mmappet": str(pmsms_dir), "tof2mz.mmappet": str(tof2mz_dir)})
    cmd = (
        f"{python_bin} -c \"import mmappet; mmappet.open_dataset('{precursors_dir}')"
        f".to_parquet('{out}/precursors.parquet', index=False)\""
    )
    subprocess.run(cmd, shell=True, check=True)


def main():
    p = argparse.ArgumentParser(
        description=(
            "Stage a pmsms dataset, tof2mz dataset, and precursors mmappet (converted "
            "to parquet) into one directory, ready to pass to Sage as its spectra input."
        )
    )
    p.add_argument("pmsms_dir", type=Path, help="Path to the pmsms.mmappet dataset directory.")
    p.add_argument("tof2mz_dir", type=Path, help="Path to the tof2mz.mmappet dataset directory.")
    p.add_argument("precursors_dir", type=Path, help="Path to the precursors mmappet dataset directory.")
    p.add_argument("output_dir", type=Path, help="Output directory to stage the Sage spectra input into.")
    p.add_argument("--python-bin", type=Path, default=Path("venvs/common/bin/python"), help="Python interpreter with mmappet installed, used for the parquet conversion.")
    args = p.parse_args()
    stage_sage_input(args.pmsms_dir, args.tof2mz_dir, args.precursors_dir, args.output_dir, args.python_bin)


if __name__ == "__main__":
    main()
