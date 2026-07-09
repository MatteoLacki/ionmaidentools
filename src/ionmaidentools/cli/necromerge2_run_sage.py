#!/usr/bin/env python3
"""CLI: necromerge2-run-sage — run the Sage search engine and write a run_info sidecar.

Absorbs the procedural body of the old Snakemake `sage` rule: prints the Sage version,
runs the search against an already-materialized JSON config (written upstream by
necroflow's `write_sage_config` built-in text-file rule -- Sage's compiled Rust binary
genuinely requires a JSON file on disk, out of scope to rebuild it to accept flags
directly), asserts the expected output files exist, and writes a run_info.json
describing the invocation (kept for future regression-DB wiring — see TODO in the
necromerge2 necroflow pipeline module).
"""
import argparse
import json
import subprocess
from pathlib import Path


def run_sage(spectra_dir: Path, fasta_path: Path, config_path: Path, outdir: Path, run_info_path: Path, sage_bin: Path) -> None:
    subprocess.run([str(sage_bin), "--version"], check=True)
    cmd = (
        f"{sage_bin} -f {fasta_path} --annotate-matches --write-pin"
        f" --output_directory {outdir} {config_path} {spectra_dir}"
    )
    Path(outdir).mkdir(parents=True, exist_ok=True)
    subprocess.run(cmd, shell=True, check=True)

    out = Path(outdir)
    for rel in ("results.json", "results.sage.pin", "results.sage.tsv", "matched_fragments.sage.tsv"):
        if not (out / rel).exists():
            raise RuntimeError(f"expected sage output missing: {out / rel}")

    run_info = {
        "search_engine": "sage",
        "search_tool_call": cmd,
        "search_config_path": str(config_path),
    }
    Path(run_info_path).parent.mkdir(parents=True, exist_ok=True)
    Path(run_info_path).write_text(json.dumps(run_info, indent=2) + "\n")


def main():
    p = argparse.ArgumentParser(
        description="Run Sage against a staged spectra directory and record a run_info.json sidecar."
    )
    p.add_argument("spectra_dir", type=Path, help="Directory staged by necromerge2-stage-sage-input.")
    p.add_argument("fasta_path", type=Path, help="Path to the FASTA database.")
    p.add_argument("config_path", type=Path, help="Path to Sage's own JSON config file.")
    p.add_argument("outdir", type=Path, help="Output directory for raw Sage results.")
    p.add_argument("run_info_path", type=Path, help="Output path for the run_info.json sidecar.")
    p.add_argument("--sage-bin", type=Path, default=Path("software/sage/devel_fixed/sage"), help="Path to the sage executable.")
    args = p.parse_args()

    run_sage(args.spectra_dir, args.fasta_path, args.config_path, args.outdir, args.run_info_path, args.sage_bin)


if __name__ == "__main__":
    main()
