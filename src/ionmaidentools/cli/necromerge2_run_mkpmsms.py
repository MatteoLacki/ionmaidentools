#!/usr/bin/env python3
"""CLI: necromerge2-run-mkpmsms — run the mkpmsms C++ tool and its post-processing.

Absorbs the procedural body of the old Snakemake `mkpmsms` rule: reads the
`[pseudomsms]` section of a pipeline config to build the tool's CLI flags, runs the
mkpmsms binary, asserts its outputs exist, then chains the peak-picking stats plot
and the precursor cut-and-index step.
"""
import argparse
import subprocess
import tomllib
from pathlib import Path


def _to_str(x) -> str:
    if isinstance(x, list):
        return " ".join(map(str, x))
    s = str(x)
    return f'"{s}"' if " " in s else s


def run_mkpmsms(
    staged_dir: Path,
    config_path: Path,
    pmsms_out: Path,
    precursors_out: Path,
    mkpmsms_bin: Path,
    cut_and_index_bin: Path,
    plot_stats_bin: Path,
    threads: int,
    batch_size: int,
) -> None:
    staged = Path(staged_dir)
    ms2 = staged / "events.ms2"
    transprec = staged / "transmitted_precursors.mmappet"
    filter_mm = staged / "precursor_filter.mmappet"

    with open(config_path, "rb") as f:
        settings = tomllib.load(f)["pseudomsms"]
    method = settings["tofs_extraction_method"]
    kwargs = " ".join(
        f"--{name} {_to_str(param)}"
        for name, param in settings.get("tofs_extraction_params", {}).items()
    )

    pmsms_dir = Path(pmsms_out)
    cmd = (
        f"{mkpmsms_bin} --fragments {ms2} --transmitted-precursors {transprec}"
        f" --precursors {filter_mm} --output {pmsms_dir} --method {method}"
        f" --threads {threads} --batch {batch_size} {kwargs}"
    )
    subprocess.run(cmd, shell=True, check=True)

    for rel in (
        "schema.txt", "0.bin", "1.bin", "2.bin",
        "dataindex.mmappet", "dataindex.mmappet/schema.txt",
        "dataindex.mmappet/0.bin", "dataindex.mmappet/1.bin", "dataindex.mmappet/2.bin",
    ):
        if not (pmsms_dir / rel).exists():
            raise RuntimeError(f"expected mkpmsms output missing: {pmsms_dir / rel}")

    subprocess.run([str(plot_stats_bin), str(pmsms_dir)], check=True)
    subprocess.run(
        [str(cut_and_index_bin), str(filter_mm), str(pmsms_dir / "dataindex.mmappet"), str(precursors_out)],
        check=True,
    )


def main():
    p = argparse.ArgumentParser(
        description=(
            "Run mkpmsms against a staged input directory (events.ms2, "
            "transmitted_precursors.mmappet, precursor_filter.mmappet), then chain "
            "peak-picking stats plotting and precursor cut-and-index."
        )
    )
    p.add_argument("staged_dir", type=Path, help="Directory staged by necromerge2-stage-dir.")
    p.add_argument("config_path", type=Path, help="Pipeline config TOML containing a [pseudomsms] section.")
    p.add_argument("pmsms_out", type=Path, help="Output directory for the pmsms dataset.")
    p.add_argument("precursors_out", type=Path, help="Output path for the cut-and-indexed precursors mmappet.")
    p.add_argument("--mkpmsms-bin", type=Path, default=Path("git/ionmaidenmetal/build/mkpmsms"), help="Path to the mkpmsms binary.")
    p.add_argument("--cut-and-index-bin", type=Path, default=Path("venvs/common/bin/cut_and_index_precursors"), help="Path to the cut_and_index_precursors executable.")
    p.add_argument("--plot-stats-bin", type=Path, default=Path("venvs/common/bin/plot_ms2peakpicking_stats"), help="Path to the plot_ms2peakpicking_stats executable.")
    p.add_argument("--threads", type=int, default=1, help="Number of threads to pass to mkpmsms.")
    p.add_argument("--batch-size", type=int, default=1024, help="Batch size to pass to mkpmsms.")
    args = p.parse_args()
    run_mkpmsms(
        args.staged_dir, args.config_path, args.pmsms_out, args.precursors_out,
        args.mkpmsms_bin, args.cut_and_index_bin, args.plot_stats_bin,
        args.threads, args.batch_size,
    )


if __name__ == "__main__":
    main()
