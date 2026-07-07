#!/usr/bin/env python3
"""CLI: necromerge2-run-sage — run the Sage search engine and write a run_info sidecar.

Absorbs the procedural body of the old Snakemake `sage` rule: builds Sage's own JSON
config from explicit flags (Sage's compiled Rust binary genuinely requires a JSON file --
out of scope to rebuild it to accept flags directly, unlike the Python tools this
pipeline also shells out to), prints the Sage version, runs the search, asserts the
expected output files exist, and writes a run_info.json describing the invocation (kept
for future regression-DB wiring — see TODO in the necromerge2 necroflow pipeline module).
"""
import argparse
import json
import subprocess
from pathlib import Path


def _bool(s: str) -> bool:
    return s.lower() in ("true", "1", "yes")


def build_sage_config(args) -> dict:
    return {
        "database": {
            "bucket_size": args.database_bucket_size,
            "fragment_min_mz": args.database_fragment_min_mz,
            "fragment_max_mz": args.database_fragment_max_mz,
            "peptide_min_mass": args.database_peptide_min_mass,
            "peptide_max_mass": args.database_peptide_max_mass,
            "min_ion_index": args.database_min_ion_index,
            "enzyme": {
                "missed_cleavages": args.database_enzyme_missed_cleavages,
                "cleave_at": args.database_enzyme_cleave_at,
                "restrict": args.database_enzyme_restrict,
                "min_len": args.database_enzyme_min_len,
                "max_len": args.database_enzyme_max_len,
                "c_terminal": args.database_enzyme_c_terminal,
            },
            "static_mods": json.loads(args.static_mods_json),
            "variable_mods": json.loads(args.variable_mods_json),
            "max_variable_mods": args.database_max_variable_mods,
            "decoy_tag": args.database_decoy_tag,
            "generate_decoys": args.database_generate_decoys,
        },
        "deisotope": args.deisotope,
        "chimera": args.chimera,
        "wide_window": args.wide_window,
        "predict_rt": args.predict_rt,
        "min_peaks": args.min_peaks,
        "max_peaks": args.max_peaks,
        "min_matched_peaks": args.min_matched_peaks,
        "max_fragment_charge": args.max_fragment_charge,
        "ignore_precursor_charge": args.ignore_precursor_charge,
        "parallel": args.parallel,
        "report_psms": args.report_psms,
        "precursor_tol": {"ppm": [args.precursor_tol_ppm_lo, args.precursor_tol_ppm_hi]},
        "fragment_tol": {"ppm": [args.fragment_tol_ppm_lo, args.fragment_tol_ppm_hi]},
        "isotope_errors": [args.isotope_errors_lo, args.isotope_errors_hi],
    }


def run_sage(spectra_dir: Path, fasta_path: Path, config_path: Path, outdir: Path, run_info_path: Path, sage_bin: Path) -> None:
    subprocess.run([str(sage_bin), "--version"], check=True)
    cmd = (
        f"{sage_bin} -f {fasta_path} --annotate-matches --write-pin"
        f" --output_directory {outdir} {config_path} {spectra_dir}"
    )
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
    p.add_argument("outdir", type=Path, help="Output directory for raw Sage results.")
    p.add_argument("run_info_path", type=Path, help="Output path for the run_info.json sidecar.")
    p.add_argument("--sage-bin", type=Path, default=Path("software/sage/devel_fixed/sage"), help="Path to the sage executable.")

    p.add_argument("--deisotope", type=_bool, required=True)
    p.add_argument("--chimera", type=_bool, required=True)
    p.add_argument("--wide-window", type=_bool, required=True)
    p.add_argument("--predict-rt", type=_bool, required=True)
    p.add_argument("--min-peaks", type=int, required=True)
    p.add_argument("--max-peaks", type=int, required=True)
    p.add_argument("--min-matched-peaks", type=int, required=True)
    p.add_argument("--max-fragment-charge", type=int, required=True)
    p.add_argument("--ignore-precursor-charge", type=_bool, required=True)
    p.add_argument("--parallel", type=_bool, required=True)
    p.add_argument("--report-psms", type=int, required=True)
    p.add_argument("--isotope-errors-lo", type=int, required=True)
    p.add_argument("--isotope-errors-hi", type=int, required=True)
    p.add_argument("--database-bucket-size", type=int, required=True)
    p.add_argument("--database-fragment-min-mz", type=float, required=True)
    p.add_argument("--database-fragment-max-mz", type=float, required=True)
    p.add_argument("--database-peptide-min-mass", type=float, required=True)
    p.add_argument("--database-peptide-max-mass", type=float, required=True)
    p.add_argument("--database-min-ion-index", type=int, required=True)
    p.add_argument("--database-max-variable-mods", type=int, required=True)
    p.add_argument("--database-decoy-tag", type=str, required=True)
    p.add_argument("--database-generate-decoys", type=_bool, required=True)
    p.add_argument("--database-enzyme-missed-cleavages", type=int, required=True)
    p.add_argument("--database-enzyme-cleave-at", type=str, required=True)
    p.add_argument("--database-enzyme-restrict", type=str, required=True)
    p.add_argument("--database-enzyme-min-len", type=int, required=True)
    p.add_argument("--database-enzyme-max-len", type=int, required=True)
    p.add_argument("--database-enzyme-c-terminal", type=_bool, required=True)
    p.add_argument("--precursor-tol-ppm-lo", type=float, required=True)
    p.add_argument("--precursor-tol-ppm-hi", type=float, required=True)
    p.add_argument("--fragment-tol-ppm-lo", type=float, required=True)
    p.add_argument("--fragment-tol-ppm-hi", type=float, required=True)
    p.add_argument("--static-mods-json", type=str, required=True)
    p.add_argument("--variable-mods-json", type=str, required=True)
    args = p.parse_args()

    sage_config = build_sage_config(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    config_path = args.outdir / "sage_config.json"
    config_path.write_text(json.dumps(sage_config, indent=2) + "\n")

    run_sage(args.spectra_dir, args.fasta_path, config_path, args.outdir, args.run_info_path, args.sage_bin)


if __name__ == "__main__":
    main()
