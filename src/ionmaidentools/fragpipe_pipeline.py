"""Necroflow FragPipe workflows: the main ionmaiden pipeline followed by a FragPipe search
of its search mzML (`ionmaiden_fragpipe_pipeline`, for jobs with `[fragpipe]`), and the
FragPipe smoke test on simulated peptides (`fragpipe_synthetic_pipeline`). Builds on
`pipelines.py`, which knows nothing of FragPipe.
"""

from __future__ import annotations

import tomlkit

from dictodot import DotDict
from necroflow import Node, Pipeline, Resource, command, outputs, shell, workflow

from ionmaidentools.pipelines import (
    CORES,
    Fasta,
    Mgf,
    MmappetDataset,
    Pmsms,
    TofFilteredMzml,
    TofFilteredMzmlIdmap,
    build_ionmaiden_pipeline,
    source_fasta,
    write_pipeline_config,
)


class FragpipeWorkflow(Node):
    """A hand-maintained FragPipe `.workflow` file -- treated as a data file
    like Fasta/BrukerD, symlinked in as-is. Never generated or patched by
    the pipeline (unlike Sage's JSON config)."""

    filename = "fragpipe_workflow.workflow"


class FragpipeManifest(Node):
    filename = "fragpipe_manifest.tsv"


class FragpipeResultsDir(Node):
    filename = "fragpipe_results"


class FragpipeLog(Node):
    filename = "fragpipe_full_log.txt"


class FragpipeSummary(Node):
    filename = "fragpipe_summary.txt"


class FragpipeDecoyFasta(Node):
    """Target+decoy FASTA via FragPipe's own Philosopher tool (`database
    --custom ... --prefix rev_`) -- unlike Sage, FragPipe/Philosopher has no
    generate-decoys-at-search-time option, so the database must already
    contain them. Standalone output, not wired into `run_fragpipe`: the
    `.workflow` file stays a plain, hand-maintained data file (see
    `FragpipeWorkflow`), so pointing `database.db-path=` at this file is the
    user's responsibility, same as everything else in that file."""

    filename = "decoy_database.fas"


class SyntheticPmsms(Pmsms):
    """Fragment ions for peptides simulated from a FASTA (Koina-predicted, no
    real acquisition) -- see scripts/simulate_peptides_to_pmsms.py. mz is
    baked in directly, with no `tof` (unlike real pmsms, whose m/z every
    consumer reads through a Tof2MzTable)."""

    filename = "synthetic_pmsms.mmappet"


class SyntheticPrecursors(MmappetDataset):
    filename = "synthetic_precursors.mmappet"


class SyntheticMzml(TofFilteredMzml):
    """Accepted wherever TofFilteredMzml is, e.g. write_fragpipe_manifest's
    mzml input -- same file shape, different (simulated) origin."""


class SyntheticMzmlIdmap(TofFilteredMzmlIdmap):
    pass


class SyntheticMgf(Mgf):
    pass


@command
def source_fragpipe_workflow(path: str):
    workflow = outputs(FragpipeWorkflow)
    shell("ln -s $(realpath {path}) {workflow}")
    return workflow


@command
def generate_fragpipe_decoy_fasta(fasta: Fasta):
    decoy_fasta = outputs(FragpipeDecoyFasta)
    # `cd` changes the meaning of every relative path substituted into this
    # template, including the philosopher binary's own repo-relative path and
    # {decoy_fasta} itself -- so all of them are frozen to absolute via shell
    # variables *before* the cd, same reasoning as source_fasta's realpath.
    shell(
        "REPO_ROOT=$(pwd) && FASTA_ABS=$(realpath {fasta})"
        " && DECOY_ABS=$(realpath -m {decoy_fasta})"
        " && cd {workdir}"
        ' && "$REPO_ROOT/software/fragpipe/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9"'
        " workspace --init --nocheck"
        ' && "$REPO_ROOT/software/fragpipe/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9"'
        ' database --custom "$FASTA_ABS" --prefix rev_'
        ' && mv $(ls *decoys*.fas | head -n 1) "$DECOY_ABS"'
        " && rm -rf .meta"
        ' && test -f "$DECOY_ABS"'
    )
    return decoy_fasta


@command
def patch_fragpipe_workflow(workflow: FragpipeWorkflow, decoy_fasta: FragpipeDecoyFasta):
    patched_workflow = outputs(FragpipeWorkflow)
    shell(
        "sed -e 's|^database\\.db-path=.*|database.db-path={decoy_fasta}|'"
        " {workflow} > {patched_workflow}"
    )
    return patched_workflow


@command
def simulate_peptides_to_pmsms(
    fasta: Fasta,
    charges: str,
    max_peptides_per_protein: int,
    seed: int,
):
    pmsms, precursors = outputs(SyntheticPmsms, SyntheticPrecursors)
    shell(
        "venvs/common/bin/python scripts/simulate_peptides_to_pmsms.py"
        " {fasta} {pmsms} {precursors}"
        " --charges {charges} --max-peptides-per-protein {max_peptides_per_protein}"
        " --seed {seed}"
        " && test -f {pmsms}/schema.txt && test -f {pmsms}/dataindex.mmappet/schema.txt"
        " && test -f {precursors}/schema.txt"
    )
    return pmsms, precursors


@command
def convert_synthetic_pmsms_to_mzml(
    pmsms: SyntheticPmsms, precursors: SyntheticPrecursors,
    threads: Resource[int] = CORES,
):
    mzml, idmap = outputs(SyntheticMzml, SyntheticMzmlIdmap)
    shell(
        "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
        " --threads {threads} --numpress --zlib-level 9"
        " && test -f {mzml} && test -f {idmap}/schema.txt"
    )
    return mzml, idmap


@command
def convert_synthetic_pmsms_to_mgf(
    pmsms: SyntheticPmsms,
    precursors: SyntheticPrecursors,
    config_path: str,
):
    mgf = outputs(SyntheticMgf)
    shell(
        "venvs/common/bin/msms2mgf {pmsms} {precursors} {config_path} {mgf}"
        " && test -f {mgf}"
    )
    return mgf


@command
def write_fragpipe_manifest(mzml: TofFilteredMzml):
    manifest = outputs(FragpipeManifest)
    shell('printf "%s\\tA\\t1\\tDDA" "$(realpath {mzml})" > {manifest}')
    return manifest


@command
def run_fragpipe(
    manifest: FragpipeManifest,
    workflow: FragpipeWorkflow,
    ram: int = 0,
    threads: Resource[int] = CORES,
):
    dir = outputs(FragpipeResultsDir)
    shell(
        "software/fragpipe/fragpipe-24.0/bin/fragpipe --headless"
        " --workflow {workflow} --manifest {manifest} --workdir {dir}"
        " --config-tools-folder software/fragpipe/fragpipe-24.0/tools"
        " --threads {threads} --ram {ram}"
        ' && test -n "$(ls {dir}/log_*.txt 2>/dev/null)"'
    )
    return dir


@command
def extract_fragpipe_log(dir: FragpipeResultsDir):
    log = outputs(FragpipeLog)
    shell("cp $(ls {dir}/log_*.txt | head -n 1) {log}")
    return log


@command
def summarize_fragpipe(log: FragpipeLog):
    summary = outputs(FragpipeSummary)
    shell(
        "grep -m 1 temp {log} > {summary}"
        " && grep -n 'MASS CALIBRATION' {log} -A 8 >> {summary}"
        " && grep 'Final report numbers after FDR filtering, and post-processing' {log} >> {summary}"
    )
    return summary

def fragpipe_search(P: Pipeline, fasta, mzml, fragpipe_cfg: DotDict) -> None:
    """FragPipe on the `mzml` Node against `fasta` plus Philosopher's decoys, labelled
    `fragpipe_*` on `P`, inside an active workflow. `fragpipe_cfg` is the job's
    `[fragpipe]`: `workflow_path`, optional `ram` (FragPipe's `--ram`, default 0)."""
    P.fragpipe_workflow = source_fragpipe_workflow(path=fragpipe_cfg.workflow_path)
    P.fragpipe_decoy_fasta = generate_fragpipe_decoy_fasta(fasta)
    P.fragpipe_workflow_patched = patch_fragpipe_workflow(
        P.fragpipe_workflow, P.fragpipe_decoy_fasta
    )
    P.fragpipe_manifest = write_fragpipe_manifest(mzml)
    P.fragpipe_results_dir = run_fragpipe(
        P.fragpipe_manifest,
        P.fragpipe_workflow_patched,
        ram=fragpipe_cfg.get("ram", 0),
    )
    P.fragpipe_log = extract_fragpipe_log(P.fragpipe_results_dir)
    P.fragpipe_summary = summarize_fragpipe(P.fragpipe_log)


@workflow
def ionmaiden_fragpipe_pipeline(P: Pipeline, config: dict) -> None:
    """The main IonMaiden pipeline, then FragPipe on its search mzML (`[fragpipe]`
    required) -- a side-by-side comparison with SAGE on the same spectra."""
    build_ionmaiden_pipeline(P, config)
    fragpipe_search(P, P.fasta, P.search_mzml, DotDict.Recursive(config).fragpipe)


@workflow
def fragpipe_synthetic_pipeline(P: Pipeline, config: dict) -> None:
    """FragPipe smoke test on simulated data -- no Bruker .d input, no Sage.

    1. Simulate tryptic peptides from cfg.fasta_path (a small multi-protein FASTA)
       via Koina-predicted fragment ions, written as a pmsms dataset
       (`simulate_peptides_to_pmsms`).
    2. Convert that pmsms to a spectrum file: mzML (`convert_synthetic_pmsms_to_mzml`,
       via `git/pmsms2mzml`) or MGF (`convert_synthetic_pmsms_to_mgf`, via
       `venvs/common/bin/msms2mgf` + `[mgf].config_path`), picked by
       cfg.output_format ("mzml", the default, or "mgf").
    3. mzML only: run FragPipe on it (`fragpipe_search`, the same search
       `ionmaiden_fragpipe_pipeline` runs).
       MGF has no search engine wired to it here (msms2mgf's plain-MGF output isn't
       what run_fragpipe's manifest step expects) -- P.synthetic_mgf is terminal.

    Verified end to end against a real FragPipe 24.0 + Philosopher install with a
    10-protein/~200-peptide synthetic FASTA: both convert_synthetic_pmsms_to_mzml
    (-> FragPipe, all 10 proteins recovered) and convert_synthetic_pmsms_to_mgf
    produced valid output from the same simulate_peptides_to_pmsms run.
    """
    cfg: DotDict = DotDict.Recursive(config)
    P.pipeline_config = write_pipeline_config(text=tomlkit.dumps(config))

    # Peptide Simulation
    P.fasta = source_fasta(path=cfg.fasta_path)
    P.synthetic_pmsms, P.synthetic_precursors = simulate_peptides_to_pmsms(
        P.fasta,
        charges=cfg.get("simulation", {}).get("charges", "2,3"),
        max_peptides_per_protein=cfg.get("simulation", {}).get(
            "max_peptides_per_protein", 12
        ),
        seed=cfg.get("simulation", {}).get("seed", 20260715),
    )

    # Spectrum File Creation
    output_format = cfg.get("output_format", "mzml")
    if output_format == "mzml":
        P.synthetic_mzml, P.synthetic_mzml_idmap = convert_synthetic_pmsms_to_mzml(
            P.synthetic_pmsms,
            P.synthetic_precursors,
        )

        fragpipe_search(P, P.fasta, P.synthetic_mzml, cfg.fragpipe)
    elif output_format == "mgf":
        mgf_config_path = cfg.get("mgf", {}).get("config_path")
        if not mgf_config_path:
            raise ValueError(
                "cfg.mgf.config_path is required when cfg.output_format == 'mgf'"
            )
        P.synthetic_mgf = convert_synthetic_pmsms_to_mgf(
            P.synthetic_pmsms,
            P.synthetic_precursors,
            config_path=mgf_config_path,
        )
    else:
        raise ValueError(
            f"cfg.output_format must be 'mzml' or 'mgf', got {output_format!r}"
        )
