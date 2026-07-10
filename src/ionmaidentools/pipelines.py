"""necroflow migration of necromerge2's Snakemake `short_test` chain.

First-pass scope only: tof-filtered branch, through `sage_summarize`. Plain-MGF/mzML,
FragPipe, `studies.smk`, and the `configs.smk` config generator are not ported.

Config fields travel as explicit named kwargs, one per fixed config field, straight from
the job TOML's own parsed dict -- no config file, no blob -- except where a tool's own
config table is already shaped like the file the tool wants to read (Sage's JSON config,
timstofu's scale-estimation TOML): there, necroflow's built-in `text_file` rule writes
the table straight to a cached node file instead of exploding it into a flag per field.
The external tools this pipeline shells out to (timstofu/quadops/boxing) were patched to
accept fields as flags directly, with their old config-file argument made optional (and
still positionally compatible with the old Snakemake pipeline) -- `ms1_find_argmaxes`/
`ms1_fit_scale_estimates` keep that optional TOML `settings_path` argument, which is what
`write_scale_estimation_config` now feeds. `git/ionmaidenmetal`'s `mkpmsms` C++ binary got
the same treatment the other direction: it previously had no config-file option at all
(`[pseudomsms]`'s algorithm-dependent flag set was hand-exploded via a `_cli_flags` dict-
to-flags helper, ported verbatim from the old Snakemake rule), so `--config <path.toml>`
was added to its shared `ParamTree` CLI parser (`param_tree.hpp`) -- `write_pseudomsms_config`
now feeds it `cfg.pseudomsms` directly, same as the Python tools, and `_cli_flags` is gone.
Every rule invokes a real installed CLI directly (no custom Python wrapper CLIs) except
`necromerge2-run-sage`, which runs Sage's compiled binary and writes a run_info.json
sidecar (its own JSON config file is generated separately by `write_sage_config` straight
from `cfg.sage`, which already matches Sage's config schema 1:1), and
`necromerge2-stage-sage-input`, which stages Sage's spectra input directory (a real
co-location requirement of Sage's own input format).

TODO(regression-db): the old Snakemake `sage_summarize`/`short_test` rules recorded
results into a SQLite regression DB and did an interactive baseline comparison. Neither
is ported here -- `sage_run_info`/`sage_summary` are the pipeline's terminal outputs.
"""
from __future__ import annotations

import json
import os
import shlex
import tomlkit

from dictodot import DotDict
from necroflow import NodeType, Pipeline, Rules
from pathlib import Path

R = Rules()

# Matches the old Snakemake rules' `threads: workflow.cores` -- these 4 tools
# (mkpmsms, precursor_neighbors_csr, tof_filter, score_based_pmsms_filter) are each
# meant to claim the whole machine while they run. `{threads}` in their command
# strings resolves from this constraint automatically (necroflow's built-in
# constraint-placeholder mechanism), so it no longer needs to travel through job
# config.
CORES = os.cpu_count() or 1


# --- source node types ---
class BrukerD(NodeType):
    filename = "input.d"


class Fasta(NodeType):
    filename = "fasta.fasta"


class MmappetDataset(NodeType):
    """Base type for outputs that are mmappet directories (see CLAUDE.md's
    "Precursor table format" convention). No `filename` of its own -- every
    concrete mmappet output subclasses this and sets its own."""


# --- compute artifact node types ---
class Ms1Events(NodeType):
    filename = "events.ms1"


class Ms2Events(NodeType):
    filename = "events.ms2"


class Tof2Mz(MmappetDataset):
    filename = "tof2mz.mmappet"


class ScaleEstimationConfig(NodeType):
    filename = "scale_estimation_config.toml"


class ArgmaxSample(MmappetDataset):
    filename = "argmax_sample.mmappet"


class ArgmaxSieveStats(NodeType):
    filename = "sieve_stats.toml"


class SampleTensors(MmappetDataset):
    filename = "sample_tensors.mmappet"


class ScaleEstimates(NodeType):
    filename = "scale_estimates"


class PrecursorCandidateSelectionConfig(NodeType):
    filename = "precursor_candidate_selection_config.toml"


class RawPrecursorClusters(MmappetDataset):
    filename = "raw_precursor_clusters.mmappet"


class PostprocessingConfig(NodeType):
    filename = "postprocessing_config.toml"


class PostprocessedPrecursorClusters(MmappetDataset):
    filename = "postprocessed_precursor_clusters.mmappet"


class PrecursorTransmissionConfig(NodeType):
    filename = "precursor_transmission_config.toml"


class TransmittedMs1Events(MmappetDataset):
    filename = "transmitted_ms1events.mmappet"


class TransmittedPrecursorClusters(MmappetDataset):
    filename = "transmitted_precursors.mmappet"


class FirstFilterPrecursors(MmappetDataset):
    filename = "first_filter_precursors.mmappet"


class PseudomsmsConfig(NodeType):
    filename = "pseudomsms_config.toml"


class Pmsms(MmappetDataset):
    filename = "pmsms.mmappet"


class MkpmsmsStats(NodeType):
    filename = "ms2peakpicking_stats.txt"


class Ms2IndexedPrecursors(MmappetDataset):
    filename = "ms2indexed_precursors.mmappet"


class PreSageFilteredPrecursors(MmappetDataset):
    filename = "pre_sage_filtered_precursors.mmappet"


class PrecursorNeighborsConfig(NodeType):
    filename = "precursor_neighbors_config.toml"


class PrecursorGridIndex(NodeType):
    filename = "precursor_grid_index"


class PrecursorNeighborsCsr(NodeType):
    filename = "precursor_neighbors_csr"


class NeighborScore(MmappetDataset):
    filename = "neighbor_score.mmappet"


class TofFilteredPmsms(MmappetDataset):
    filename = "tof_filtered_pmsms.mmappet"


class TofFilteredPrecursors(MmappetDataset):
    filename = "tof_filtered_precursors.mmappet"


class SageInputStaged(NodeType):
    filename = "sage_input.pmsms"


class SageConfig(NodeType):
    filename = "sage_config.json"


class SageRawOutdir(NodeType):
    filename = "sage_raw_outdir"


class SageRunInfo(NodeType):
    filename = "run_info.json"


class SageSummary(NodeType):
    filename = "results.sage.summary.tsv"


# --- source rules (symlink pre-existing files/dirs, no validation) ---
@R.command("ln -s $(realpath {path}) {tdf}")
def source_bruker_d(path: str):
    return BrukerD[tdf]


@R.command("ln -s $(realpath {path}) {fasta}")
def source_fasta(path: str):
    return Fasta[fasta]


# --- compute rules ---
@R.command(
    "venvs/common/bin/d2ms1 {tdf} {ms1}"
    " && test -f {ms1}/tof_row_starts.dat"
    " && test -f {ms1}/tof_urt_diff_index.dat"
    " && test -f {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
)
def tdf2ms1(tdf: BrukerD):
    return Ms1Events[ms1]


@R.command("git/ionmaidenmetal/build/tdf2ms ms2 {tdf} {ms2} --overwrite")
def tdf2ms2(tdf: BrukerD):
    return Ms2Events[ms2]


@R.command("venvs/common/bin/python scripts/tdf2tof2mz.py {tdf} {ms2} {tof2mz}")
def tdf2tof2mz(tdf: BrukerD, ms2: Ms2Events):
    return Tof2Mz[tof2mz]


R.text_file("write_scale_estimation_config", ScaleEstimationConfig)


@R.command("venvs/common/bin/ms1_find_argmaxes {ms1} {config} {argmaxes} {stats} --dataset_name {dataset}")
def find_ms1_argmaxes(ms1: Ms1Events, config: ScaleEstimationConfig, dataset: str):
    return ArgmaxSample[argmaxes], ArgmaxSieveStats[stats]


@R.command(
    "venvs/common/bin/ms1_extract_sample_tensors {ms1} {argmaxes} {tensors}"
    " --radii-tof {radii_tof} --radii-urt {radii_urt} --radii-scan {radii_scan}"
)
def extract_ms1_sample_tensors(
    ms1: Ms1Events, argmaxes: ArgmaxSample,
    radii_tof: int, radii_urt: int, radii_scan: int,
):
    return SampleTensors[tensors]


@R.command(
    "venvs/common/bin/ms1_fit_scale_estimates {argmaxes} {stats} {tensors} {config} {scales}"
    " --dataset_name {dataset}"
)
def fit_ms1_scale_estimates(
    argmaxes: ArgmaxSample, stats: ArgmaxSieveStats, tensors: SampleTensors,
    config: ScaleEstimationConfig, dataset: str,
):
    return ScaleEstimates[scales]


R.text_file("write_precursor_candidate_selection_config", PrecursorCandidateSelectionConfig)


@R.command("venvs/common/bin/ms1_select_candidates {ms1} {scale_estimates} {config} {clusters} --dataset_name {dataset}")
def select_precursor_candidates(
    ms1: Ms1Events, scale_estimates: ScaleEstimates, config: PrecursorCandidateSelectionConfig, dataset: str,
):
    return RawPrecursorClusters[clusters]


R.text_file("write_postprocessing_config", PostprocessingConfig)


@R.command(
    "venvs/common/bin/ms1_postprocess_candidates {tdf} {ms1} {candidates} {scale_estimates} {config} {clusters}"
    " --dataset_name {dataset}"
)
def postprocess_precursor_candidates(
    tdf: BrukerD, ms1: Ms1Events, candidates: RawPrecursorClusters, scale_estimates: ScaleEstimates,
    config: PostprocessingConfig, dataset: str,
):
    return PostprocessedPrecursorClusters[clusters]


R.text_file("write_precursor_transmission_config", PrecursorTransmissionConfig)


@R.command(
    "venvs/common/bin/transmit_precursors {tdf} {clusters} {config} {transpec}"
    " --output-precursors {precursors} --verbose"
    " && test -f {transpec}/schema.txt"
)
def transmit_precursors_into_fragment_space(
    tdf: BrukerD, clusters: PostprocessedPrecursorClusters, config: PrecursorTransmissionConfig,
):
    return TransmittedMs1Events[transpec], TransmittedPrecursorClusters[precursors]


@R.command("venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}")
def filter_first_precursors(precursors: TransmittedPrecursorClusters, filter: str):
    return FirstFilterPrecursors[filtered]


R.text_file("write_pseudomsms_config", PseudomsmsConfig)


@R.command(
    "git/ionmaidenmetal/build/mkpmsms --fragments {ms2} --transmitted-precursors {transprec}"
    " --precursors {filter_mm} --output {pmsms} --config {config}"
    " --threads {threads} --batch 1024"
    " && test -f {pmsms}/schema.txt && test -f {pmsms}/0.bin && test -f {pmsms}/1.bin && test -f {pmsms}/2.bin"
    " && test -f {pmsms}/dataindex.mmappet/schema.txt"
    " && test -f {pmsms}/dataindex.mmappet/0.bin && test -f {pmsms}/dataindex.mmappet/1.bin"
    " && test -f {pmsms}/dataindex.mmappet/2.bin"
    # peak-picking stats is a pure stdout diagnostic (no consumer downstream, so it'd be
    # pruned as a dead-end if it were its own rule) -- captured in this node's own
    # .rip/job.log automatically, same as it landed in the Snakemake job log before.
    " && venvs/common/bin/plot_ms2peakpicking_stats {pmsms}",
    threads=CORES,
)
def run_mkpmsms_binary(
    ms2: Ms2Events, transprec: TransmittedMs1Events, filter_mm: FirstFilterPrecursors,
    config: PseudomsmsConfig,
):
    return Pmsms[pmsms]


@R.command("venvs/common/bin/cut_and_index_precursors {filter_mm} {pmsms}/dataindex.mmappet {precursors}")
def cut_and_index_precursors(filter_mm: FirstFilterPrecursors, pmsms: Pmsms):
    return Ms2IndexedPrecursors[precursors]


@R.command("venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}")
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, filter: str):
    return PreSageFilteredPrecursors[filtered]


R.text_file("write_precursor_neighbors_config", PrecursorNeighborsConfig)


@R.command("venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid} --config {config}")
def build_precursor_grid_index(precursors: PreSageFilteredPrecursors, tdf: BrukerD, config: PrecursorNeighborsConfig):
    return PrecursorGridIndex[grid]


@R.command(
    "git/ionmaidenmetal/build/precursor_neighbors_csr"
    " --boxes-input {grid_index}/boxes.mmappet --index-input {grid_index} --output {csr}"
    " $(venvs/common/bin/precursor-neighbors-params {config} {tdf})"
    " --n-threads {threads}",
    threads=CORES,
)
def compute_precursor_neighbors(
    grid_index: PrecursorGridIndex, tdf: BrukerD, config: PrecursorNeighborsConfig,
):
    return PrecursorNeighborsCsr[csr]


@R.command(
    "git/ionmaidenmetal/build/tof_filter --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
    " --out-path {score} --n-threads {threads}",
    threads=CORES,
)
def tof_score_filter(pmsms: Pmsms, neighbors_csr: PrecursorNeighborsCsr):
    return NeighborScore[score]   # no config arg -- confirmed vestigial in the Snakemake rule


@R.command(
    "venvs/common/bin/python -m timstofu.cli.score_based_pmsms_filter"
    " {pmsms} {precursors} {neighbor_score} {pmsms_out} {precursors_out}"
    " --threads {threads} --score-margin {score_margin}",
    threads=CORES,
)
def materialize_tof_filtered_pmsms(
    pmsms: Pmsms, precursors: PreSageFilteredPrecursors, neighbor_score: NeighborScore,
    score_margin: int | float,
):
    return TofFilteredPmsms[pmsms_out], TofFilteredPrecursors[precursors_out]


@R.command(".venv/bin/necromerge2-stage-sage-input {pmsms} {tof2mz} {precursors} {staged}")
def make_tof_filtered_sage_input(pmsms: TofFilteredPmsms, tof2mz: Tof2Mz, precursors: TofFilteredPrecursors):
    return SageInputStaged[staged]


R.text_file("write_sage_config", SageConfig)


@R.command(".venv/bin/necromerge2-run-sage {spectra} {fasta} {sage_config} {outdir} {run_info}")
def run_sage(spectra: SageInputStaged, fasta: Fasta, sage_config: SageConfig):
    return SageRawOutdir[outdir], SageRunInfo[run_info]


@R.command("venvs/common/bin/sage-summarize-raw {sage_dir}/results.sage.tsv {summary} --fdr 0.01")
def sage_summarize(sage_dir: SageRawOutdir):
    return SageSummary[summary]


def sage_pipeline(cfg: dict) -> Pipeline:
    """tof-filtered Sage search chain reproducing the old short_test Snakemake target."""
    cfg = DotDict.Recursive(cfg)
    P = Pipeline()
    dataset = Path(cfg.tdf_path).stem

    P.tdf = R.source_bruker_d(path=cfg.tdf_path)
    P.fasta = R.source_fasta(path=cfg.fasta_path)

    P.ms1_events = R.tdf2ms1(P.tdf)
    P.ms2_events = R.tdf2ms2(P.tdf)
    P.tof2mz = R.tdf2tof2mz(P.tdf, P.ms2_events)

    se = cfg.scale_estimation
    P.scale_estimation_config = R.write_scale_estimation_config(text=tomlkit.dumps(se))
    P.argmaxes, P.argmax_sieve_stats = R.find_ms1_argmaxes(P.ms1_events, P.scale_estimation_config, dataset=dataset)
    P.sample_tensors = R.extract_ms1_sample_tensors(
        P.ms1_events, P.argmaxes,
        radii_tof=se.radii.tof, radii_urt=se.radii.urt, radii_scan=se.radii.scan,
    )
    P.scale_estimates = R.fit_ms1_scale_estimates(
        P.argmaxes, P.argmax_sieve_stats, P.sample_tensors, P.scale_estimation_config, dataset=dataset,
    )

    P.precursor_candidate_selection_config = R.write_precursor_candidate_selection_config(
        text=tomlkit.dumps(cfg.precursor_candidate_selection)
    )
    P.raw_precursor_clusters = R.select_precursor_candidates(
        P.ms1_events, P.scale_estimates, P.precursor_candidate_selection_config, dataset=dataset,
    )

    P.postprocessing_config = R.write_postprocessing_config(
        text=tomlkit.dumps(cfg.postprocessing_of_precursors)
    )
    P.postprocessed_precursor_clusters = R.postprocess_precursor_candidates(
        P.tdf, P.ms1_events, P.raw_precursor_clusters, P.scale_estimates, P.postprocessing_config, dataset=dataset,
    )

    P.precursor_transmission_config = R.write_precursor_transmission_config(
        text=tomlkit.dumps(cfg.precursor_transmission)
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = R.transmit_precursors_into_fragment_space(
        P.tdf, P.postprocessed_precursor_clusters, P.precursor_transmission_config,
    )

    P.first_filter_precursors = R.filter_first_precursors(
        P.transmitted_precursor_clusters,
        filter=shlex.quote(cfg.precursor_filters.mkpmsms.get("filter", "")),
    )

    P.pseudomsms_config = R.write_pseudomsms_config(text=tomlkit.dumps(cfg.pseudomsms))
    P.pmsms = R.run_mkpmsms_binary(
        P.ms2_events, P.transmitted_ms1events, P.first_filter_precursors, P.pseudomsms_config,
    )
    P.mkpmsms_stats = R.plot_mkpmsms_stats(P.pmsms)
    P.ms2indexed_precursors = R.cut_and_index_precursors(P.first_filter_precursors, P.pmsms)

    P.pre_sage_filtered_precursors = R.filter_pre_sage_precursors(
        P.ms2indexed_precursors,
        filter=shlex.quote(cfg.precursor_filters.pre_sage.get("filter", "")),
    )

    P.precursor_neighbors_config = R.write_precursor_neighbors_config(
        text=tomlkit.dumps(cfg.precursor_neighbors)
    )
    P.precursor_grid_index = R.build_precursor_grid_index(
        P.pre_sage_filtered_precursors, P.tdf, P.precursor_neighbors_config,
    )
    P.precursor_neighbors_csr = R.compute_precursor_neighbors(
        P.precursor_grid_index, P.tdf, P.precursor_neighbors_config,
    )
    P.neighbor_score = R.tof_score_filter(P.pmsms, P.precursor_neighbors_csr)

    P.tof_filtered_pmsms, P.tof_filtered_precursors = R.materialize_tof_filtered_pmsms(
        P.pmsms, P.pre_sage_filtered_precursors, P.neighbor_score,
        score_margin=cfg.tof_score_filter.score_margin,
    )
    P.sage_input = R.make_tof_filtered_sage_input(P.tof_filtered_pmsms, P.tof2mz, P.tof_filtered_precursors)

    P.sage_config = R.write_sage_config(
        text=json.dumps(cfg.sage, sort_keys=True, indent=2) + "\n"
    )
    P.sage_raw_outdir, P.sage_run_info = R.run_sage(P.sage_input, P.fasta, P.sage_config)
    P.sage_summary = R.sage_summarize(P.sage_raw_outdir)
    return P
