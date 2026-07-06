"""necroflow migration of necromerge2's Snakemake `short_test` chain.

First-pass scope only: tof-filtered branch, through `sage_summarize`. Plain-MGF/mzML,
FragPipe, `studies.smk`, and the `configs.smk` config generator are not ported.

TODO(regression-db): the old Snakemake `sage_summarize`/`short_test` rules recorded
results into a SQLite regression DB and did an interactive baseline comparison. Neither
is ported here — `sage_run_info`/`sage_summary` are the pipeline's terminal outputs.
"""
from __future__ import annotations

from necroflow import NodeType, Pipeline, Rules

R = Rules()


def with_config(section: str, peg_output: str, cmd: str, ext: str = "toml") -> str:
    """Extract `section` from the unified config into a scratch file living in the
    private directory necroflow already allocates for this rule's own `peg_output`
    (any one of the rule's declared outputs — the directory is shared by all of them),
    then run `cmd`. `cmd` should reference the extracted file as the shell variable
    "$CFG" (a plain shell variable, not a necroflow placeholder — untouched by
    necroflow's own `{...}` substitution).
    """
    return (
        f'CFG="$(dirname {{{peg_output}}})/extracted_config.{ext}"'
        f' && .venv/bin/necromerge2-extract-config {{config}} {section} "$CFG"'
        f' && {cmd}'
    )


# --- source node types ---
class BrukerD(NodeType):
    filename = "input.d"


class Fasta(NodeType):
    filename = "fasta.fasta"


class UnifiedConfig(NodeType):
    filename = "config.toml"


# --- compute artifact node types ---
class Ms1Events(NodeType):
    filename = "events.ms1"


class Ms2Events(NodeType):
    filename = "events.ms2"


class Tof2Mz(NodeType):
    filename = "tof2mz.mmappet"


class ScaleEstimates(NodeType):
    filename = "scale_estimates"


class RawPrecursorClusters(NodeType):
    filename = "raw_precursor_clusters.mmappet"


class PostprocessedPrecursorClusters(NodeType):
    filename = "postprocessed_precursor_clusters.mmappet"


class TransmittedMs1Events(NodeType):
    filename = "transmitted_ms1events"


class TransmittedPrecursorClusters(NodeType):
    filename = "transmitted_precursors.mmappet"


class FirstFilterPrecursors(NodeType):
    filename = "first_filter_precursors.mmappet"


class MkpmsmsInputStaged(NodeType):
    filename = "mkpmsms_input"


class Pmsms(NodeType):
    filename = "pmsms.mmappet"


class Ms2IndexedPrecursors(NodeType):
    filename = "ms2indexed_precursors.mmappet"


class PreSageFilteredPrecursors(NodeType):
    filename = "pre_sage_filtered_precursors.mmappet"


class PrecursorGridIndex(NodeType):
    filename = "precursor_grid_index"


class PrecursorNeighborsCsr(NodeType):
    filename = "precursor_neighbors_csr"


class NeighborScore(NodeType):
    filename = "neighbor_score.mmappet"


class TofFilteredPmsms(NodeType):
    filename = "tof_filtered_pmsms"


class TofFilteredPrecursors(NodeType):
    filename = "tof_filtered_precursors.mmappet"


class SageInputStaged(NodeType):
    filename = "sage_input.pmsms"


class SageRawOutdir(NodeType):
    filename = "sage_raw_outdir"


class SageRunInfo(NodeType):
    filename = "run_info.json"


class SageResultsJson(NodeType):
    filename = "results.json"


class SagePrecursorsParquet(NodeType):
    filename = "results.sage.parquet"


class SageMatchedFragmentsParquet(NodeType):
    filename = "matched_fragments.sage.parquet"


class SagePin(NodeType):
    filename = "results.sage.pin"


class SagePrecursorsAfterFdr(NodeType):
    filename = "results.sage.fdr.parquet"


class SageSummary(NodeType):
    filename = "results.sage.summary.tsv"


# --- source rules (symlink pre-existing files/dirs, no validation) ---
@R.command("ln -s $(realpath {path}) {tdf}")
def source_bruker_d(path: str):
    return BrukerD[tdf]


@R.command("ln -s $(realpath {path}) {fasta}")
def source_fasta(path: str):
    return Fasta[fasta]


@R.command("ln -s $(realpath {path}) {config}")
def source_unified_config(path: str):
    return UnifiedConfig[config]


# --- compute rules ---
@R.command(
    "venvs/common/bin/d2ms1 {tdf} {ms1}"
    " && test -f {ms1}/tof_row_starts.dat"
    " && test -f {ms1}/tof_urt_diff_index.dat"
    " && test -f {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
)
def tdf2ms1(tdf: BrukerD):
    return Ms1Events[ms1]


@R.command(
    "git/ionmaidenmetal/build/tdf2ms ms2 {tdf} {ms2} --overwrite"
    " && venvs/common/bin/python scripts/tdf2tof2mz.py {tdf} {ms2} {tof2mz}"
)
def tdf2ms2(tdf: BrukerD):
    return Ms2Events[ms2], Tof2Mz[tof2mz]


@R.command(with_config(
    "scale_estimation", "scales",
    'venvs/common/bin/ms1_estimate_scales {ms1} "$CFG" {scales} --dataset_name {dataset}',
))
def estimate_precursor_scales(ms1: Ms1Events, config: UnifiedConfig, dataset: str):
    return ScaleEstimates[scales]


@R.command(with_config(
    "pipeline", "clusters",
    'venvs/common/bin/ms1_select_candidates {ms1} {scale_estimates} "$CFG" {clusters} --dataset_name {dataset}',
))
def select_precursor_candidates(ms1: Ms1Events, scale_estimates: ScaleEstimates, config: UnifiedConfig, dataset: str):
    return RawPrecursorClusters[clusters]


@R.command(with_config(
    "pipeline", "clusters",
    'venvs/common/bin/ms1_postprocess_candidates {tdf} {ms1} {candidates} {scale_estimates} "$CFG" {clusters}'
    ' --dataset_name {dataset}',
))
def postprocess_precursor_candidates(
    tdf: BrukerD, ms1: Ms1Events, candidates: RawPrecursorClusters,
    scale_estimates: ScaleEstimates, config: UnifiedConfig, dataset: str,
):
    return PostprocessedPrecursorClusters[clusters]


@R.command(with_config(
    "precursor_transmission", "transpec",
    'venvs/common/bin/transmit_precursors {tdf} {clusters} "$CFG" {transpec}'
    ' --output-precursors {precursors} --verbose'
    ' && test -f {transpec}/schema.txt',
))
def transmit_precursors_into_fragment_space(tdf: BrukerD, clusters: PostprocessedPrecursorClusters, config: UnifiedConfig):
    return TransmittedMs1Events[transpec], TransmittedPrecursorClusters[precursors]


@R.command(with_config(
    "precursor_filters.mkpmsms", "filtered",
    'venvs/common/bin/filter_mmappet {precursors} "$CFG" {filtered} --verbose',
))
def filter_first_precursors(precursors: TransmittedPrecursorClusters, config: UnifiedConfig):
    return FirstFilterPrecursors[filtered]


@R.command(
    ".venv/bin/necromerge2-stage-dir {staged}"
    " --link events.ms2={ms2}"
    " --link transmitted_precursors.mmappet={transprec}"
    " --link precursor_filter.mmappet={filter_mm}"
)
def stage_mkpmsms_input(ms2: Ms2Events, transprec: TransmittedMs1Events, filter_mm: FirstFilterPrecursors):
    return MkpmsmsInputStaged[staged]


@R.command(with_config(
    "pipeline", "pmsms",
    '.venv/bin/necromerge2-run-mkpmsms {staged} "$CFG" {pmsms} {precursors} --threads {threads}',
), threads=1)
def run_mkpmsms(staged: MkpmsmsInputStaged, config: UnifiedConfig, threads: int):
    return Pmsms[pmsms], Ms2IndexedPrecursors[precursors]


@R.command(with_config(
    "precursor_filters.pre_sage", "filtered",
    'venvs/common/bin/filter_mmappet {precursors} "$CFG" {filtered} --verbose',
))
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, config: UnifiedConfig):
    return PreSageFilteredPrecursors[filtered]


@R.command(with_config(
    "precursor_neighbors", "grid",
    'venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid} --config "$CFG"',
))
def build_precursor_grid_index(precursors: PreSageFilteredPrecursors, tdf: BrukerD, config: UnifiedConfig):
    return PrecursorGridIndex[grid]


@R.command(with_config(
    "precursor_neighbors", "csr",
    "git/ionmaidenmetal/build/precursor_neighbors_csr"
    ' --boxes-input {grid_index}/boxes.mmappet --index-input {grid_index} --output {csr}'
    ' $(venvs/common/bin/precursor-neighbors-params "$CFG" {tdf}) --n-threads {threads}',
), threads=1)
def compute_precursor_neighbors(grid_index: PrecursorGridIndex, tdf: BrukerD, config: UnifiedConfig, threads: int):
    return PrecursorNeighborsCsr[csr]


@R.command(
    "git/ionmaidenmetal/build/tof_filter --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
    " --out-path {score} --n-threads {threads}",
    threads=1,
)
def tof_score_filter(pmsms: Pmsms, neighbors_csr: PrecursorNeighborsCsr, threads: int):
    return NeighborScore[score]   # no config arg — confirmed vestigial in the Snakemake rule


@R.command(with_config(
    "tof_score_filter", "pmsms_out",
    "venvs/common/bin/python -m timstofu.cli.score_based_pmsms_filter"
    ' {pmsms} {precursors} {neighbor_score} "$CFG" {pmsms_out} {precursors_out} --threads {threads}',
), threads=1)
def materialize_tof_filtered_pmsms(pmsms: Pmsms, precursors: PreSageFilteredPrecursors, neighbor_score: NeighborScore, config: UnifiedConfig, threads: int):
    return TofFilteredPmsms[pmsms_out], TofFilteredPrecursors[precursors_out]


@R.command(".venv/bin/necromerge2-stage-sage-input {pmsms} {tof2mz} {precursors} {staged}")
def make_tof_filtered_sage_input(pmsms: TofFilteredPmsms, tof2mz: Tof2Mz, precursors: TofFilteredPrecursors):
    return SageInputStaged[staged]


@R.command(with_config(
    "sage", "outdir",
    '.venv/bin/necromerge2-run-sage {spectra} {fasta} "$CFG" {outdir} {run_info}',
    ext="json",
))
def run_sage(spectra: SageInputStaged, fasta: Fasta, config: UnifiedConfig):
    return SageRawOutdir[outdir], SageRunInfo[run_info]


@R.command(
    "cp {sage_dir}/results.json {results_json}"
    " && cp {sage_dir}/results.sage.pin {pin}"
    " && venvs/common/bin/tsv2parquet {sage_dir}/results.sage.tsv {precursors}"
    " && venvs/common/bin/df_head {precursors}"
    " && venvs/common/bin/tsv2parquet {sage_dir}/matched_fragments.sage.tsv {fragments}"
    " && venvs/common/bin/df_head {fragments}"
)
def extract_sage_results(sage_dir: SageRawOutdir):
    return (
        SageResultsJson[results_json],
        SagePrecursorsParquet[precursors],
        SageMatchedFragmentsParquet[fragments],
        SagePin[pin],
    )


@R.command("venvs/common/bin/sage-filter {precursors} {filtered} --fdr 0.01 && venvs/common/bin/df_head {filtered}")
def sage_filter(precursors: SagePrecursorsParquet):
    return SagePrecursorsAfterFdr[filtered]


@R.command("venvs/common/bin/sage-summarize {filtered} {summary}")
def sage_summarize(filtered: SagePrecursorsAfterFdr):
    return SageSummary[summary]


def sage_pipeline(cfg: dict) -> Pipeline:
    """tof-filtered Sage search chain reproducing the old short_test Snakemake target."""
    P = Pipeline()
    threads = int(cfg.get("threads", 1))
    dataset = str(cfg["dataset"])

    P.tdf = R.source_bruker_d(path=str(cfg["tdf_path"]))
    P.fasta = R.source_fasta(path=str(cfg["fasta_path"]))
    P.config = R.source_unified_config(path=str(cfg["config_path"]))

    P.ms1_events = R.tdf2ms1(P.tdf)
    P.ms2_events, P.tof2mz = R.tdf2ms2(P.tdf)

    P.scale_estimates = R.estimate_precursor_scales(P.ms1_events, P.config, dataset=dataset)
    P.raw_precursor_clusters = R.select_precursor_candidates(P.ms1_events, P.scale_estimates, P.config, dataset=dataset)
    P.postprocessed_precursor_clusters = R.postprocess_precursor_candidates(
        P.tdf, P.ms1_events, P.raw_precursor_clusters, P.scale_estimates, P.config, dataset=dataset
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = R.transmit_precursors_into_fragment_space(
        P.tdf, P.postprocessed_precursor_clusters, P.config
    )
    P.first_filter_precursors = R.filter_first_precursors(P.transmitted_precursor_clusters, P.config)
    P.mkpmsms_input = R.stage_mkpmsms_input(P.ms2_events, P.transmitted_ms1events, P.first_filter_precursors)
    P.pmsms, P.ms2indexed_precursors = R.run_mkpmsms(P.mkpmsms_input, P.config, threads=threads)
    P.pre_sage_filtered_precursors = R.filter_pre_sage_precursors(P.ms2indexed_precursors, P.config)
    P.precursor_grid_index = R.build_precursor_grid_index(P.pre_sage_filtered_precursors, P.tdf, P.config)
    P.precursor_neighbors_csr = R.compute_precursor_neighbors(P.precursor_grid_index, P.tdf, P.config, threads=threads)
    P.neighbor_score = R.tof_score_filter(P.pmsms, P.precursor_neighbors_csr, threads=threads)
    P.tof_filtered_pmsms, P.tof_filtered_precursors = R.materialize_tof_filtered_pmsms(
        P.pmsms, P.pre_sage_filtered_precursors, P.neighbor_score, P.config, threads=threads
    )
    P.sage_input = R.make_tof_filtered_sage_input(P.tof_filtered_pmsms, P.tof2mz, P.tof_filtered_precursors)
    P.sage_raw_outdir, P.sage_run_info = R.run_sage(P.sage_input, P.fasta, P.config)
    P.results_json, P.sage_precursors, P.matched_fragments, P.sage_pin = R.extract_sage_results(P.sage_raw_outdir)
    P.sage_precursors_after_fdr = R.sage_filter(P.sage_precursors)
    P.sage_summary = R.sage_summarize(P.sage_precursors_after_fdr)
    return P
