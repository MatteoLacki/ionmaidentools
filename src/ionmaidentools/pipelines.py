"""Necroflow ionmaiden pipeline.

TODO(regression-db): the old Snakemake `sage_summarize`/`short_test` rules recorded
results into a SQLite regression DB and did an interactive baseline comparison. Neither
is ported here -- `sage_summary` is the pipeline's terminal output.

Fingerprint v3 used.
"""

from __future__ import annotations

import json
import os
import shlex
import tomlkit

from dictodot import DotDict
from necroflow import (
    CommandArgs,
    NodeType,
    Pipeline,
    command,
    output,
    symlink_file,
    text_file,
)

CORES = os.cpu_count() or 1


# --- source node types ---
class BrukerD(NodeType):
    filename = "input.d"


class Fasta(NodeType):
    filename = "fasta.fasta"


class MkpmsmsBinary(NodeType):
    filename = "mkpmsms"


class SageBinary(NodeType):
    filename = "sage"


class SageSummarizeModule(NodeType):
    filename = "sage.py"


class MmappetDataset(NodeType):
    """Base type for outputs that are mmappet directories (see CLAUDE.md's
    "Precursor table format" convention). No `filename` of its own -- every
    concrete mmappet output subclasses this and sets its own."""


# --- compute artifact node types ---
class Ms1Events(NodeType):
    filename = "events.ms1"


class Ms2Events(NodeType):
    filename = "events.ms2"


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


class PrecursorNeighborCountConfig(NodeType):
    filename = "precursor_neighbor_count_config.toml"


class PrecursorCandidateSelectionConfig(NodeType):
    filename = "precursor_candidate_selection_config.toml"


class RawCandidateFeatures(MmappetDataset):
    filename = "raw_candidate_features.mmappet"


class RawPrecursorClusters(MmappetDataset):
    filename = "raw_precursor_clusters.mmappet"


class PrecursorAnnotationConfig(NodeType):
    filename = "precursor_annotation_config.toml"


class PostprocessingConfig(NodeType):
    filename = "postprocessing_config.toml"


class AnnotatedPrecursorClusters(MmappetDataset):
    filename = "annotated_precursor_clusters.mmappet"


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


class TofFilteredPmsms(Pmsms):
    filename = "tof_filtered_pmsms.mmappet"


class TofFilteredPrecursors(PreSageFilteredPrecursors):
    filename = "tof_filtered_precursors.mmappet"


class MzPmsms(Pmsms):
    """A pmsms dataset with a materialized `mz` column (uncorrected), produced by
    timstofu's materialize_pmsms_mz. Accepted wherever Pmsms is."""

    filename = "mz_pmsms.mmappet"


class MzRecalibration(NodeType):
    """Grid-sampled fragment m/z ppm-correction artifact (dimension "mz"),
    produced by recalibrate_pmsms_mz alongside the pmsms it already applied the
    correction to -- kept for inspection/reuse, not re-consumed downstream."""

    filename = "mz_recalibration.mzcalib"


class RecalibratedPmsms(MzPmsms):
    """MzPmsms with its `mz` column corrected in place by recalibrate_pmsms_mz.
    Accepted wherever MzPmsms (or, transitively, Pmsms) is."""

    filename = "recalibrated_pmsms.mmappet"


class SageConfig(NodeType):
    filename = "sage_config.json"


class SageResultsJson(NodeType):
    filename = "results.json"


class Pin(NodeType):
    """A Percolator-IN TSV: SpecId, Label, ScanNr, <features...>, Peptide,
    Proteins. Format is generic (nothing Sage-specific) -- `SageResultsPin`
    and `SagepyRescorePin` both satisfy it, so `mokapot` accepts either."""


class SageResultsPin(Pin):
    filename = "results.sage.pin"


class SageResultsTsv(NodeType):
    filename = "results.sage.tsv"


class SageMatchedFragments(NodeType):
    filename = "matched_fragments.sage.tsv"


class DumpPeptidesBinary(NodeType):
    filename = "dump_peptides"


class DumpPeptidesConfig(NodeType):
    """Just the `database` subdictionary of a `SageConfig` — the only
    settings `dump_peptides` depends on (fasta, enzyme, mods, mass bounds,
    decoy_tag). Sliced out of the full sage config by `extract_dump_peptides_config`
    so `dump_peptides` itself never needs to read (or be valid against)
    unrelated search settings like `precursor_tol`."""

    filename = "dump_peptides_config.json"


class DumpedPeptides(NodeType):
    filename = "peptides.parquet"


class SagepyRescoreConfig(NodeType):
    filename = "sagepy_rescore_config.toml"


class SagepyRescorePredictions(NodeType):
    filename = "psms_with_predictions.parquet"


class SagepyRescorePin(Pin):
    filename = "sagepy_rescore.pin"


class MokapotUsedPin(NodeType):
    filename = "used.pin"


class MokapotPeptides(NodeType):
    filename = "mokapot.peptides.txt"


class MokapotPsms(NodeType):
    filename = "mokapot.psms.txt"


class ConfidentPsmsParquet(NodeType):
    filename = "confident_psms.parquet"


class SagePmsmsMapping(NodeType):
    filename = "sage_mapped_to_pmsms"


class ScoreComparisonPlots(NodeType):
    filename = "score_comparison"


class Mzml(NodeType):
    filename = "mzml.mzML"


class MzmlIdmap(MmappetDataset):
    filename = "idmap.mmappet"


class Mgf(NodeType):
    filename = "spectra.mgf"


class TofFilteredMzml(Mzml):
    pass


class TofFilteredMzmlIdmap(MzmlIdmap):
    pass


class TofFilteredMgf(Mgf):
    filename = "tof_filtered.mgf"


class TarGz(NodeType):
    filename = "archive.tar.gz"


class Png(NodeType):
    filename = "plot.png"


class FragpipeWorkflow(NodeType):
    """A hand-maintained FragPipe `.workflow` file -- treated as a data file
    like Fasta/BrukerD, symlinked in as-is. Never generated or patched by
    the pipeline (unlike Sage's JSON config)."""

    filename = "fragpipe_workflow.workflow"


class FragpipeManifest(NodeType):
    filename = "fragpipe_manifest.tsv"


class FragpipeResultsDir(NodeType):
    filename = "fragpipe_results"


class FragpipeLog(NodeType):
    filename = "fragpipe_full_log.txt"


class FragpipeSummary(NodeType):
    filename = "fragpipe_summary.txt"


class FragpipeDecoyFasta(NodeType):
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
    baked in directly (unlike TofFilteredPmsms, which needs
    materialize_pmsms_mz run on it first to gain an mz column)."""

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


class SageSummary(NodeType):
    filename = "results.sage.summary.tsv"


class RecalibrationPrecursorSelectionConfig(NodeType):
    filename = "recalibration_precursor_selection_config.toml"


class RecalibrationPrecursors(TofFilteredPrecursors):
    """A subset of TofFilteredPrecursors -- accepted wherever it is, e.g. run_sage's
    `precursors` input for the calibration-only pass."""

    filename = "recalibration_precursors.mmappet"


class RecalibrationConfig(NodeType):
    filename = "recalibration_config.toml"


class FragmentTolerance(NodeType):
    """Fragment-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    recalibrate_pmsms_mz from its own fitted model's residuals."""

    filename = "fragment_tolerance.json"


class PrecursorTolerance(NodeType):
    """Precursor-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    recalibrate_precursors from its own fitted model's residuals."""

    filename = "precursor_tolerance.json"


class SerializedMzModel(NodeType):
    """A fitted searchops.models.MzCorrectionModel, dumped via its own save()
    (inspection only -- recalibrate_pmsms_mz/recalibrate_precursors already
    apply the model themselves, nothing downstream reads this back in)."""

    filename = "mz_model"


class RecalibratedPrecursors(TofFilteredPrecursors):
    """search_precursors with `mz` corrected by recalibrate_precursors's own
    fitted model -- accepted wherever TofFilteredPrecursors (or, transitively,
    PreSageFilteredPrecursors) is, e.g. run_sage's final-pass `precursors` input."""

    filename = "recalibrated_precursors.mmappet"


class RecalibratedPpmPlot(Png):
    filename = "recalibrated_ppm.png"


class Predictions(NodeType):
    """sequence,charge -> rt,iim parquet from `git/featureprediction`'s
    `predict_rt_iim` -- the `--predicted-properties` cache SAGE1 needs to
    look up a candidate's predicted RT/IIM during search (it can't derive
    these itself). See plans/better_sage_filtering.md's B.4/B.5."""

    filename = "predictions.parquet"


class RtTolerance(NodeType):
    """rt_tol_sec tolerance JSON (`ValueTolSpline`-shaped, wrapped as
    `{"rt_tol_sec": {...}}` -- see `feature_prediction.tolerance
    .write_value_tol_spline_json`). Produced by both `predict_rt_iim`
    (pre-correction) and `correct_precursors_rt_iim` (post-correction,
    tighter, the one actually fed to SAGE1 when RT/IIM correction runs)."""

    filename = "rt_tolerance.json"


class MobilityTolerance(NodeType):
    """Same as `RtTolerance`, for `mobility_tol`/IIM."""

    filename = "mobility_tolerance.json"


class RtIimFitPlot(Png):
    filename = "rt_iim_fit.png"


class RtIimCorrectedPrecursors(RecalibratedPrecursors):
    """RecalibratedPrecursors (mz-corrected) with `rt`/`inv_ion_mobility`
    also corrected by `feature_prediction.precursor_correction
    .correct_precursors_rt_iim` -- accepted wherever RecalibratedPrecursors
    (or, transitively, TofFilteredPrecursors/PreSageFilteredPrecursors) is,
    e.g. run_sage's final-pass `precursors` input in mode 3. See
    plans/better_sage_filtering.md's B.6."""

    filename = "rt_iim_corrected_precursors.mmappet"


class PrecursorCorrectionRtModel(NodeType):
    filename = "precursor_correction_rt_model.json"


class PrecursorCorrectionIimModels(NodeType):
    """A directory of per-charge XGBoost models (`{charge}.json`, plus a
    pooled fallback `0.json`) -- see `precursor_correction.fit_iim_correction`."""

    filename = "precursor_correction_iim_models"


class PrecursorCorrectionFitPlot(Png):
    filename = "precursor_correction_fit_plot.png"


# --- source rules (symlink pre-existing files/dirs, no validation) ---
@command("ln -s $(realpath {path}) {tdf}")
def source_bruker_d(path: str):
    tdf = output(BrukerD)
    return tdf


@command("ln -s $(realpath {path}) {fasta}")
def source_fasta(path: str):
    fasta = output(Fasta)
    return fasta


@command("ln -s $(realpath {path}) {mkpmsms}")
def source_mkpmsms_binary(path: str):
    mkpmsms = output(MkpmsmsBinary)
    return mkpmsms


@symlink_file
def source_sage_binary(path: str):
    sage_binary = output(SageBinary)
    return sage_binary


@symlink_file
def source_dump_peptides_binary(path: str):
    dump_peptides_binary = output(DumpPeptidesBinary)
    return dump_peptides_binary


# Symlinked so edits to the installed-editable module invalidate sage_summarize.
@symlink_file
def source_sage_summarize_module(path: str):
    sage_summarize_module = output(SageSummarizeModule)
    return sage_summarize_module


# --- compute rules ---
@command(
    "venvs/common/bin/d2ms1 {tdf} {ms1}"
    " && test -f {ms1}/tof_row_starts.dat"
    " && test -f {ms1}/tof_urt_diff_index.dat"
    " && test -f {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
)
def tdf2ms1(tdf: BrukerD):
    ms1 = output(Ms1Events)
    return ms1


@command("git/ionmaidenmetal/build/tdf2ms ms2 {tdf} {ms2} --overwrite")
def tdf2ms2(tdf: BrukerD):
    ms2 = output(Ms2Events)
    return ms2


@text_file
def write_scale_estimation_config(text: str):
    config = output(ScaleEstimationConfig)
    return config


@command("venvs/common/bin/ms1_find_argmaxes {ms1} {config} {argmaxes} {stats}")
def find_ms1_argmaxes(ms1: Ms1Events, config: ScaleEstimationConfig):
    argmaxes = output(ArgmaxSample)
    stats = output(ArgmaxSieveStats)
    return argmaxes, stats


@command(
    "venvs/common/bin/ms1_extract_sample_tensors {ms1} {argmaxes} {config} {tensors}"
)
def extract_ms1_sample_tensors(
    ms1: Ms1Events, argmaxes: ArgmaxSample, config: ScaleEstimationConfig
):
    tensors = output(SampleTensors)
    return tensors


@command(
    "venvs/common/bin/ms1_fit_scale_estimates {argmaxes} {stats} {tensors} {config} {scales}"
)
def fit_ms1_scale_estimates(
    argmaxes: ArgmaxSample,
    stats: ArgmaxSieveStats,
    tensors: SampleTensors,
    config: ScaleEstimationConfig,
):
    scales = output(ScaleEstimates)
    return scales


@text_file
def write_precursor_neighbor_count_config(text: str):
    config = output(PrecursorNeighborCountConfig)
    return config


@text_file
def write_precursor_candidate_selection_config(text: str):
    config = output(PrecursorCandidateSelectionConfig)
    return config


@command(
    "venvs/common/bin/ms1_count_candidate_neighbors {ms1} {scale_estimates} {config} {features}"
)
def count_candidate_neighbors(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    config: PrecursorNeighborCountConfig,
):
    features = output(RawCandidateFeatures)
    return features


@command(
    "venvs/common/bin/ms1_score_candidates {ms1} {scale_estimates} {features} {config} {clusters}"
)
def score_candidates(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    features: RawCandidateFeatures,
    config: PrecursorCandidateSelectionConfig,
):
    clusters = output(RawPrecursorClusters)
    return clusters


@text_file
def write_precursor_annotation_config(text: str):
    config = output(PrecursorAnnotationConfig)
    return config


@text_file
def write_postprocessing_config(text: str):
    config = output(PostprocessingConfig)
    return config


@command(
    "venvs/common/bin/ms1_annotate_candidates {tdf} {ms1} {candidates} {scale_estimates} {config} {annotated}"
)
def annotate_precursor_clusters(
    tdf: BrukerD,
    ms1: Ms1Events,
    candidates: RawPrecursorClusters,
    scale_estimates: ScaleEstimates,
    config: PrecursorAnnotationConfig,
):
    annotated = output(AnnotatedPrecursorClusters)
    return annotated


@command(
    "venvs/common/bin/ms1_decharge_candidates {tdf} {ms1} {annotated} {config} {clusters}"
)
def decharge_precursor_clusters(
    tdf: BrukerD,
    ms1: Ms1Events,
    annotated: AnnotatedPrecursorClusters,
    config: PostprocessingConfig,
):
    clusters = output(PostprocessedPrecursorClusters)
    return clusters


@text_file
def write_precursor_transmission_config(text: str):
    config = output(PrecursorTransmissionConfig)
    return config


@command(
    "venvs/common/bin/transmit_precursors {tdf} {clusters} {config} {transpec}"
    " --output-precursors {precursors} --verbose"
    " && test -f {transpec}/schema.txt"
)
def transmit_precursors_into_fragment_space(
    tdf: BrukerD,
    clusters: PostprocessedPrecursorClusters,
    config: PrecursorTransmissionConfig,
):
    transpec = output(TransmittedMs1Events)
    precursors = output(TransmittedPrecursorClusters)
    return transpec, precursors


@command(
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
)
def filter_first_precursors(precursors: TransmittedPrecursorClusters, filter: str):
    filtered = output(FirstFilterPrecursors)
    return filtered


@text_file
def write_pseudomsms_config(text: str):
    config = output(PseudomsmsConfig)
    return config


@command(
    "{mkpmsms} --fragments {ms2} --transmitted-precursors {transprec}"
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
    mkpmsms: MkpmsmsBinary,
    ms2: Ms2Events,
    transprec: TransmittedMs1Events,
    filter_mm: FirstFilterPrecursors,
    config: PseudomsmsConfig,
):
    pmsms = output(Pmsms)
    return pmsms


@command(
    "venvs/common/bin/cut_and_index_precursors {filter_mm} {pmsms}/dataindex.mmappet {precursors}"
)
def cut_and_index_precursors(filter_mm: FirstFilterPrecursors, pmsms: Pmsms):
    precursors = output(Ms2IndexedPrecursors)
    return precursors


@command(
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
)
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, filter: str):
    filtered = output(PreSageFilteredPrecursors)
    return filtered


@text_file
def write_precursor_neighbors_config(text: str):
    config = output(PrecursorNeighborsConfig)
    return config


@command(
    "venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid} --config {config}"
)
def build_precursor_grid_index(
    precursors: PreSageFilteredPrecursors,
    tdf: BrukerD,
    config: PrecursorNeighborsConfig,
):
    grid = output(PrecursorGridIndex)
    return grid


@command(
    "git/ionmaidenmetal/build/precursor_neighbors_csr"
    " --boxes-input {grid_index}/boxes.mmappet --index-input {grid_index} --output {csr}"
    " $(venvs/common/bin/precursor-neighbors-params {config} {tdf})"
    " --n-threads {threads}",
    threads=CORES,
)
def compute_precursor_neighbors(
    grid_index: PrecursorGridIndex,
    tdf: BrukerD,
    config: PrecursorNeighborsConfig,
):
    csr = output(PrecursorNeighborsCsr)
    return csr


@command(
    "git/ionmaidenmetal/build/tof_filter --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
    " --out-path {score} --n-threads {threads}",
    threads=CORES,
)
def tof_score_filter(pmsms: Pmsms, neighbors_csr: PrecursorNeighborsCsr):
    score = output(
        NeighborScore
    )  # no config arg -- confirmed vestigial in the Snakemake rule
    return score


@command(
    "venvs/common/bin/python -m timstofu.cli.score_based_pmsms_filter"
    " {pmsms} {precursors} {neighbor_score} {pmsms_out} {precursors_out}"
    " --threads {threads} --score-margin {score_margin}",
    threads=CORES,
)
def materialize_tof_filtered_pmsms(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    neighbor_score: NeighborScore,
    score_margin: int | float,
):
    pmsms_out = output(TofFilteredPmsms)
    precursors_out = output(TofFilteredPrecursors)
    return pmsms_out, precursors_out


@text_file
def write_recalibration_precursor_selection_config(text: str):
    config = output(RecalibrationPrecursorSelectionConfig)
    return config


@command(
    "venvs/common/bin/python -m timstofu.cli.select_recalibration_precursors"
    " {precursors} {config} {selected}"
)
def select_recalibration_precursors(
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationPrecursorSelectionConfig,
):
    selected = output(RecalibrationPrecursors)
    return selected


@text_file
def write_recalibration_config(text: str):
    config = output(RecalibrationConfig)
    return config


@command(
    "venvs/common/bin/materialize_pmsms_mz {input_pmsms} {tdf} {output_pmsms}"
)
def materialize_pmsms_mz(input_pmsms: Pmsms, tdf: BrukerD):
    output_pmsms = output(MzPmsms)
    return output_pmsms


@command(
    "venvs/common/bin/recalibrate-pmsms-mz {sage_results_tsv} {matched_fragments} {mz_pmsms}"
    " {output_pmsms} {mz_recalibration} {tolerance} {plot} --config {config} --fdr {fdr}"
)
def recalibrate_pmsms_mz(
    sage_results_tsv: SageResultsTsv,
    matched_fragments: SageMatchedFragments,
    mz_pmsms: MzPmsms,
    config: RecalibrationConfig,
    fdr: int | float,
):
    output_pmsms = output(RecalibratedPmsms)
    mz_recalibration = output(MzRecalibration)
    tolerance = output(FragmentTolerance)
    plot = output(Png)
    return output_pmsms, mz_recalibration, tolerance, plot


@command(
    "venvs/common/bin/recalibrate-precursors {sage_results_tsv} {precursors}"
    " {output_precursors} {tolerance} {plot} {model} --config {config} --fdr {fdr}"
)
def recalibrate_precursors(
    sage_results_tsv: SageResultsTsv,
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationConfig,
    fdr: int | float,
):
    output_precursors = output(RecalibratedPrecursors)
    tolerance = output(PrecursorTolerance)
    plot = output(Png)
    model = output(SerializedMzModel)
    return output_precursors, tolerance, plot, model


@command(
    "venvs/common/bin/plot-recalibrated-ppm {initial_sage_results_tsv} {sage_results_tsv}"
    " {initial_matched_fragments} {matched_fragments} {precursor_tolerance} {fragment_tolerance}"
    " {plot} --fdr {fdr}"
)
def plot_recalibrated_ppm(
    initial_sage_results_tsv: SageResultsTsv,
    sage_results_tsv: SageResultsTsv,
    initial_matched_fragments: SageMatchedFragments,
    matched_fragments: SageMatchedFragments,
    precursor_tolerance: PrecursorTolerance,
    fragment_tolerance: FragmentTolerance,
    fdr: int | float,
):
    plot = output(RecalibratedPpmPlot)
    return plot


@command(
    ".venv/bin/python -m necroflow.tools.config_set"
    " {sage_config} {workdir}/precursor_tol_updated.json"
    " --target precursor_tol.ppm --source {precursor_tolerance} --source-field ppm"
    " && .venv/bin/python -m necroflow.tools.config_set"
    " {workdir}/precursor_tol_updated.json {recalibrated_sage_config}"
    " --target fragment_tol.ppm --source {fragment_tolerance} --source-field ppm"
)
def update_sage_config(
    sage_config: SageConfig,
    precursor_tolerance: PrecursorTolerance,
    fragment_tolerance: FragmentTolerance,
):
    recalibrated_sage_config = output(SageConfig)
    return recalibrated_sage_config


@text_file
def write_sage_config(text: str):
    config = output(SageConfig)
    return config


@text_file
def write_dump_peptides_config(text: str):
    dump_peptides_config = output(DumpPeptidesConfig)
    return dump_peptides_config


@command(
    "{dump_peptides_binary} -f {fasta} -c {dump_peptides_config} -o {peptides}"
)
def dump_peptides(
    fasta: Fasta,
    dump_peptides_config: DumpPeptidesConfig,
    dump_peptides_binary: DumpPeptidesBinary,
):
    peptides = output(DumpedPeptides)
    return peptides


@command(
    "{sage_binary} --version"
    " && {sage_binary} -f {fasta} --annotate-matches --write-pin"
    " --output_directory {workdir} --pmsms {pmsms} --precursors {precursors}"
    " {sage_config}"
    " && test -f {results_json} && test -f {results_pin}"
    " && test -f {results_tsv} && test -f {matched_fragments}",
    threads=CORES,
)
def run_sage(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
    sage_binary: SageBinary,
):
    results_json = output(SageResultsJson)
    results_pin = output(SageResultsPin)
    results_tsv = output(SageResultsTsv)
    matched_fragments = output(SageMatchedFragments)
    return results_json, results_pin, results_tsv, matched_fragments


@command(
    "{sage_binary} --version"
    " && {sage_binary} -f {fasta} --annotate-matches --write-pin"
    " --output_directory {workdir} --pmsms {pmsms} --precursors {precursors}"
    " --predicted-properties {predictions}"
    " {sage_config}"
    " && test -f {results_json} && test -f {results_pin}"
    " && test -f {results_tsv} && test -f {matched_fragments}",
    threads=CORES,
)
def run_sage_with_predicted_properties(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
    sage_binary: SageBinary,
    predictions: Predictions,
):
    """Same as `run_sage`, plus `--predicted-properties` -- a distinct rule
    rather than a conditional flag on `run_sage` itself (necroflow's
    `@command` templates are static; see `_mokapot_command` below for the
    Python-callback alternative, not used here since `predictions` would
    need to become an optional DAG input either way). Only used for mode 3
    (mz+RT+IIM recalibration)'s final pass -- `sage_config` must already
    have `rt_tol_sec`/`mobility_tol` set (`update_sage_config_rt_iim`),
    same requirement `Input::build()` enforces Rust-side. See
    plans/better_sage_filtering.md's B.4-B.6.
    """
    results_json = output(SageResultsJson)
    results_pin = output(SageResultsPin)
    results_tsv = output(SageResultsTsv)
    matched_fragments = output(SageMatchedFragments)
    return results_json, results_pin, results_tsv, matched_fragments


@command(
    "venvs/featureprediction/bin/feature-prediction-generate"
    " {dumped_peptides} {sage_results_tsv} {predictions} {rt_tolerance} {mobility_tolerance} {plot}"
    " --min-charge {min_charge} --max-charge {max_charge}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi} --fdr {fdr}"
)
def predict_rt_iim(
    dumped_peptides: DumpedPeptides,
    sage_results_tsv: SageResultsTsv,
    min_charge: int,
    max_charge: int,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    fdr: int | float,
):
    predictions = output(Predictions)
    rt_tolerance = output(RtTolerance)
    mobility_tolerance = output(MobilityTolerance)
    plot = output(RtIimFitPlot)
    return predictions, rt_tolerance, mobility_tolerance, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-correct-precursors"
    " {sage_results_tsv} {predictions} {mz_corrected_precursors}"
    " {output_precursors} {rt_tolerance} {mobility_tolerance}"
    " {rt_model} {iim_models} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi} --fdr {fdr}"
)
def correct_precursors_rt_iim(
    sage_results_tsv: SageResultsTsv,
    predictions: Predictions,
    mz_corrected_precursors: RecalibratedPrecursors,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    fdr: int | float,
):
    output_precursors = output(RtIimCorrectedPrecursors)
    rt_tolerance = output(RtTolerance)
    mobility_tolerance = output(MobilityTolerance)
    rt_model = output(PrecursorCorrectionRtModel)
    iim_models = output(PrecursorCorrectionIimModels)
    plot = output(PrecursorCorrectionFitPlot)
    return output_precursors, rt_tolerance, mobility_tolerance, rt_model, iim_models, plot


@command(
    ".venv/bin/python -m necroflow.tools.config_set"
    " {sage_config} {workdir}/rt_tol_updated.json"
    " --target rt_tol_sec --source {rt_tolerance} --source-field rt_tol_sec"
    " && .venv/bin/python -m necroflow.tools.config_set"
    " {workdir}/rt_tol_updated.json {recalibrated_sage_config}"
    " --target mobility_tol --source {mobility_tolerance} --source-field mobility_tol"
)
def update_sage_config_rt_iim(
    sage_config: SageConfig,
    rt_tolerance: RtTolerance,
    mobility_tolerance: MobilityTolerance,
):
    """Chains onto `update_sage_config`'s output (mz's `precursor_tol`/
    `fragment_tol` already patched) -- `sage_config` here is that rule's
    `recalibrated_sage_config`, not the plain `write_sage_config` output.
    `predicted_properties` itself is *not* set here: like `--pmsms`/
    `--precursors`/`-f {fasta}`, it's a path-valued override passed
    directly as a CLI flag (`run_sage_with_predicted_properties`), not
    embedded into the config JSON via `config_set`."""
    recalibrated_sage_config = output(SageConfig)
    return recalibrated_sage_config


def _mokapot_command(args: CommandArgs) -> str:
    """Python command callback, not a static template -- lets `--plugin`
    be added conditionally (empty for the plain-Sage-PIN call, `--plugin
    xgboost` for the sagepy-rescore call) without necroflow's string-
    template placeholders needing to express a conditional substring.
    """
    pin = shlex.quote(str(args.inputs.pin))
    used_pin = shlex.quote(str(args.outputs.used_pin))
    peptides = shlex.quote(str(args.outputs.peptides))
    psms = shlex.quote(str(args.outputs.psms))
    workdir = shlex.quote(str(args.workdir))
    plugin_flag = f" --plugin {args.config.plugin}" if args.config.plugin else ""
    return (
        f"venvs/mokapot/bin/python scripts/mokapot_pin_adapter.py -i {pin} -o {used_pin}"
        f" && venvs/mokapot/bin/mokapot {used_pin} --dest_dir {workdir}"
        f" --train_fdr {args.config.train_fdr} --test_fdr {args.config.test_fdr}"
        f"{plugin_flag}"
        f" && test -f {peptides} && test -f {psms}"
    )


@command(_mokapot_command)
def mokapot(
    pin: Pin,
    train_fdr: float = 0.05,
    test_fdr: float = 0.01,
    plugin: str | None = None,
):
    used_pin = output(MokapotUsedPin)
    peptides = output(MokapotPeptides)
    psms = output(MokapotPsms)
    return used_pin, peptides, psms


@text_file
def write_sagepy_rescore_config(text: str):
    config = output(SagepyRescoreConfig)
    return config


def _sagepy_rescore_prediction_config(config: dict) -> dict:
    """Remove settings owned by the downstream mokapot rule."""
    prediction_config = dict(config)
    prediction_config.pop("train_fdr", None)
    prediction_config.pop("test_fdr", None)
    return prediction_config


@command(
    "venvs/sagepy_rescore/bin/sagepy-rescore-from-sage"
    " --psms-parquet {sage_results_tsv} --matched-fragments-parquet {sage_matched_fragments}"
    " --output {workdir}"
    " --with-predictors --predict-only"
    " --config {config}"
    " && test -f {predictions}"
)
def run_sagepy_rescore_predict(
    sage_results_tsv: SageResultsTsv,
    sage_matched_fragments: SageMatchedFragments,
    config: SagepyRescoreConfig,
):
    predictions = output(SagepyRescorePredictions)
    return predictions


@command(
    "venvs/sagepy_rescore/bin/python scripts/write_sagepy_rescore_pin.py"
    " -i {predictions} -o {pin}"
)
def write_sagepy_rescore_pin(predictions: SagepyRescorePredictions):
    pin = output(SagepyRescorePin)
    return pin


@command("venvs/common/bin/sage-summarize-raw {sage_results_tsv} {summary} --fdr {fdr}")
def sage_summarize(
    sage_results_tsv: SageResultsTsv,
    sage_summarize_module: SageSummarizeModule,
    fdr: int | float,
):
    summary = output(SageSummary)
    return summary


@command("venvs/common/bin/sage-filter {sage_results_tsv} {confident_psms} --fdr {fdr}")
def filter_sage_results(sage_results_tsv: SageResultsTsv, fdr: int | float):
    confident_psms = output(ConfidentPsmsParquet)
    return confident_psms


@command(
    "venvs/common/bin/sage-pmsms-mapper {confident_psms} {matched_fragments}"
    " {precursors} {pmsms} {mapped}"
)
def sage_map_to_pmsms(
    confident_psms: ConfidentPsmsParquet,
    matched_fragments: SageMatchedFragments,
    precursors: PreSageFilteredPrecursors,
    pmsms: Pmsms,
):
    mapped = output(SagePmsmsMapping)
    return mapped


@command(
    "venvs/common/bin/sage_score_mapper {precursors} {pmsms}"
    " {mapping}/precursors.parquet {mapping}/mapping.parquet"
    " --config {config} -o {plots}"
)
def score_comparison(
    precursors: PreSageFilteredPrecursors,
    pmsms: Pmsms,
    mapping: SagePmsmsMapping,
    config: PseudomsmsConfig,
):
    plots = output(ScoreComparisonPlots)
    return plots


@command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --threads {threads} --numpress --zlib-level 9"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_search_pmsms_to_mzml(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
):
    mzml = output(TofFilteredMzml)
    idmap = output(TofFilteredMzmlIdmap)
    return mzml, idmap


@command(
    "venvs/common/bin/msms2mgf_multicharge {pmsms} {precursors} {config_path} {mgf}"
    " --threads_cnt {threads}"
    " && test -f {mgf}",
    threads=CORES,
)
def convert_search_pmsms_to_mgf(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    config_path: str,
):
    mgf = output(TofFilteredMgf)
    return mgf


@command(
    "tar -cf - -C $(dirname {path}) $(basename {path}) | pigz -p {threads} > {archive}"
    " && test -f {archive}",
    threads=CORES,
)
def compress_with_pigz(path: NodeType):
    archive = output(TarGz)
    return archive


@command("ln -s $(realpath {path}) {workflow}")
def source_fragpipe_workflow(path: str):
    workflow = output(FragpipeWorkflow)
    return workflow


@command(
    # `cd` changes the meaning of every relative path substituted into this
    # template, including the philosopher binary's own repo-relative path and
    # {decoy_fasta} itself -- so all of them are frozen to absolute via shell
    # variables *before* the cd, same reasoning as source_fasta's realpath.
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
def generate_fragpipe_decoy_fasta(fasta: Fasta):
    decoy_fasta = output(FragpipeDecoyFasta)
    return decoy_fasta


@command(
    "sed -e 's|^database\\.db-path=.*|database.db-path={decoy_fasta}|'"
    " {workflow} > {patched_workflow}"
)
def patch_fragpipe_workflow(workflow: FragpipeWorkflow, decoy_fasta: FragpipeDecoyFasta):
    patched_workflow = output(FragpipeWorkflow)
    return patched_workflow


@command(
    "venvs/common/bin/python scripts/simulate_peptides_to_pmsms.py"
    " {fasta} {pmsms} {precursors}"
    " --charges {charges} --max-peptides-per-protein {max_peptides_per_protein}"
    " --seed {seed}"
    " && test -f {pmsms}/schema.txt && test -f {pmsms}/dataindex.mmappet/schema.txt"
    " && test -f {precursors}/schema.txt"
)
def simulate_peptides_to_pmsms(
    fasta: Fasta,
    charges: str,
    max_peptides_per_protein: int,
    seed: int,
):
    pmsms = output(SyntheticPmsms)
    precursors = output(SyntheticPrecursors)
    return pmsms, precursors


@command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --threads {threads} --numpress --zlib-level 9"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_synthetic_pmsms_to_mzml(
    pmsms: SyntheticPmsms, precursors: SyntheticPrecursors
):
    mzml = output(SyntheticMzml)
    idmap = output(SyntheticMzmlIdmap)
    return mzml, idmap


@command(
    "venvs/common/bin/msms2mgf {pmsms} {precursors} {config_path} {mgf}"
    " && test -f {mgf}"
)
def convert_synthetic_pmsms_to_mgf(
    pmsms: SyntheticPmsms,
    precursors: SyntheticPrecursors,
    config_path: str,
):
    mgf = output(SyntheticMgf)
    return mgf


@command('printf "%s\\tA\\t1\\tDDA" "$(realpath {mzml})" > {manifest}')
def write_fragpipe_manifest(mzml: TofFilteredMzml):
    manifest = output(FragpipeManifest)
    return manifest


@command(
    "software/fragpipe/fragpipe-24.0/bin/fragpipe --headless"
    " --workflow {workflow} --manifest {manifest} --workdir {dir}"
    " --config-tools-folder software/fragpipe/fragpipe-24.0/tools"
    " --threads {threads} --ram {ram}"
    ' && test -n "$(ls {dir}/log_*.txt 2>/dev/null)"',
    threads=CORES,
)
def run_fragpipe(manifest: FragpipeManifest, workflow: FragpipeWorkflow, ram: int = 0):
    dir = output(FragpipeResultsDir)
    return dir


@command("cp $(ls {dir}/log_*.txt | head -n 1) {log}")
def extract_fragpipe_log(dir: FragpipeResultsDir):
    log = output(FragpipeLog)
    return log


@command(
    "grep -m 1 temp {log} > {summary}"
    " && grep -n 'MASS CALIBRATION' {log} -A 8 >> {summary}"
    " && grep 'Final report numbers after FDR filtering, and post-processing' {log} >> {summary}"
)
def summarize_fragpipe(log: FragpipeLog):
    summary = output(FragpipeSummary)
    return summary


def ionmaiden_pipeline(P: Pipeline, config: dict) -> None:
    """Main IonMaiden pipeline."""
    cfg = DotDict.Recursive(config)

    # Acquisition
    P.tdf = source_bruker_d(P, path=cfg.tdf_path)
    P.fasta = source_fasta(P, path=cfg.fasta_path)
    P.mkpmsms_binary = source_mkpmsms_binary(
        P, path="git/ionmaidenmetal/build/mkpmsms"
    )
    P.sage_summarize_module = source_sage_summarize_module(
        P, path="git/searchops/src/searchops/sage.py"
    )
    P.sage_binary = source_sage_binary(P, path="software/sage/devel_fixed/sage")
    P.dump_peptides_binary = source_dump_peptides_binary(
        P, path="software/sage/devel_fixed/dump_peptides"
    )

    # Raw Extraction
    P.ms1_events = tdf2ms1(P, P.tdf)
    P.ms2_events = tdf2ms2(P, P.tdf)

    # MS1 Scale Calibration
    P.scale_estimation_config = write_scale_estimation_config(
        P, text=tomlkit.dumps(cfg.scale_estimation)
    )
    P.argmaxes, P.argmax_sieve_stats = find_ms1_argmaxes(
        P, P.ms1_events, P.scale_estimation_config
    )
    P.sample_tensors = extract_ms1_sample_tensors(
        P, P.ms1_events, P.argmaxes, P.scale_estimation_config
    )
    P.scale_estimates = fit_ms1_scale_estimates(
        P,
        P.argmaxes,
        P.argmax_sieve_stats,
        P.sample_tensors,
        P.scale_estimation_config,
    )

    # Precursor Selection
    P.precursor_neighbor_count_config = write_precursor_neighbor_count_config(
        P, text=tomlkit.dumps(cfg.precursor_neighbor_count)
    )
    P.precursor_candidate_selection_config = write_precursor_candidate_selection_config(
        P, text=tomlkit.dumps(cfg.precursor_candidate_selection)
    )
    P.raw_candidate_features = count_candidate_neighbors(
        P,
        P.ms1_events,
        P.scale_estimates,
        P.precursor_neighbor_count_config,
    )
    P.raw_precursor_clusters = score_candidates(
        P,
        P.ms1_events,
        P.scale_estimates,
        P.raw_candidate_features,
        P.precursor_candidate_selection_config,
    )

    # Precursor Postprocessing
    P.precursor_annotation_config = write_precursor_annotation_config(
        P, text=tomlkit.dumps(cfg.precursor_annotation)
    )
    P.postprocessing_config = write_postprocessing_config(
        P, text=tomlkit.dumps(cfg.postprocessing_of_precursors)
    )
    P.annotated_precursor_clusters = annotate_precursor_clusters(
        P,
        P.tdf,
        P.ms1_events,
        P.raw_precursor_clusters,
        P.scale_estimates,
        P.precursor_annotation_config,
    )
    P.postprocessed_precursor_clusters = decharge_precursor_clusters(
        P,
        P.tdf,
        P.ms1_events,
        P.annotated_precursor_clusters,
        P.postprocessing_config,
    )

    # Precursor Transmission
    P.precursor_transmission_config = write_precursor_transmission_config(
        P, text=tomlkit.dumps(cfg.precursor_transmission)
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = (
        transmit_precursors_into_fragment_space(
            P,
            P.tdf,
            P.postprocessed_precursor_clusters,
            P.precursor_transmission_config,
        )
    )

    P.first_filter_precursors = filter_first_precursors(
        P,
        P.transmitted_precursor_clusters,
        filter=cfg.precursor_filters.mkpmsms.get("filter", ""),
    )

    # Pseudo-MS/MS Assembly
    P.pseudomsms_config = write_pseudomsms_config(P, text=tomlkit.dumps(cfg.pseudomsms))
    P.pmsms = run_mkpmsms_binary(
        P,
        P.mkpmsms_binary,
        P.ms2_events,
        P.transmitted_ms1events,
        P.first_filter_precursors,
        P.pseudomsms_config,
    )

    # Precursor Indexing
    P.ms2indexed_precursors = cut_and_index_precursors(
        P, P.first_filter_precursors, P.pmsms
    )

    P.pre_sage_filtered_precursors = filter_pre_sage_precursors(
        P,
        P.ms2indexed_precursors,
        filter=cfg.precursor_filters.pre_sage.get("filter", ""),
    )

    # Neighbor Graph
    if "precursor_neighbors" in cfg:
        P.precursor_neighbors_config = write_precursor_neighbors_config(
            P, text=tomlkit.dumps(cfg.precursor_neighbors)
        )
        P.precursor_grid_index = build_precursor_grid_index(
            P,
            P.pre_sage_filtered_precursors,
            P.tdf,
            P.precursor_neighbors_config,
        )
        P.precursor_neighbors_csr = compute_precursor_neighbors(
            P,
            P.precursor_grid_index,
            P.tdf,
            P.precursor_neighbors_config,
        )
    elif "tof_score_filter" in cfg:
        raise ValueError(
            "[tof_score_filter] is configured but [precursor_neighbors] is not -- "
            "ToF Score Filtering requires a neighbor graph. Add [precursor_neighbors] "
            "or remove [tof_score_filter]."
        )

    # ToF Score Filtering
    if "tof_score_filter" in cfg:
        P.neighbor_score = tof_score_filter(P, P.pmsms, P.precursor_neighbors_csr)

        P.search_pmsms, P.search_precursors = materialize_tof_filtered_pmsms(
            P,
            P.pmsms,
            P.pre_sage_filtered_precursors,
            P.neighbor_score,
            score_margin=cfg.tof_score_filter.score_margin,
        )
    else:
        P.search_pmsms = P.pmsms
        P.search_precursors = P.pre_sage_filtered_precursors

    P.search_mz_pmsms = materialize_pmsms_mz(P, P.search_pmsms, P.tdf)
    final_mz_pmsms = P.search_mz_pmsms
    final_precursors = P.search_precursors

    # Search

    if "sage" in cfg:
        P.sage_config = write_sage_config(
            P, text=json.dumps(cfg.sage, sort_keys=True, indent=2) + "\n"
        )

        # Derived straight from `cfg.sage.database` (not sliced out of the
        # `sage_config` node above) so two job configs that only differ in
        # unrelated sage settings (precursor_tol, report_psms, ...) but
        # share the same digestion settings produce the same
        # content-addressed node here, reusing one `dump_peptides` run.
        P.dump_peptides_config = write_dump_peptides_config(
            P, text=json.dumps(cfg.sage.database, sort_keys=True, indent=2) + "\n"
        )
        P.dumped_peptides = dump_peptides(
            P, P.fasta, P.dump_peptides_config, P.dump_peptides_binary
        )

        if "recalibration" in cfg:
            P.recalibration_precursor_selection_config = (
                write_recalibration_precursor_selection_config(
                    P, text=tomlkit.dumps(cfg.recalibration_precursor_selection)
                )
            )
            P.recalibration_precursors = select_recalibration_precursors(
                P,
                P.search_precursors,
                P.recalibration_precursor_selection_config,
            )
            (
                P.filtered_sage_results_json,
                P.filtered_sage_results_pin,
                P.filtered_sage_results_tsv,
                P.filtered_sage_matched_fragments,
            ) = run_sage(
                P,
                P.search_mz_pmsms,
                P.recalibration_precursors,
                P.fasta,
                P.sage_config,
                P.sage_binary,
            )
            P.recalibration_config = write_recalibration_config(
                P, text=tomlkit.dumps(cfg.recalibration)
            )
            (
                P.recalibrated_mz_pmsms,
                P.fragment_mz_recalibration,
                P.fragment_mz_search_tolerance,
                P.fragment_mz_recalibration_fit_plot,
            ) = recalibrate_pmsms_mz(
                P,
                P.filtered_sage_results_tsv,
                P.filtered_sage_matched_fragments,
                P.search_mz_pmsms,
                P.recalibration_config,
                fdr=cfg.sage_summarize.fdr,
            )
            (
                P.recalibrated_precursors,
                P.precursor_mz_search_tolerance,
                P.precursor_mz_recalibration_fit_plot,
                P.precursro_mz_recalibration_model_serialization,
            ) = recalibrate_precursors(
                P,
                P.filtered_sage_results_tsv,
                P.search_precursors,
                P.recalibration_config,
                fdr=cfg.sage_summarize.fdr,
            )
            P.recalibrated_sage_config = update_sage_config(
                P, P.sage_config, P.precursor_mz_search_tolerance, P.fragment_mz_search_tolerance,
            )

            # Mode 3 (mz + RT + IIM) nests inside mode 2 (mz alone) since it
            # chains onto mz's already-corrected outputs -- gated by its own
            # key so mode 2 stays the default when "rt_iim" is absent, not
            # silently upgraded. See plans/better_sage_filtering.md's B.6.
            if "rt_iim" in cfg.recalibration:
                # min_charge/max_charge: explicit `[recalibration.rt_iim]`
                # override if given, else mirror whatever charge range the
                # SAGE run this feeds will itself search -- `cfg.sage`'s own
                # `precursor_charge` if set, else SAGE's own compiled-in
                # default (2, 4) (`crates/sage-cli/src/input.rs`:
                # `precursor_charge.unwrap_or((2, 4))`). Real job configs
                # (e.g. jobs/f9477_gam_test.toml) often leave
                # `sage.precursor_charge` unset entirely, relying on
                # per-precursor charges instead -- this must not crash on
                # that, and should match SAGE's real default when it does.
                default_precursor_charge = cfg.sage.get("precursor_charge", (2, 4))
                rt_iim_min_charge = cfg.recalibration.rt_iim.get(
                    "min_charge", default_precursor_charge[0]
                )
                rt_iim_max_charge = cfg.recalibration.rt_iim.get(
                    "max_charge", default_precursor_charge[1]
                )
                # tolerance_percentiles: required explicitly in
                # `[recalibration.rt_iim]`, never inherited from
                # `cfg.recalibration`'s own (mz-scoped) value -- each
                # dimension gets its own percentiles, deliberately not
                # shared just because the values happen to look similar.
                rt_iim_tolerance_lo, rt_iim_tolerance_hi = cfg.recalibration.rt_iim[
                    "tolerance_percentiles"
                ]
                (
                    P.predictions,
                    P.rt_tolerance,
                    P.mobility_tolerance,
                    P.rt_iim_fit_plot,
                ) = predict_rt_iim(
                    P,
                    P.dumped_peptides,
                    P.filtered_sage_results_tsv,
                    min_charge=rt_iim_min_charge,
                    max_charge=rt_iim_max_charge,
                    tolerance_lo=rt_iim_tolerance_lo,
                    tolerance_hi=rt_iim_tolerance_hi,
                    fdr=cfg.sage_summarize.fdr,
                )
                (
                    P.rt_iim_corrected_precursors,
                    P.precursor_correction_rt_tolerance,
                    P.precursor_correction_mobility_tolerance,
                    P.precursor_correction_rt_model,
                    P.precursor_correction_iim_models,
                    P.precursor_correction_fit_plot,
                ) = correct_precursors_rt_iim(
                    P,
                    P.filtered_sage_results_tsv,
                    P.predictions,
                    P.recalibrated_precursors,
                    tolerance_lo=rt_iim_tolerance_lo,
                    tolerance_hi=rt_iim_tolerance_hi,
                    fdr=cfg.sage_summarize.fdr,
                )
                P.recalibrated_sage_config_rt_iim = update_sage_config_rt_iim(
                    P,
                    P.recalibrated_sage_config,
                    P.precursor_correction_rt_tolerance,
                    P.precursor_correction_mobility_tolerance,
                )
                (
                    P.sage_results_json,
                    P.sage_results_pin,
                    P.sage_results_tsv,
                    P.sage_matched_fragments,
                ) = run_sage_with_predicted_properties(
                    P,
                    P.recalibrated_mz_pmsms,
                    P.rt_iim_corrected_precursors,
                    P.fasta,
                    P.recalibrated_sage_config_rt_iim,
                    P.sage_binary,
                    P.predictions,
                )
                final_precursors = P.rt_iim_corrected_precursors
            else:
                (
                    P.sage_results_json,
                    P.sage_results_pin,
                    P.sage_results_tsv,
                    P.sage_matched_fragments,
                ) = run_sage(
                    P,
                    P.recalibrated_mz_pmsms,
                    P.recalibrated_precursors,
                    P.fasta,
                    P.recalibrated_sage_config,
                    P.sage_binary,
                )
                final_precursors = P.recalibrated_precursors

            P.recalibrated_ppm_plot = plot_recalibrated_ppm(
                P,
                P.filtered_sage_results_tsv,
                P.sage_results_tsv,
                P.filtered_sage_matched_fragments,
                P.sage_matched_fragments,
                P.precursor_mz_search_tolerance,
                P.fragment_mz_search_tolerance,
                fdr=cfg.sage_summarize.fdr,
            )
            P.confident_psms = filter_sage_results(
                P, P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
            )
            P.sage_pmsms_mapping = sage_map_to_pmsms(
                P,
                P.confident_psms,
                P.sage_matched_fragments,
                P.search_precursors,
                P.recalibrated_mz_pmsms,
            )
            P.score_comparison = score_comparison(
                P,
                P.search_precursors,
                P.search_pmsms,
                P.sage_pmsms_mapping,
                P.pseudomsms_config,
            )
            final_mz_pmsms = P.recalibrated_mz_pmsms
        else:
            (
                P.sage_results_json,
                P.sage_results_pin,
                P.sage_results_tsv,
                P.sage_matched_fragments,
            ) = run_sage(
                P,
                P.search_mz_pmsms,
                P.search_precursors,
                P.fasta,
                P.sage_config,
                P.sage_binary,
            )
            P.confident_psms = filter_sage_results(
                P, P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
            )
            P.sage_pmsms_mapping = sage_map_to_pmsms(
                P,
                P.confident_psms,
                P.sage_matched_fragments,
                P.search_precursors,
                P.search_mz_pmsms,
            )
            P.score_comparison = score_comparison(
                P,
                P.search_precursors,
                P.search_pmsms,
                P.sage_pmsms_mapping,
                P.pseudomsms_config,
            )

        P.mokapot_used_pin, P.mokapot_peptides, P.mokapot_psms = mokapot(
            P, P.sage_results_pin
        )

        if "sagepy_rescore" in cfg:
            prediction_config = _sagepy_rescore_prediction_config(
                cfg.sagepy_rescore
            )
            P.sagepy_rescore_config = write_sagepy_rescore_config(
                P, text=tomlkit.dumps(prediction_config)
            )
            P.sagepy_rescore_predictions = run_sagepy_rescore_predict(
                P,
                P.sage_results_tsv,
                P.sage_matched_fragments,
                P.sagepy_rescore_config,
            )
            P.sagepy_rescore_pin = write_sagepy_rescore_pin(
                P, P.sagepy_rescore_predictions
            )
            (
                P.sagepy_rescore_used_pin,
                P.sagepy_rescore_peptides,
                P.sagepy_rescore_psms,
            ) = mokapot(
                P,
                P.sagepy_rescore_pin,
                train_fdr=cfg.sagepy_rescore.get("train_fdr", 0.01),
                test_fdr=cfg.sagepy_rescore.get("test_fdr", 0.01),
                plugin="xgboost",
            )

        # FDR Summary
        P.sage_summary = sage_summarize(
            P, P.sage_results_tsv, P.sage_summarize_module, fdr=cfg.sage_summarize.fdr
        )

    # Exports -- final_mz_pmsms/final_precursors are the mz/rt/iim-corrected
    # outputs when recalibration ran (mode 2: mz only; mode 3: mz+RT+IIM),
    # otherwise the plain (uncorrected) MzPmsms/search_precursors from
    # materialize_pmsms_mz above (mode 1). Both must come from the same mode
    # SAGE1 actually searched against, or the exported headers and the SAGE
    # results disagree on what a peak's mz/rt/iim actually was -- see
    # plans/better_sage_filtering.md's B.6.
    P.search_mzml, P.search_mzml_idmap = convert_search_pmsms_to_mzml(
        P, final_mz_pmsms, final_precursors,
    )
    mgf_config_path = cfg.get("mgf", {}).get("config_path")
    if mgf_config_path:
        P.search_mgf = convert_search_pmsms_to_mgf(
            P,
            final_mz_pmsms,
            final_precursors,
            config_path=mgf_config_path,
        )

    if "fragpipe" in cfg:
        P.fragpipe_workflow = source_fragpipe_workflow(
            P, path=cfg.fragpipe.workflow_path
        )
        P.fragpipe_decoy_fasta = generate_fragpipe_decoy_fasta(P, P.fasta)
        P.fragpipe_workflow_patched = patch_fragpipe_workflow(
            P, P.fragpipe_workflow, P.fragpipe_decoy_fasta
        )
        P.fragpipe_manifest = write_fragpipe_manifest(P, P.search_mzml)
        P.fragpipe_results_dir = run_fragpipe(
            P,
            P.fragpipe_manifest,
            P.fragpipe_workflow_patched,
            ram=cfg.fragpipe.get("ram", 0),
        )
        P.fragpipe_log = extract_fragpipe_log(P, P.fragpipe_results_dir)
        P.fragpipe_summary = summarize_fragpipe(P, P.fragpipe_log)


def fragpipe_synthetic_pipeline(P: Pipeline, config: dict) -> None:
    """FragPipe smoke test on simulated data -- no Bruker .d input, no Sage.

    1. Simulate tryptic peptides from cfg.fasta_path (a small multi-protein FASTA)
       via Koina-predicted fragment ions, written as a pmsms dataset
       (`simulate_peptides_to_pmsms`).
    2. Convert that pmsms to a spectrum file: mzML (`convert_synthetic_pmsms_to_mzml`,
       via `git/pmsms2mzml`) or MGF (`convert_synthetic_pmsms_to_mgf`, via
       `venvs/common/bin/msms2mgf` + `[mgf].config_path`), picked by
       cfg.output_format ("mzml", the default, or "mgf").
    3. mzML only: run FragPipe on it (reusing ionmaiden_pipeline's FragPipe rules
       verbatim -- source_fragpipe_workflow/generate_fragpipe_decoy_fasta/
       write_fragpipe_manifest/run_fragpipe/extract_fragpipe_log/summarize_fragpipe).
       MGF has no search engine wired to it here (msms2mgf's plain-MGF output isn't
       what run_fragpipe's manifest step expects) -- P.synthetic_mgf is terminal.

    Verified end to end against a real FragPipe 24.0 + Philosopher install with a
    10-protein/~200-peptide synthetic FASTA: both convert_synthetic_pmsms_to_mzml
    (-> FragPipe, all 10 proteins recovered) and convert_synthetic_pmsms_to_mgf
    produced valid output from the same simulate_peptides_to_pmsms run.
    """
    cfg: DotDict = DotDict.Recursive(config)

    # Peptide Simulation
    P.fasta = source_fasta(P, path=cfg.fasta_path)
    P.synthetic_pmsms, P.synthetic_precursors = simulate_peptides_to_pmsms(
        P,
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
            P,
            P.synthetic_pmsms,
            P.synthetic_precursors,
        )

        # FragPipe Search
        P.fragpipe_workflow = source_fragpipe_workflow(
            P, path=cfg.fragpipe.workflow_path
        )
        P.fragpipe_decoy_fasta = generate_fragpipe_decoy_fasta(P, P.fasta)
        P.fragpipe_workflow_patched = patch_fragpipe_workflow(
            P, P.fragpipe_workflow, P.fragpipe_decoy_fasta
        )
        P.fragpipe_manifest = write_fragpipe_manifest(P, P.synthetic_mzml)
        P.fragpipe_results_dir = run_fragpipe(
            P, P.fragpipe_manifest, P.fragpipe_workflow_patched
        )
        P.fragpipe_log = extract_fragpipe_log(P, P.fragpipe_results_dir)
        P.fragpipe_summary = summarize_fragpipe(P, P.fragpipe_log)
    elif output_format == "mgf":
        mgf_config_path = cfg.get("mgf", {}).get("config_path")
        if not mgf_config_path:
            raise ValueError(
                "cfg.mgf.config_path is required when cfg.output_format == 'mgf'"
            )
        P.synthetic_mgf = convert_synthetic_pmsms_to_mgf(
            P,
            P.synthetic_pmsms,
            P.synthetic_precursors,
            config_path=mgf_config_path,
        )
    else:
        raise ValueError(
            f"cfg.output_format must be 'mzml' or 'mgf', got {output_format!r}"
        )
