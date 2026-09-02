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
class PipelineConfig(NodeType):
    filename = "pipeline_config.toml"


class Ms1Events(NodeType):
    filename = "events.ms1"


class Ms2Events(NodeType):
    filename = "events.ms2"


class Ms2TfsEvents(MmappetDataset):
    filename = "events_ms2_tfs.mmappet"


class Ms2TsfEvents(MmappetDataset):
    filename = "events_ms2_tsf.mmappet"


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


class PredictedRt(NodeType):
    """sequence -> rt parquet from `git/featureprediction`'s `predict_rt`
    -- the `--predicted-rt` cache SAGE1 needs to look up a candidate's
    predicted RT during search (it can't derive this itself). Independent
    of `PredictedIim` -- see plans/rt_iim_independent_dimensions.md
    (original combined design: plans/better_sage_filtering.md's B.4/B.5)."""

    filename = "predicted_rt.parquet"


class PredictedIim(NodeType):
    """sequence,charge -> iim parquet from `git/featureprediction`'s
    `predict_iim` -- the `--predicted-iim` cache. Independent of
    `PredictedRt`."""

    filename = "predicted_iim.parquet"


class RtPredictionCache(NodeType):
    """`mmappeteer.PredictionCache` directory (raw Chronologer HI, keyed by
    sequence) from `git/featureprediction`'s `fill_rt_cache` -- append-only,
    growing across every real (non-forced) rerun with the same
    `dumped_peptides` input. `mutable=True` for the same reason as
    `FragmentIntensityCache`: growth between runs must not invalidate
    `predict_rt`, which only ever reads it. Split out from `predict_rt`
    itself because necroflow's `mutable=True` requires a single-output rule
    and `predict_rt` has three outputs (`predicted_rt`/`rt_tolerance`/
    `plot`). Added 2026-09-01 -- before this, `predict_rt` had no way to
    pass `--cache-path` at all (the CLI flag didn't exist), so every job's
    RT prediction made a real, uncached Chronologer call for every sequence
    regardless of overlap with a previous job's dumped peptides."""

    filename = "rt_prediction_cache"


class IimPredictionCache(NodeType):
    """`mmappeteer.PredictionCache` directory (converted 1/K0, keyed by
    sequence+charge) from `git/featureprediction`'s `fill_iim_cache`. Same
    `mutable=True`/split rationale as `RtPredictionCache`. Kept as its own
    node (not sharing `RtPredictionCache`'s directory even though
    `cache.py`'s `PredictionCache` class supports both RT and IIM tables in
    one instance) so an RT-only job never has to depend on anything
    IIM/IM2Deep-shaped at all -- consistent with this pipeline's existing
    RT/IIM independence convention."""

    filename = "iim_prediction_cache"


class FragmentIntensityCache(NodeType):
    """`mmappeteer.PredictionCache` directory from `git/featureprediction`'s
    `fragment_intensity.py` -- append-only, growing across every real
    (non-forced) rerun of `predict_fragment_intensity` with the same
    inputs (lookup-before-append: a rerun makes no duplicate Koina calls
    for keys already cached). Its producer rule is `mutable=True`
    (necroflow's own "persistent single-output state whose external byte
    changes should not invalidate consumers" mechanism, `docs/rules.md`'s
    "Mutable Rules") specifically so that the cache growing between runs
    never stales anything downstream, and its workdir is exempt from
    necroflow's autoclean. See `plans/fragment_intensity_cache.md`."""

    filename = "fragment_intensity_cache"


class FragmentIntensityForSage(NodeType):
    """Job-scoped *index* from `git/featureprediction`'s
    `export_fragment_intensity_for_sage` -- `sequence, charge, start, end`
    pointers (not a copy of the sparse payload) for exactly this job's
    `dumped_peptides` x charge range, resolved against the (much bigger,
    shared, ever-growing) `FragmentIntensityCache`. A future SAGE reader
    keeps that cache's `arrays.mmappet` mmapped and uses this file only to
    look up which `[start, end)` range belongs to which `(sequence,
    charge)` -- see `git/featureprediction`'s AI.md. Not yet consumed by
    SAGE itself."""

    filename = "fragment_intensity_for_sage.parquet"


class NoPrediction(NodeType):
    """Zero-cost sentinel standing in for `RtTolerance`/`MobilityTolerance`
    when that dimension isn't active (`[recalibration.rt]`/
    `[recalibration.iim]` absent, 2026-08-25 -- see the pipeline factory).
    `update_sage_config_rt_iim` always takes both RT and IIM
    `rt_tolerance`/`mobility_tolerance` params (plain required `NodeType`s,
    not mixed Node/`None`); this stands in for whichever one wasn't
    computed. `run_sage`'s `predicted_rt`/`predicted_iim` used to need this
    too (`run_sage_with_predicted`, before necroflow's "mixed Node/value
    inputs" support) but no longer do -- `run_sage` now takes a true
    `PredictedRt | None`/`PredictedIim | None`, no sentinel required
    (2026-08-26). Only `update_sage_config_rt_iim`'s inputs still need this
    sentinel -- out of scope for that merge, see `run_sage`'s docstring.
    (necroflow *does* support skipping a rule entirely via plain `if/else`
    branching in the pipeline factory -- see its `docs/rules.md`'s
    "Conditional pipelines" -- but here RT/IIM are independently optional,
    so properly avoiding this sentinel too would mean up to 3 separate
    `update_sage_config_rt_iim` rule variants instead of one with a
    Python-callback command; not done yet, see
    plans/rt_iim_independent_dimensions.md.) A `@text_file` rule with zero
    `NodeType` inputs and a fixed literal `text` default -- same
    fingerprint every time, so it materializes once and every call site
    reuses the same node."""

    filename = "no_prediction.marker"


class RtTolerance(NodeType):
    """rt_tol_sec tolerance JSON (`ValueTolSpline`-shaped, wrapped as
    `{"rt_tol_sec": {...}}` -- see `feature_prediction.tolerance
    .write_value_tol_spline_json`). Produced by both `predict_rt`
    (pre-correction) and `correct_precursors_rt` (post-correction,
    tighter, the one actually fed to SAGE1 when RT correction runs)."""

    filename = "rt_tolerance.json"


class MobilityTolerance(NodeType):
    """Same as `RtTolerance`, for `mobility_tol`/IIM. Produced by both
    `predict_iim` and `correct_precursors_iim`."""

    filename = "mobility_tolerance.json"


class RtFitPlot(Png):
    filename = "rt_fit.png"


class IimFitPlot(Png):
    filename = "iim_fit.png"


class RtCorrectedPrecursors(RecalibratedPrecursors):
    """RecalibratedPrecursors (mz-corrected) with `rt` also corrected by
    `feature_prediction.precursor_correction.correct_precursors_rt` --
    `inv_ion_mobility` untouched. Accepted wherever RecalibratedPrecursors
    is, including `correct_precursors_iim`'s own `mz_corrected_precursors`
    input when both dimensions are active (chained, not recombined -- see
    plans/rt_iim_independent_dimensions.md)."""

    filename = "rt_corrected_precursors.mmappet"


class RtIimCorrectedPrecursors(RecalibratedPrecursors):
    """RecalibratedPrecursors (mz-corrected, optionally already rt-corrected
    via `RtCorrectedPrecursors` when both dimensions are active) with
    `inv_ion_mobility` also corrected by
    `feature_prediction.precursor_correction.correct_precursors_iim` --
    accepted wherever RecalibratedPrecursors (or, transitively,
    TofFilteredPrecursors/PreSageFilteredPrecursors) is, e.g. run_sage's
    final-pass `precursors` input in mode 3. Name kept from the pre-split
    design (the "final, however-corrected precursors" type) -- whether
    `rt` was actually corrected too depends on which precursors table was
    fed in as `mz_corrected_precursors`, not on this type alone. See
    plans/better_sage_filtering.md's B.6."""

    filename = "rt_iim_corrected_precursors.mmappet"


class PrecursorCorrectionRtModel(NodeType):
    filename = "precursor_correction_rt_model.json"


class PrecursorCorrectionIimModels(NodeType):
    """A directory of per-charge XGBoost models (`{charge}.json`, plus a
    pooled fallback `0.json`) -- see `precursor_correction.fit_iim_correction`."""

    filename = "precursor_correction_iim_models"


class PrecursorCorrectionRtFitPlot(Png):
    filename = "precursor_correction_rt_fit_plot.png"


class PrecursorCorrectionIimFitPlot(Png):
    filename = "precursor_correction_iim_fit_plot.png"


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
@text_file
def write_pipeline_config(text: str):
    config = output(PipelineConfig)
    return config


@command(
    "git/ionmaidenmetal/build/tdf2ms ms1 {tdf} {ms1}"
    " --threads {threads} --overwrite"
    " && test -f {ms1}/tof_row_starts.mmappet/schema.txt"
    " && test -f {ms1}/tof_urt_diff_index.mmappet/schema.txt"
    " && test -f {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt",
    threads=CORES,
)
def tdf2ms1(tdf: BrukerD):
    ms1 = output(Ms1Events)
    return ms1


@command("git/ionmaidenmetal/build/tdf2ms ms2 {tdf} {ms2} --overwrite")
def tdf2ms2(tdf: BrukerD):
    ms2 = output(Ms2Events)
    return ms2


@command(
    "git/ionmaidenmetal/build/tdf2ms ms2-tfs {tdf} {ms2_tfs}"
    " --threads {threads} --overwrite"
    " && test -f {ms2_tfs}/data.mmappet/schema.txt"
    " && test -f {ms2_tfs}/packed_scans.mmappet/schema.txt"
    " && test -f {ms2_tfs}/packed_scans.mmappet/shape.txt"
    " && test -f {ms2_tfs}/tof_index.mmappet/schema.txt"
    " && test -f {ms2_tfs}/tof_index.mmappet/shape.txt"
    " && test -f {ms2_tfs}/frame_index.mmappet/schema.txt"
    " && test -f {ms2_tfs}/stats.json",
    threads=CORES,
)
def tdf2ms2_tfs(tdf: BrukerD):
    ms2_tfs = output(Ms2TfsEvents)
    return ms2_tfs


@command(
    "git/ionmaidenmetal/build/tdf2ms ms2-tsf {tdf} {ms2_tsf}"
    " --threads {threads} --overwrite"
    " && test -f {ms2_tsf}/data.mmappet/schema.txt"
    " && test -f {ms2_tsf}/tof_scan_row_starts.mmappet/schema.txt"
    " && test -f {ms2_tsf}/tof_scan_row_starts.mmappet/shape.txt"
    " && test -f {ms2_tsf}/stats.json",
    threads=CORES,
)
def tdf2ms2_tsf(tdf: BrukerD):
    ms2_tsf = output(Ms2TsfEvents)
    return ms2_tsf


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


# Duplicated from `git/featureprediction`'s `fragment_intensity` module --
# can't import across venvs (same reasoning as `_DEFAULT_KOINA_*_SERVER_URL`
# above). Not yet exposed as job config: no real job has needed a different
# value, and `DEFAULT_COLLISION_ENERGY` is already documented as a
# placeholder, not a measured value -- see that repo's AI.md.
_DEFAULT_FRAGMENT_COLLISION_ENERGY = 30.0
_DEFAULT_FRAGMENT_FRAGMENTATION_TYPE = "HCD"


@command(
    "venvs/featureprediction/bin/feature-prediction-generate-fragments"
    " {dumped_peptides} {fragment_intensity_cache}"
    " --min-charge {min_charge} --max-charge {max_charge}"
    " --collision-energy {collision_energy} --fragmentation-type {fragmentation_type}",
    mutable=True,
)
def predict_fragment_intensity(
    dumped_peptides: DumpedPeptides,
    min_charge: int,
    max_charge: int,
    collision_energy: float,
    fragmentation_type: str,
):
    """Populate the persistent fragment-intensity `PredictionCache` for
    `dumped_peptides`' sequences x `[min_charge, max_charge]`.
    `mutable=True` (see `FragmentIntensityCache`'s docstring): the cache
    grows in place across reruns with the same inputs rather than
    starting fresh each time, and its workdir is exempt from autoclean.
    Independently requestable (like `ms2_tsf_events`/`ms2_tfs_events`) --
    no downstream consumer in the DAG yet (SAGE doesn't read this cache),
    and it's a real, expensive, network-bound operation (a full
    human-proteome fill took ~43 minutes against the live Koina server),
    so it must never run just because `dumped_peptides` exists -- only
    when explicitly requested via `.requests`. See
    `plans/fragment_intensity_cache.md`.
    """
    fragment_intensity_cache = output(FragmentIntensityCache)
    return fragment_intensity_cache


@command(
    "venvs/featureprediction/bin/feature-prediction-export-fragments-for-sage"
    " {dumped_peptides} {fragment_intensity_cache} {fragment_intensity_for_sage}"
    " --min-charge {min_charge} --max-charge {max_charge} --collision-energy {collision_energy}"
)
def export_fragment_intensity_for_sage(
    dumped_peptides: DumpedPeptides,
    fragment_intensity_cache: FragmentIntensityCache,
    min_charge: int,
    max_charge: int,
    collision_energy: float,
):
    """Scope `fragment_intensity_cache` (shared, ever-growing across every
    job) down to exactly this job's `dumped_peptides` x `[min_charge,
    max_charge]`, via `git/featureprediction`'s DuckDB-based export -- see
    that repo's AI.md, "`export_fragment_intensity.py`". Ordinary
    (non-mutable) rule despite depending on a `mutable=True` parent: per
    `docs/rules.md`'s "Mutable Rules", "if a mutable call executes during
    the current run, every consumer replays" -- so this always sees the
    cache state left by `predict_fragment_intensity`'s own most recent
    (real, non-forced) run in the same invocation, never a stale one.
    Independently requestable, same reasoning as `predict_fragment_intensity`
    itself -- no downstream consumer in the DAG yet (SAGE doesn't read this
    export), so it must never run just because `dumped_peptides`/
    `fragment_intensity_cache` exist, only when explicitly requested via
    `.requests`.
    """
    fragment_intensity_for_sage = output(FragmentIntensityForSage)
    return fragment_intensity_for_sage


def _run_sage_command(args: CommandArgs) -> str:
    """Python command callback, not a static template -- lets
    `--predicted-rt`/`--predicted-iim` be added independently, purely by
    whether `predicted_rt`/`predicted_iim` resolved to a real Node or the
    mixed Node/`None` input's `None` default (necroflow's "mixed Node/value
    inputs" support, added 2026-08-21) -- `args.inputs.*` preserves a plain
    `None` verbatim (only managed Nodes resolve to `Path`), so this needs
    neither a separate `dimensions` scalar nor the `NoPrediction` sentinel
    `run_sage_with_predicted` used to require just to keep the parameter a
    real DAG edge. One rule now covers pass-1/mode-1/mode-2 (neither
    prediction) and mode 3's final pass (either/both), replacing the
    previous `run_sage`/`run_sage_with_predicted` split. See
    plans/rt_iim_independent_dimensions.md for the split this undoes.
    """
    sage_binary = shlex.quote(str(args.inputs.sage_binary))
    fasta = shlex.quote(str(args.inputs.fasta))
    pmsms = shlex.quote(str(args.inputs.pmsms))
    precursors = shlex.quote(str(args.inputs.precursors))
    sage_config = shlex.quote(str(args.inputs.sage_config))
    workdir = shlex.quote(str(args.workdir))
    results_json = shlex.quote(str(args.outputs.results_json))
    results_pin = shlex.quote(str(args.outputs.results_pin))
    results_tsv = shlex.quote(str(args.outputs.results_tsv))
    matched_fragments = shlex.quote(str(args.outputs.matched_fragments))

    flags = ""
    if args.inputs.predicted_rt is not None:
        flags += f" --predicted-rt {shlex.quote(str(args.inputs.predicted_rt))}"
    if args.inputs.predicted_iim is not None:
        flags += f" --predicted-iim {shlex.quote(str(args.inputs.predicted_iim))}"
    # Both-or-neither, mirroring Rust's own `Input::build()` validation --
    # `predicted_fragment_intensity_cache` is the *directory* the
    # PredictionCache node produced; the Rust reader only ever wants its
    # `arrays.mmappet` subdirectory (never `index.sqlite3`/`write.lock`),
    # see `docs/ai/predicted_fragment_intensity.md`.
    if args.inputs.predicted_fragment_intensity_index is not None:
        frag_index = shlex.quote(str(args.inputs.predicted_fragment_intensity_index))
        frag_cache = shlex.quote(
            str(args.inputs.predicted_fragment_intensity_cache / "arrays.mmappet")
        )
        flags += (
            f" --predicted-fragment-intensity-index {frag_index}"
            f" --predicted-fragment-intensity-cache {frag_cache}"
        )

    return (
        f"{sage_binary} --version && {sage_binary} -f {fasta}"
        f" --annotate-matches --write-pin --output_directory {workdir}"
        f" --pmsms {pmsms} --precursors {precursors}{flags} {sage_config}"
        f" && test -f {results_json} && test -f {results_pin}"
        f" && test -f {results_tsv} && test -f {matched_fragments}"
    )


@command(_run_sage_command, threads=CORES)
def run_sage(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
    sage_binary: SageBinary,
    predicted_rt: PredictedRt | None = None,
    predicted_iim: PredictedIim | None = None,
    predicted_fragment_intensity_index: FragmentIntensityForSage | None = None,
    predicted_fragment_intensity_cache: FragmentIntensityCache | None = None,
):
    """Run Sage. `predicted_rt`/`predicted_iim` are optional (mixed
    Node/`None` inputs) -- omitted for pass-1 and mode 1/2's plain search,
    passed as real Nodes for whichever dimension mode 3 has active.
    `sage_config` must already carry `rt_tol_sec`+`rt_sigma_sec` and/or
    `mobility_tol`+`iim_sigma` for whichever dimensions are passed here
    (`update_sage_config_rt_iim`), same requirement `Input::build()`
    enforces Rust-side. See plans/better_sage_filtering.md's B.4-B.6 and
    plans/rt_iim_independent_dimensions.md.

    `predicted_fragment_intensity_index`/`_cache` are likewise both-or-
    neither (mixed Node/`None`), gated by `"fragment_intensity" in cfg` at
    the pipeline-factory call site, not by anything in this function --
    feature-only (`ms2_*` scoring columns), no hard eviction, independent
    of predicted_rt/predicted_iim. See
    `git/sage/docs/ai/predicted_fragment_intensity.md`.
    """
    results_json = output(SageResultsJson)
    results_pin = output(SageResultsPin)
    results_tsv = output(SageResultsTsv)
    matched_fragments = output(SageMatchedFragments)
    return results_json, results_pin, results_tsv, matched_fragments


@text_file
def write_no_prediction_marker(text: str):
    marker = output(NoPrediction)
    return marker


_NO_PREDICTION_TEXT = "no prediction requested for this dimension\n"

# Duplicated from `git/featureprediction`'s `koina_client` -- can't import
# across venvs (this repo shells out to a separate `venvs/featureprediction`
# install, doesn't depend on that package directly), same reasoning as
# `default_precursor_charge` duplicating SAGE's own compiled-in (2, 4)
# default just above.
#
# Two different defaults, not one shared value: `predict_iim` still goes
# through `koinapy.Koina` (IM2Deep), which speaks **gRPC** on :8500;
# `predict_rt`'s direct-HTTP Chronologer client needs the *separate* HTTP
# REST port, :8501 -- verified live (2026-08-25) that POSTing plain HTTP to
# :8500 gets back raw HTTP/2 SETTINGS-frame bytes (`BadStatusLine`), since
# gRPC always runs over HTTP/2; :8501 is real Triton HTTP on the same host.
_DEFAULT_KOINA_GRPC_SERVER_URL = "192.168.1.222:8500"  # predict_iim (koinapy/IM2Deep)
_DEFAULT_KOINA_HTTP_SERVER_URL = "192.168.1.222:8501"  # predict_rt (direct HTTP/Chronologer)


def _server_url_arg(value: str | list[str] | None, default: str) -> str:
    """`[recalibration.rt]`/`[recalibration.iim]`'s `server_url` (a plain
    string or a TOML array -- either is valid config shape) into the single
    comma-separated string `feature-prediction-generate-{rt,iim}`'s
    `--server-url` flag expects (it splits on comma back into a priority
    list -- ip0 tried first, falling back to ip1 etc. on failure). `default`
    must match the caller's own protocol/port (see the two constants
    above -- RT and IIM are not interchangeable here)."""
    if value is None:
        return default
    if isinstance(value, str):
        return value
    return ",".join(value)




@command(
    "venvs/featureprediction/bin/feature-prediction-fill-rt-cache"
    " {dumped_peptides} {rt_prediction_cache} --server-url {server_url}",
    mutable=True,
)
def fill_rt_prediction_cache(dumped_peptides: DumpedPeptides, server_url: str):
    """Populate the RT-prediction cache for `dumped_peptides`' sequences.
    `mutable=True` for the same reason as `predict_fragment_intensity` --
    see `RtPredictionCache`'s docstring. Independently requestable, same
    reasoning as that rule too: this is a real, network-bound Koina call
    (though far cheaper than the fragment-intensity fill -- Chronologer's
    direct-HTTP path runs the full F9477 dump in ~90s), so it must never
    run just because `dumped_peptides` exists, only when `predict_rt`
    actually needs it.
    """
    rt_prediction_cache = output(RtPredictionCache)
    return rt_prediction_cache


@command(
    "venvs/featureprediction/bin/feature-prediction-generate-rt"
    " {dumped_peptides} {sage_results_tsv} {predicted_rt} {rt_tolerance} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --server-url {server_url} --fdr {fdr}"
    " --cache-path {rt_prediction_cache}"
)
def predict_rt(
    dumped_peptides: DumpedPeptides,
    sage_results_tsv: SageResultsTsv,
    rt_prediction_cache: RtPredictionCache,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    server_url: str,
    fdr: int | float,
):
    """`rt_prediction_cache` should be `fill_rt_prediction_cache`'s output,
    filled with the same `dumped_peptides` beforehand -- given that, this
    call becomes a pure cache lookup (`predict_hi_cached` finds zero
    missing sequences), no real Koina call. Still works correctly (just
    slower, filling on demand) if pointed at a cold/partial cache."""
    predicted_rt = output(PredictedRt)
    rt_tolerance = output(RtTolerance)
    plot = output(RtFitPlot)
    return predicted_rt, rt_tolerance, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-fill-iim-cache"
    " {dumped_peptides} {iim_prediction_cache}"
    " --min-charge {min_charge} --max-charge {max_charge} --server-url {server_url}",
    mutable=True,
)
def fill_iim_prediction_cache(
    dumped_peptides: DumpedPeptides, min_charge: int, max_charge: int, server_url: str
):
    """Populate the IIM-prediction cache for `dumped_peptides`' sequences x
    `[min_charge, max_charge]`. Same rationale as `fill_rt_prediction_cache`
    -- see `IimPredictionCache`'s docstring. This one is the genuinely
    expensive Koina call (IM2Deep, synchronous, historically the ~45-minute
    bottleneck on a full F9477 dump), so caching it properly matters more
    than for RT.
    """
    iim_prediction_cache = output(IimPredictionCache)
    return iim_prediction_cache


@command(
    "venvs/featureprediction/bin/feature-prediction-generate-iim"
    " {dumped_peptides} {sage_results_tsv} {predicted_iim} {mobility_tolerance} {plot}"
    " --min-charge {min_charge} --max-charge {max_charge}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --server-url {server_url} --fdr {fdr}"
    " --cache-path {iim_prediction_cache}"
)
def predict_iim(
    dumped_peptides: DumpedPeptides,
    sage_results_tsv: SageResultsTsv,
    iim_prediction_cache: IimPredictionCache,
    min_charge: int,
    max_charge: int,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    server_url: str,
    fdr: int | float,
):
    """`iim_prediction_cache` should be `fill_iim_prediction_cache`'s output,
    filled with the same `dumped_peptides`/charge range beforehand -- same
    "becomes a pure cache lookup" reasoning as `predict_rt`'s docstring."""
    predicted_iim = output(PredictedIim)
    mobility_tolerance = output(MobilityTolerance)
    plot = output(IimFitPlot)
    return predicted_iim, mobility_tolerance, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-correct-precursors-rt"
    " {sage_results_tsv} {predicted_rt} {mz_corrected_precursors}"
    " {output_precursors} {rt_tolerance} {rt_model} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}"
)
def correct_precursors_rt(
    sage_results_tsv: SageResultsTsv,
    predicted_rt: PredictedRt,
    mz_corrected_precursors: RecalibratedPrecursors,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    output_precursors = output(RtCorrectedPrecursors)
    rt_tolerance = output(RtTolerance)
    rt_model = output(PrecursorCorrectionRtModel)
    plot = output(PrecursorCorrectionRtFitPlot)
    return output_precursors, rt_tolerance, rt_model, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-correct-precursors-iim"
    " {sage_results_tsv} {predicted_iim} {mz_corrected_precursors}"
    " {output_precursors} {mobility_tolerance} {iim_models} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}"
)
def correct_precursors_iim(
    sage_results_tsv: SageResultsTsv,
    predicted_iim: PredictedIim,
    mz_corrected_precursors: RecalibratedPrecursors,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    """`mz_corrected_precursors` accepts a plain `RecalibratedPrecursors`
    (IIM-only mode) or an `RtCorrectedPrecursors` (both dimensions active,
    chained -- `RtCorrectedPrecursors` is a `RecalibratedPrecursors`
    subclass, satisfies this contract either way). See
    plans/rt_iim_independent_dimensions.md."""
    output_precursors = output(RtIimCorrectedPrecursors)
    mobility_tolerance = output(MobilityTolerance)
    iim_models = output(PrecursorCorrectionIimModels)
    plot = output(PrecursorCorrectionIimFitPlot)
    return output_precursors, mobility_tolerance, iim_models, plot


def _update_sage_config_rt_iim_command(args: CommandArgs) -> str:
    """Python command callback -- chains 0, 2, or 4 `config_set` calls
    based on `args.config.dimensions` (RT/IIM independently optional).
    Unlike `run_sage` (which now takes `predicted_rt`/`predicted_iim` as
    true mixed Node/`None` inputs, no sentinel needed), `rt_tolerance`/
    `mobility_tolerance` here are still always real DAG edges (resolve to
    the `NoPrediction` sentinel when inactive) -- out of scope for the
    2026-08-26 `run_sage` merge, see that function's docstring. Only the active
    dimension's `config_set` calls actually run -- decided by the scalar
    `dimensions` config, never by inspecting file content. In practice
    `dimensions` is never empty here (the pipeline factory only calls this
    rule inside `if "rt" in cfg.recalibration or "iim" in cfg.recalibration:`,
    which always sets at least one dimension), but a plain copy-through is
    a safe fallback.

    Each active dimension contributes *two* steps, not one: `rt_tol_sec`/
    `mobility_tol` (the hard-eviction window) and `rt_sigma_sec`/`iim_sigma`
    (the LDA/ranking z² scale) -- both read from the same `rt_tolerance`/
    `mobility_tolerance` artifact (`git/featureprediction`'s
    `correct_precursors_rt`/`correct_precursors_iim` write them as sibling
    top-level keys in one file, see that repo's `AI.md`), since SAGE's
    `Input::build` rejects `predicted_rt`/`rt_tol_sec` being set without
    `rt_sigma_sec` (and vice versa for IIM) -- see `plans/lda_external_rt_iim_features.md`.
    """
    sage_config = shlex.quote(str(args.inputs.sage_config))
    recalibrated_sage_config = shlex.quote(str(args.outputs.recalibrated_sage_config))
    dimensions = set(args.config.dimensions)

    current = sage_config
    steps = []
    if "rt" in dimensions:
        rt_tolerance = shlex.quote(str(args.inputs.rt_tolerance))
        tol_out = shlex.quote(str(args.workdir / "rt_tol_updated.json"))
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {tol_out}"
            f" --target rt_tol_sec --source {rt_tolerance} --source-field rt_tol_sec"
        )
        current = tol_out
        sigma_out = shlex.quote(str(args.workdir / "rt_sigma_updated.json"))
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {sigma_out}"
            f" --target rt_sigma_sec --source {rt_tolerance} --source-field rt_sigma_sec"
        )
        current = sigma_out
    if "iim" in dimensions:
        mobility_tolerance = shlex.quote(str(args.inputs.mobility_tolerance))
        tol_out = shlex.quote(str(args.workdir / "mobility_tol_updated.json"))
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {tol_out}"
            f" --target mobility_tol --source {mobility_tolerance} --source-field mobility_tol"
        )
        current = tol_out
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {recalibrated_sage_config}"
            f" --target iim_sigma --source {mobility_tolerance} --source-field iim_sigma"
        )
        current = recalibrated_sage_config
    if current != recalibrated_sage_config:
        steps.append(f"cp {current} {recalibrated_sage_config}")
    return " && ".join(steps)


@command(_update_sage_config_rt_iim_command)
def update_sage_config_rt_iim(
    sage_config: SageConfig,
    rt_tolerance: RtTolerance | NoPrediction,
    mobility_tolerance: MobilityTolerance | NoPrediction,
    dimensions: tuple[str, ...],
):
    """Chains onto `update_sage_config`'s output (mz's `precursor_tol`/
    `fragment_tol` already patched) -- `sage_config` here is that rule's
    `recalibrated_sage_config`, not the plain `write_sage_config` output.
    `--predicted-rt`/`--predicted-iim` themselves are *not* set here: like
    `--pmsms`/`--precursors`/`-f {fasta}`, they're path-valued overrides
    passed directly as CLI flags (`run_sage`), not embedded into the
    config JSON via `config_set`."""
    recalibrated_sage_config = output(SageConfig)
    return recalibrated_sage_config


def _mokapot_command(args: CommandArgs) -> str:
    """Python command callback, not a static template -- lets `--plugin`
    be added conditionally (empty for the plain-Sage-PIN call, `--plugin
    xgboost` for the sagepy-rescore call) without necroflow's string-
    template placeholders needing to express a conditional substring.
    Same reasoning for `--mode`/`--rt-source`/`--iim-source`: only the
    plain-Sage-PIN call passes `rt_source`/`iim_source`, which selects
    `scripts/mokapot_pin_adapter.py --mode sage` (its leakage-safe feature
    registry, see `plans/mokapot_leakage_safe_pin.md`); the sagepy-rescore
    call passes neither, so the adapter defaults to `--mode passthrough`
    (its original, unchanged, FileName-drop-only behavior) -- that PIN is
    already leakage-filtered upstream by
    `sagepy_rescore.features.build_feature_frame`.
    """
    pin = shlex.quote(str(args.inputs.pin))
    used_pin = shlex.quote(str(args.outputs.used_pin))
    peptides = shlex.quote(str(args.outputs.peptides))
    psms = shlex.quote(str(args.outputs.psms))
    workdir = shlex.quote(str(args.workdir))
    plugin_flag = f" --plugin {args.config.plugin}" if args.config.plugin else ""
    adapter_flags = ""
    if args.config.rt_source is not None or args.config.iim_source is not None:
        rt_source = shlex.quote(args.config.rt_source or "none")
        iim_source = shlex.quote(args.config.iim_source or "none")
        adapter_flags = f" --mode sage --rt-source {rt_source} --iim-source {iim_source}"
    return (
        f"venvs/mokapot/bin/python scripts/mokapot_pin_adapter.py -i {pin} -o {used_pin}"
        f"{adapter_flags}"
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
    rt_source: str | None = None,
    iim_source: str | None = None,
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
    P.pipeline_config = write_pipeline_config(P, text=tomlkit.dumps(config))

    # Acquisition
    P.tdf = source_bruker_d(P, path=cfg.tdf_path)
    P.fasta = source_fasta(P, path=cfg.fasta_path)
    P.mkpmsms_binary = source_mkpmsms_binary(
        P, path="git/ionmaidenmetal/build/mkpmsms"
    )
    P.sage_summarize_module = source_sage_summarize_module(
        P, path="git/searchops/src/searchops/sage.py"
    )
    P.sage_binary = source_sage_binary(P, path="git/sage/target/release/sage")
    P.dump_peptides_binary = source_dump_peptides_binary(
        P, path="git/sage/target/release/dump_peptides"
    )

    # Raw Extraction
    P.ms1_events = tdf2ms1(P, P.tdf)
    P.ms2_events = tdf2ms2(P, P.tdf)
    P.ms2_tfs_events = tdf2ms2_tfs(P, P.tdf)
    P.ms2_tsf_events = tdf2ms2_tsf(P, P.tdf)

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

        # Independently requestable (see `predict_fragment_intensity`'s
        # docstring) -- min_charge/max_charge mirror the same
        # cfg.sage.precursor_charge-or-SAGE's-own-(2,4)-default derivation
        # `[recalibration.iim]` uses below, computed here too since this
        # runs unconditionally (not nested inside `"recalibration" in cfg`).
        _fragment_min_charge, _fragment_max_charge = cfg.sage.get("precursor_charge", (2, 4))
        P.predicted_fragment_intensity = predict_fragment_intensity(
            P,
            P.dumped_peptides,
            min_charge=_fragment_min_charge,
            max_charge=_fragment_max_charge,
            collision_energy=_DEFAULT_FRAGMENT_COLLISION_ENERGY,
            fragmentation_type=_DEFAULT_FRAGMENT_FRAGMENTATION_TYPE,
        )
        # Same independently-requestable reasoning (see
        # `export_fragment_intensity_for_sage`'s docstring) -- reuses the
        # same charge range/collision_energy the cache above was filled
        # with, so a request for this never sees a coverage gap from a
        # mismatched range.
        P.fragment_intensity_for_sage = export_fragment_intensity_for_sage(
            P,
            P.dumped_peptides,
            P.predicted_fragment_intensity,
            min_charge=_fragment_min_charge,
            max_charge=_fragment_max_charge,
            collision_energy=_DEFAULT_FRAGMENT_COLLISION_ENERGY,
        )

        # `[fragment_intensity]` (presence-only table) is the on/off switch
        # for actually feeding the cache above into search -- without this
        # gate, threading `P.fragment_intensity_for_sage`/
        # `P.predicted_fragment_intensity` into `run_sage` unconditionally
        # would make *every* sage job transitively depend on
        # `predict_fragment_intensity`, breaking its own documented
        # "never runs just because `dumped_peptides` exists" invariant.
        # Mirrors `"recalibration" in cfg`/`"rt"/"iim" in cfg.recalibration`'s
        # own table-presence-as-flag convention.
        if "fragment_intensity" in cfg:
            _final_pass_fragment_intensity_index = P.fragment_intensity_for_sage
            _final_pass_fragment_intensity_cache = P.predicted_fragment_intensity
        else:
            _final_pass_fragment_intensity_index = None
            _final_pass_fragment_intensity_cache = None

        def _finalize_confident_psms(search_precursors, mz_pmsms, search_pmsms):
            """confident_psms -> sage_pmsms_mapping -> score_comparison, the
            same three calls needed after any final `run_sage` call
            (mode-1/2/3 alike, see `recalibration_modes.md`) -- only which
            precursors/pmsms Nodes get passed differs: raw `search_*` Nodes
            when there's no recalibration at all, the recalibrated ones
            otherwise. Kept as a plain closure over `P`/`cfg`, not a new
            necroflow rule -- it only groups three existing rule calls."""
            P.confident_psms = filter_sage_results(
                P, P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
            )
            P.sage_pmsms_mapping = sage_map_to_pmsms(
                P,
                P.confident_psms,
                P.sage_matched_fragments,
                search_precursors,
                mz_pmsms,
            )
            P.score_comparison = score_comparison(
                P,
                search_precursors,
                search_pmsms,
                P.sage_pmsms_mapping,
                P.pseudomsms_config,
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

            # RT and IIM are independent, sibling optional steps, each gated
            # purely on its own `[recalibration.rt]`/`.iim]` table presence
            # (2026-08-25 rationale -- mirrors mz's own `"recalibration" in
            # cfg` pattern, no separate flag needed; see
            # plans/rt_iim_independent_dimensions.md). `dimensions` may end
            # up empty (neither active -- previously called "mode 2") --
            # `update_sage_config_rt_iim` below tolerates that as a
            # documented safe copy-through (its own command builder's
            # docstring), so mz-only and RT/IIM-active jobs share the same
            # single final `run_sage` call below instead of two
            # near-duplicate branches (2026-09-0x flattening; original
            # two-branch design: plans/better_sage_filtering.md's B.6).
            dimensions = tuple(d for d in ("rt", "iim") if d in cfg.recalibration)

            # tolerance_percentiles/tolerance_method: required explicitly in
            # `[recalibration.rt]`/`[recalibration.iim]` (separate tables,
            # 2026-08-25) -- never inherited from `cfg.recalibration`'s own
            # (mz-scoped) `[recalibration.mz]` value; each of the three
            # dimensions (mz, rt, iim) gets its own percentiles/method. See
            # `git/featureprediction`'s `tolerance.select_tolerance`.
            precursor_correction_rt_tolerance = None
            precursor_correction_mobility_tolerance = None
            precursors_for_iim_correction = P.recalibrated_precursors
            final_precursors = P.recalibrated_precursors

            if "rt" in dimensions:
                rt_tolerance_lo, rt_tolerance_hi = cfg.recalibration.rt["tolerance_percentiles"]
                rt_tolerance_method = cfg.recalibration.rt.get("tolerance_method", "theoretic")
                rt_server_url = _server_url_arg(
                    cfg.recalibration.rt.get("server_url"), _DEFAULT_KOINA_HTTP_SERVER_URL
                )
                P.rt_prediction_cache = fill_rt_prediction_cache(
                    P, P.dumped_peptides, server_url=rt_server_url
                )
                P.predicted_rt, P.rt_tolerance, P.rt_fit_plot = predict_rt(
                    P,
                    P.dumped_peptides,
                    P.filtered_sage_results_tsv,
                    P.rt_prediction_cache,
                    tolerance_lo=rt_tolerance_lo,
                    tolerance_hi=rt_tolerance_hi,
                    tolerance_method=rt_tolerance_method,
                    server_url=rt_server_url,
                    fdr=cfg.sage_summarize.fdr,
                )
                (
                    P.rt_corrected_precursors,
                    precursor_correction_rt_tolerance,
                    P.precursor_correction_rt_model,
                    P.precursor_correction_rt_fit_plot,
                ) = correct_precursors_rt(
                    P,
                    P.filtered_sage_results_tsv,
                    P.predicted_rt,
                    P.recalibrated_precursors,
                    tolerance_lo=rt_tolerance_lo,
                    tolerance_hi=rt_tolerance_hi,
                    tolerance_method=rt_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
                )
                precursors_for_iim_correction = P.rt_corrected_precursors
                final_precursors = P.rt_corrected_precursors
            else:
                # Plain Python `None`, not a `NoPrediction` sentinel node --
                # `run_sage`'s `predicted_rt` is a true mixed Node/`None`
                # input, so there's no DAG edge to satisfy when this
                # dimension is inactive.
                P.predicted_rt = None

            if "iim" in dimensions:
                # min_charge/max_charge: explicit `[recalibration.iim]`
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
                rt_iim_min_charge = cfg.recalibration.iim.get(
                    "min_charge", default_precursor_charge[0]
                )
                rt_iim_max_charge = cfg.recalibration.iim.get(
                    "max_charge", default_precursor_charge[1]
                )
                iim_tolerance_lo, iim_tolerance_hi = cfg.recalibration.iim["tolerance_percentiles"]
                iim_tolerance_method = cfg.recalibration.iim.get("tolerance_method", "theoretic")
                iim_server_url = _server_url_arg(
                    cfg.recalibration.iim.get("server_url"), _DEFAULT_KOINA_GRPC_SERVER_URL
                )
                P.iim_prediction_cache = fill_iim_prediction_cache(
                    P,
                    P.dumped_peptides,
                    min_charge=rt_iim_min_charge,
                    max_charge=rt_iim_max_charge,
                    server_url=iim_server_url,
                )
                P.predicted_iim, P.mobility_tolerance, P.iim_fit_plot = predict_iim(
                    P,
                    P.dumped_peptides,
                    P.filtered_sage_results_tsv,
                    P.iim_prediction_cache,
                    min_charge=rt_iim_min_charge,
                    max_charge=rt_iim_max_charge,
                    tolerance_lo=iim_tolerance_lo,
                    tolerance_hi=iim_tolerance_hi,
                    tolerance_method=iim_tolerance_method,
                    server_url=iim_server_url,
                    fdr=cfg.sage_summarize.fdr,
                )
                # Chained, not recombined: correct_precursors_iim's
                # mz_corrected_precursors input becomes correct_precursors_rt's
                # own output when both dimensions are active, so a single
                # final precursors table ends up with both corrections
                # applied (RtCorrectedPrecursors is a RecalibratedPrecursors
                # subclass, satisfies either call site's contract).
                (
                    P.rt_iim_corrected_precursors,
                    precursor_correction_mobility_tolerance,
                    P.precursor_correction_iim_models,
                    P.precursor_correction_iim_fit_plot,
                ) = correct_precursors_iim(
                    P,
                    P.filtered_sage_results_tsv,
                    P.predicted_iim,
                    precursors_for_iim_correction,
                    tolerance_lo=iim_tolerance_lo,
                    tolerance_hi=iim_tolerance_hi,
                    tolerance_method=iim_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
                )
                final_precursors = P.rt_iim_corrected_precursors
            else:
                P.predicted_iim = None

            # update_sage_config_rt_iim needs a real RtTolerance/
            # MobilityTolerance for whichever dimension's correction ran
            # (tighter, post-correction fit) -- the sentinel for whichever
            # dimension is inactive (including both, when `dimensions` is
            # empty). Unlike predicted_rt/predicted_iim above,
            # update_sage_config_rt_iim's rt_tolerance/mobility_tolerance
            # inputs are still plain required NodeTypes, not mixed
            # Node/None -- out of scope for this refactor -- so the
            # sentinel stays here.
            P.precursor_correction_rt_tolerance = (
                precursor_correction_rt_tolerance
                if precursor_correction_rt_tolerance is not None
                else write_no_prediction_marker(P, text=_NO_PREDICTION_TEXT)
            )
            P.precursor_correction_mobility_tolerance = (
                precursor_correction_mobility_tolerance
                if precursor_correction_mobility_tolerance is not None
                else write_no_prediction_marker(P, text=_NO_PREDICTION_TEXT)
            )

            # Unconditional -- a safe no-op copy-through when `dimensions`
            # is empty (`_update_sage_config_rt_iim_command`'s own
            # docstring documents this fallback: `cp {sage_config}
            # {recalibrated_sage_config}`, nothing else), so mz-only jobs
            # and RT/IIM-active jobs share this one call and the one final
            # `run_sage` call below.
            P.recalibrated_sage_config_rt_iim = update_sage_config_rt_iim(
                P,
                P.recalibrated_sage_config,
                P.precursor_correction_rt_tolerance,
                P.precursor_correction_mobility_tolerance,
                dimensions=dimensions,
            )
            (
                P.sage_results_json,
                P.sage_results_pin,
                P.sage_results_tsv,
                P.sage_matched_fragments,
            ) = run_sage(
                P,
                P.recalibrated_mz_pmsms,
                final_precursors,
                P.fasta,
                P.recalibrated_sage_config_rt_iim,
                P.sage_binary,
                predicted_rt=P.predicted_rt,
                predicted_iim=P.predicted_iim,
                predicted_fragment_intensity_index=_final_pass_fragment_intensity_index,
                predicted_fragment_intensity_cache=_final_pass_fragment_intensity_cache,
            )

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
            _finalize_confident_psms(
                P.search_precursors, P.recalibrated_mz_pmsms, P.search_pmsms
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
                predicted_fragment_intensity_index=_final_pass_fragment_intensity_index,
                predicted_fragment_intensity_cache=_final_pass_fragment_intensity_cache,
            )
            _finalize_confident_psms(
                P.search_precursors, P.search_mz_pmsms, P.search_pmsms
            )

        # getattr, not `P.predicted_rt` directly -- the non-recalibration
        # branch above never binds these labels at all (only the
        # recalibration branch's own rt/iim sub-blocks do, each `None` when
        # that dimension is inactive), and `Pipeline.__getattr__` raises
        # `AttributeError` for an unbound label, not a Python-level
        # default. `rt_source`/`iim_source` "external" exactly when that
        # dimension's real external-prediction Node exists for this job --
        # see `plans/mokapot_leakage_safe_pin.md`.
        # `[mokapot].plugin` (optional, e.g. `plugin = "xgboost"`) -- unlike
        # the sagepy_rescore branch below (always xgboost, unrelated call),
        # this call site's model choice is config-driven so a job can pick
        # mokapot's own default (linear SVM) or xgboost.
        P.mokapot_used_pin, P.mokapot_peptides, P.mokapot_psms = mokapot(
            P,
            P.sage_results_pin,
            # "" (not None) when unset -- necroflow's dependency-provenance
            # recording can't serialize a bare `None` scalar kwarg value to
            # TOML; "" is falsy in `_mokapot_command`'s own `if
            # args.config.plugin` check, so the effect (no --plugin flag,
            # mokapot's own default model) is identical.
            plugin=(cfg.mokapot.get("plugin", "") if "mokapot" in cfg else ""),
            rt_source=("external" if getattr(P, "predicted_rt", None) is not None else "none"),
            iim_source=("external" if getattr(P, "predicted_iim", None) is not None else "none"),
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
    P.pipeline_config = write_pipeline_config(P, text=tomlkit.dumps(config))

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
