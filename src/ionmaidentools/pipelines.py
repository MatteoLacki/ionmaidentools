"""Necroflow ionmaiden pipeline.

`run_sage` and `mokapot` count their confident spectra, peptidoforms and ions as the last
step of their own commands (`sage-fdr-levels`/`mokapot-fdr-levels`, git/searchops) and,
when `ION_MAIDEN_REGRESSION_DB` is set, record the run there (git/searchops
`docs/ai/run_db.md`).

Fingerprint v3 used.
"""

from __future__ import annotations

import json
import os
import shlex
from pathlib import Path

import tomlkit

from dictodot import DotDict
from necroflow import (
    Node,
    Pipeline,
    Resource,
    command,
    outputs,
    shell,
    symlink,
    text_file,
    workflow,
)

CORES = os.cpu_count() or 1


# --- base node types ---
class MmappetDataset(Node):
    """Base type for outputs that are mmappet directories (see CLAUDE.md's
    "Precursor table format" convention). No `filename` of its own -- every
    concrete mmappet output subclasses this and sets its own."""


class Json(Node):
    """A JSON file."""


class Toml(Node):
    """A TOML file."""


class Tsv(Node):
    """A tab-separated text file."""


class Png(Node):
    """A PNG image."""


class Executable(Node):
    """A program the pipeline runs."""


# --- source node types ---
class BrukerD(Node):
    filename = "input.d"


class Fasta(Node):
    filename = "fasta.fasta"


class MkpmsmsBinary(Executable):
    filename = "mkpmsms"


class SageBinary(Executable):
    filename = "sage"


class SageSummarizeModule(Node):
    filename = "sage.py"


class Tdf2msBinary(Executable):
    filename = "tdf2ms"


class PrecursorNeighborsCsrBinary(Executable):
    filename = "precursor_neighbors_csr"


class TofFilterBinary(Executable):
    filename = "tof_filter"


class Pmsms2mzmlBinary(Executable):
    filename = "pmsms2mzml"


class PythonModule(Node):
    """The module a Python tool's launcher imports, at `source`: an input of the steps
    running the tool, so that editing it reruns them. Edits deeper in its package are not
    tracked. One subclass per tool, so a step accepts only its own tool's module."""

    source: str


TIMSTOFU = "git/timstofu/src/timstofu"
SEARCHOPS_CLI = "git/searchops/src/searchops/cli"


class FindArgmaxesModule(PythonModule):
    filename = "find_argmaxes.py"
    source = f"{TIMSTOFU}/scale_estimation/find_argmaxes.py"


class ExtractSampleTensorsModule(PythonModule):
    filename = "extract_sample_tensors.py"
    source = f"{TIMSTOFU}/scale_estimation/extract_sample_tensors.py"


class FitScaleEstimatesModule(PythonModule):
    filename = "fit_scale_estimates.py"
    source = f"{TIMSTOFU}/scale_estimation/fit_scale_estimates.py"


class CountNeighborsModule(PythonModule):
    filename = "count_neighbors.py"
    source = f"{TIMSTOFU}/candidate_selection/count_neighbors.py"


class ScoreCandidatesModule(PythonModule):
    filename = "score_candidates.py"
    source = f"{TIMSTOFU}/candidate_selection/score_candidates.py"


class AnnotateCandidatesModule(PythonModule):
    filename = "annotate.py"
    source = f"{TIMSTOFU}/candidate_postprocessing/annotate.py"


class DechargeCandidatesModule(PythonModule):
    filename = "decharge.py"
    source = f"{TIMSTOFU}/candidate_postprocessing/decharge.py"


class FilterMmappetModule(PythonModule):
    filename = "filter_mmappet.py"
    source = f"{TIMSTOFU}/cli/filter_mmappet.py"


class WriteTof2mzTableModule(PythonModule):
    filename = "write_tof2mz_table.py"
    source = f"{TIMSTOFU}/cli/write_tof2mz_table.py"


class TopCellPrecursorsModule(PythonModule):
    filename = "top_cell_precursors.py"
    source = f"{TIMSTOFU}/cli/top_cell_precursors.py"


class ScoreBasedPmsmsFilterModule(PythonModule):
    filename = "score_based_pmsms_filter.py"
    source = f"{TIMSTOFU}/cli/score_based_pmsms_filter.py"


class SelectRecalibrationPrecursorsModule(PythonModule):
    filename = "select_recalibration_precursors.py"
    source = f"{TIMSTOFU}/cli/select_recalibration_precursors.py"


class TransmitPrecursorsModule(PythonModule):
    filename = "cli.py"
    source = "git/quadops/src/quadops/fragops/cli.py"


class CutAndIndexPrecursorsModule(PythonModule):
    filename = "cut_and_index_precursors.py"
    source = "git/mgfops/mgf_ops/cli/cut_and_index_precursors.py"


class MgfWritersModule(PythonModule):
    filename = "writers.py"
    source = "git/mgfops/mgf_ops/writers.py"


class BuildPrecursorGridIndexModule(PythonModule):
    filename = "build_precursor_grid_index.py"
    source = "git/boxing/src/boxing/cli/build_precursor_grid_index.py"


class PrecursorNeighborsParamsModule(PythonModule):
    filename = "precursor_neighbors_params.py"
    source = "git/boxing/src/boxing/cli/precursor_neighbors_params.py"


class RecalibratePmsmsMzModule(PythonModule):
    filename = "recalibrate_pmsms_mz.py"
    source = f"{SEARCHOPS_CLI}/recalibrate_pmsms_mz.py"


class RecalibratePrecursorsModule(PythonModule):
    filename = "recalibrate_precursors.py"
    source = f"{SEARCHOPS_CLI}/recalibrate_precursors.py"


class PlotRecalibratedPpmModule(PythonModule):
    filename = "plot_recalibrated_ppm.py"
    source = f"{SEARCHOPS_CLI}/plot_recalibrated_ppm.py"


class UpdatePrecursorsModule(PythonModule):
    filename = "update_precursors.py"
    source = f"{SEARCHOPS_CLI}/update_precursors.py"


class FdrLevelsModule(PythonModule):
    filename = "fdr_levels.py"
    source = f"{SEARCHOPS_CLI}/fdr_levels.py"


class SageSummarizeRawModule(PythonModule):
    filename = "sage_summarize_raw.py"
    source = f"{SEARCHOPS_CLI}/sage_summarize_raw.py"


class SageFilterModule(PythonModule):
    filename = "sage_filter.py"
    source = f"{SEARCHOPS_CLI}/sage_filter.py"


class SageScoreMapperModule(PythonModule):
    filename = "sage_score_mapper.py"
    source = f"{SEARCHOPS_CLI}/sage_score_mapper.py"


class FeaturePredictionCliModule(PythonModule):
    filename = "cli.py"
    source = "git/featureprediction/src/feature_prediction/cli.py"


class ConfigSetModule(PythonModule):
    filename = "config_set.py"
    source = "git/necroflow/src/necroflow/tools/config_set.py"


class MokapotPinAdapterModule(PythonModule):
    filename = "mokapot_pin_adapter.py"
    source = "scripts/mokapot_pin_adapter.py"


class MokapotPluginModule(PythonModule):
    """A mokapot plugin's module, by its `--plugin` name in `MOKAPOT_PLUGIN_MODULES`."""


class MokapotXgboostModule(MokapotPluginModule):
    filename = "plugin.py"
    source = "git/mokapot-xgboost-plugin/src/mokapot_xgboost/plugin.py"


class MokapotLightgbmModule(MokapotPluginModule):
    filename = "plugin.py"
    source = "git/mokapot-lightgbm-plugin/src/mokapot_lightgbm/plugin.py"


MOKAPOT_PLUGIN_MODULES = {"xgboost": MokapotXgboostModule, "lightgbm": MokapotLightgbmModule}


# --- compute artifact node types ---
class PipelineConfig(Toml):
    """The job's complete effective config, as TOML."""

    filename = "pipeline_config.toml"


class Ms1Events(Node):
    filename = "events.ms1"


class Ms2Events(Node):
    filename = "events.ms2"


class Ms2TfsEvents(MmappetDataset):
    filename = "events_ms2_tfs.mmappet"


class Ms2TsfEvents(MmappetDataset):
    filename = "events_ms2_tsf.mmappet"


class ScaleEstimationConfig(Toml):
    """`[scale_estimation]`, read by the three MS1 scale-estimation steps."""

    filename = "scale_estimation_config.toml"


class ArgmaxSample(MmappetDataset):
    filename = "argmax_sample.mmappet"


class ArgmaxSieveStats(Toml):
    filename = "sieve_stats.toml"


class SampleTensors(MmappetDataset):
    filename = "sample_tensors.mmappet"


class ScaleEstimates(Node):
    filename = "scale_estimates"


class PrecursorNeighborCountConfig(Toml):
    """`[precursor_neighbor_count]`, read by `count_candidate_neighbors`."""

    filename = "precursor_neighbor_count_config.toml"


class PrecursorCandidateSelectionConfig(Toml):
    """`[precursor_candidate_selection]`, read by `score_candidates`."""

    filename = "precursor_candidate_selection_config.toml"


class RawCandidateFeatures(MmappetDataset):
    filename = "raw_candidate_features.mmappet"


class RawPrecursorClusters(MmappetDataset):
    filename = "raw_precursor_clusters.mmappet"


class PrecursorAnnotationConfig(Toml):
    """`[precursor_annotation]`, read by `annotate_precursor_clusters`."""

    filename = "precursor_annotation_config.toml"


class PostprocessingConfig(Toml):
    """`[postprocessing_of_precursors]`, read by `decharge_precursor_clusters`."""

    filename = "postprocessing_config.toml"


class AnnotatedPrecursorClusters(MmappetDataset):
    filename = "annotated_precursor_clusters.mmappet"


class PostprocessedPrecursorClusters(MmappetDataset):
    filename = "postprocessed_precursor_clusters.mmappet"


class PrecursorTransmissionConfig(Toml):
    """`[precursor_transmission]`, read by `transmit_precursors_into_fragment_space`."""

    filename = "precursor_transmission_config.toml"


class TransmittedMs1Events(MmappetDataset):
    filename = "transmitted_ms1events.mmappet"


class TransmittedPrecursorClusters(MmappetDataset):
    filename = "transmitted_precursors.mmappet"


class FirstFilterPrecursors(MmappetDataset):
    filename = "first_filter_precursors.mmappet"


class PseudomsmsConfig(Toml):
    """`[pseudomsms]`: mkpmsms' settings."""

    filename = "pseudomsms_config.toml"


class Pmsms(MmappetDataset):
    filename = "pmsms.mmappet"


class Ms2IndexedPrecursors(MmappetDataset):
    filename = "ms2indexed_precursors.mmappet"


class PreSageFilteredPrecursors(MmappetDataset):
    filename = "pre_sage_filtered_precursors.mmappet"


class PrecursorNeighborsConfig(Toml):
    """`[score_competition]` without `score_margin`: the precursor-neighbour search
    settings."""

    filename = "precursor_neighbors_config.toml"


class PrecursorGridIndex(Node):
    filename = "precursor_grid_index"


class PrecursorNeighborsCsr(Node):
    filename = "precursor_neighbors_csr"


class NeighborScore(MmappetDataset):
    filename = "neighbor_score.mmappet"


class TofFilteredPmsms(Pmsms):
    filename = "tof_filtered_pmsms.mmappet"


class TofFilteredPrecursors(PreSageFilteredPrecursors):
    filename = "tof_filtered_precursors.mmappet"


class SageConfig(Json):
    """`[sage]` as SAGE's JSON config."""

    filename = "sage_config.json"


class SageResultsJson(Json):
    filename = "results.json"


class Pin(Node):
    """A Percolator-IN TSV: SpecId, Label, ScanNr, <features...>, Peptide,
    Proteins. Format is generic (nothing Sage-specific) -- `SageResultsPin`
    and `SagepyRescorePin` both satisfy it, so `mokapot` accepts either."""


class SageResultsPin(Pin):
    filename = "results.sage.pin"


class SageResultsTsv(Tsv):
    filename = "results.sage.tsv"


class SageMatchedFragments(Tsv):
    filename = "matched_fragments.sage.tsv"


class SageFdrLevels(Tsv):
    filename = "results.sage.fdr_levels.tsv"


class DumpPeptidesBinary(Executable):
    filename = "dump_peptides"


class DumpPeptidesConfig(Json):
    """The digestion settings from `cfg.sage.database` that `dump_peptides`
    reads (fasta, enzyme, mods, mass bounds, decoys)."""

    filename = "dump_peptides_config.json"


class CompactPeptidesTable(MmappetDataset):
    """`dump_peptides --compact`'s peptide table: the rows as mmappet columns plus
    `mass_chunks.bin`, `wide_masses.bin` and `meta.json` sidecars, read in place by
    `sage --compact-peptides` (sliced search only)."""

    filename = "compact_peptides.mmappet"


class MokapotUsedPin(Pin):
    filename = "used.pin"


class MokapotPeptides(Tsv):
    filename = "mokapot.peptides.txt"


class MokapotPsms(Tsv):
    filename = "mokapot.psms.txt"


class MokapotFdrLevels(Tsv):
    filename = "mokapot.fdr_levels.tsv"


class ConfidentPsmsParquet(Node):
    filename = "confident_psms.parquet"


class ScoreComparisonPlots(Node):
    filename = "score_comparison"


class Mzml(Node):
    filename = "mzml.mzML"


class MzmlIdmap(MmappetDataset):
    filename = "idmap.mmappet"


class Mgf(Node):
    filename = "spectra.mgf"


class TofFilteredMzml(Mzml):
    pass


class TofFilteredMzmlIdmap(MzmlIdmap):
    pass


class TofFilteredMgf(Mgf):
    filename = "tof_filtered.mgf"


class TarGz(Node):
    filename = "archive.tar.gz"


class SageSummary(Tsv):
    filename = "results.sage.summary.tsv"


class RecalibrationPrecursorSelectionConfig(Toml):
    """`[recalibration_precursor_selection]`, read by `select_recalibration_precursors`."""

    filename = "recalibration_precursor_selection_config.toml"


class RecalibrationPrecursors(TofFilteredPrecursors):
    """A subset of TofFilteredPrecursors -- accepted wherever it is, e.g. run_sage's
    `precursors` input for the calibration-only pass."""

    filename = "recalibration_precursors.mmappet"


class TopCellPrecursors(RecalibrationPrecursors):
    """Recalibration precursors whose spectrum pointers address the raw `Ms2Events`
    store: each precursor's most probable (frame, scan)."""

    filename = "top_cell_precursors.mmappet"


class Tof2MzTable(MmappetDataset):
    """Dense tof -> m/z table (float32 column `mz`): `--tof2mz` of SAGE, the mzML
    and MGF writers, which all read fragment m/z as
    `table[tof]` (divided by `1 + fragment_shift_ppm·1e-6` when the precursors
    carry that column)."""

    filename = "tof2mz.mmappet"


class RecalibratedTof2MzTable(Tof2MzTable):
    """Tof2MzTable with the fragment `f_mz` correction folded in (float64 column `mz`,
    `table[t] / (1 + f_mz(table[t])·1e-6)`), produced by fit_fragment_mz.
    Accepted wherever Tof2MzTable is."""

    filename = "recalibrated_tof2mz.mmappet"


class RecalibrationConfig(Toml):
    """`[recalibration.mz]`, read by the m/z fits."""

    filename = "recalibration_config.toml"


class FragmentTolerance(Json):
    """Fragment-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    fit_fragment_mz from its own fitted model's residuals."""

    filename = "fragment_tolerance.json"


class PrecursorTolerance(Json):
    """Precursor-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    fit_precursor_mz from its own fitted model's residuals."""

    filename = "precursor_tolerance.json"


class SerializedMzModel(Node):
    """A fitted searchops.models.MzCorrectionModel, dumped via its own save()
    (inspection only -- fit_fragment_mz/fit_precursor_mz already
    apply the model themselves, nothing downstream reads this back in)."""

    filename = "mz_model"


class PrecursorColumn(MmappetDataset):
    """One column row-aligned with search_precursors, set by update_precursors."""


class FragmentShiftPpm(PrecursorColumn):
    """Each search precursor's `fragment_shift_ppm` (`bias + f_rt(rt)`, from
    fit_fragment_mz): SAGE divides its fragment m/z by `1 + fragment_shift_ppm·1e-6`."""

    filename = "fragment_shift_ppm.mmappet"


class CorrectedMz(PrecursorColumn):
    filename = "corrected_mz.mmappet"


class CorrectedRt(PrecursorColumn):
    filename = "corrected_rt.mmappet"


class CorrectedIim(PrecursorColumn):
    filename = "corrected_inv_ion_mobility.mmappet"


class RecalibratedPrecursors(TofFilteredPrecursors):
    """search_precursors with update_precursors' columns set: each a one-column
    dataset above, row-aligned and evaluated on the raw search_precursors."""

    filename = "recalibrated_precursors.mmappet"


class FragmentMzFitPlot(Png):
    filename = "fragment_mz_fit_plot.png"


class PrecursorMzFitPlot(Png):
    filename = "precursor_mz_fit_plot.png"


class RecalibratedPpmPlot(Png):
    filename = "recalibrated_ppm.png"


class PredictedRt(MmappetDataset):
    """The SPECLIB's Chronologer HI calibrated for this run: `rt` (float64 minutes),
    one row per peptide row, SAGE's `--predicted-rt`."""

    filename = "predicted_rt.mmappet"


class SpeclibPart(Node):
    """A SPECLIB part symlinked from the job's store."""


class SpeclibFragments(SpeclibPart):
    """A SPECLIB fragments part: `index.mmappet`, `arrays.mmappet`, `annotations.mmappet`."""

    filename = "speclib_fragments"


class SpeclibRt(SpeclibPart):
    """A SPECLIB retention_time part: `rt.mmappet`, one predicted HI per peptide row."""

    filename = "speclib_rt"


class SpeclibIim(SpeclibPart):
    """A SPECLIB inverse_ion_mobility part: `iim.mmappet`, SAGE's `--predicted-iim`."""

    filename = "speclib_iim"


class RtTolerance(Json):
    """rt_tol_sec tolerance JSON (`ValueTolSpline`-shaped, wrapped as
    `{"rt_tol_sec": {...}}` -- see `feature_prediction.tolerance
    .write_value_tol_spline_json`). Produced by `fit_precursor_rt` from the
    corrected residuals."""

    filename = "rt_tolerance.json"


class MobilityTolerance(Json):
    """Same as `RtTolerance`, for `mobility_tol`/IIM. Produced by
    `fit_precursor_iim`."""

    filename = "mobility_tolerance.json"


class PrecursorCorrectionRtModel(Json):
    filename = "precursor_correction_rt_model.json"


class PrecursorCorrectionIimModels(Node):
    """A directory of per-charge XGBoost models (`{charge}.json`, plus a
    pooled fallback `0.json`) -- see `precursor_correction.fit_iim_correction`."""

    filename = "precursor_correction_iim_models"


class PrecursorCorrectionRtFitPlot(Png):
    filename = "precursor_correction_rt_fit_plot.png"


class PrecursorCorrectionIimFitPlot(Png):
    filename = "precursor_correction_iim_fit_plot.png"


def tool_module(module: type[PythonModule]) -> PythonModule:
    """`module` linked from its `source`, as an input, inside a workflow."""
    return symlink(module, module.source)


# --- compute rules ---
@command
def tdf2ms1(tdf: BrukerD, tdf2ms: Tdf2msBinary, threads: Resource[int] = CORES):
    """The `.d`'s MS1 events sorted by (tof, urt, scan), with the split TOF/URT index."""
    ms1 = outputs(Ms1Events)
    shell(
        "{tdf2ms} ms1 {tdf} {ms1}"
        " --threads {threads} --paced-writeback-mib 1024 --implicit-tof-urt --overwrite"
        " && grep -qx 'uint64 x' {ms1}/tof_row_starts.mmappet/schema.txt"
        " && test -f {ms1}/tof_urt_diff_index.mmappet/schema.txt"
        " && grep -qx 'uint32 scan' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
        " && grep -qx 'uint32 intensity' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
        " && ! grep -Eq ' (tof|urt)$' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
    )
    return ms1


@command
def tdf2ms2(tdf: BrukerD, tdf2ms: Tdf2msBinary):
    """The `.d`'s MS2 events indexed by (frame, scan): the fragments mkpmsms picks from."""
    ms2 = outputs(Ms2Events)
    shell(
        "{tdf2ms} ms2 {tdf} {ms2} --overwrite"
        " && test -f {ms2}/frame_row_starts.mmappet/0.bin"
        " && test -f {ms2}/frame_scan_diff.mmappet/0.bin"
    )
    return ms2


@command
def tdf2ms2_tfs(tdf: BrukerD, tdf2ms: Tdf2msBinary, threads: Resource[int] = CORES):
    """The MS2 events in the query-oriented TOF/frame/scan layout; a benchmark target that
    nothing in the search reads."""
    ms2_tfs = outputs(Ms2TfsEvents)
    shell(
        "{tdf2ms} ms2-tfs {tdf} {ms2_tfs}"
        " --threads {threads} --overwrite"
        " && test -f {ms2_tfs}/data.mmappet/schema.txt"
        " && test -f {ms2_tfs}/packed_scans.mmappet/schema.txt"
        " && test -f {ms2_tfs}/packed_scans.mmappet/shape.txt"
        " && test -f {ms2_tfs}/tof_index.mmappet/schema.txt"
        " && test -f {ms2_tfs}/tof_index.mmappet/shape.txt"
        " && test -f {ms2_tfs}/frame_index.mmappet/schema.txt"
        " && test -f {ms2_tfs}/stats.json"
    )
    return ms2_tfs


@command
def tdf2ms2_tsf(tdf: BrukerD, tdf2ms: Tdf2msBinary, threads: Resource[int] = CORES):
    """The MS2 events in the compact TOF/scan/frame layout; a benchmark target that nothing
    in the search reads."""
    ms2_tsf = outputs(Ms2TsfEvents)
    shell(
        "{tdf2ms} ms2-tsf {tdf} {ms2_tsf}"
        " --threads {threads} --overwrite"
        " && test -f {ms2_tsf}/data.mmappet/schema.txt"
        " && test -f {ms2_tsf}/tof_scan_row_starts.mmappet/schema.txt"
        " && test -f {ms2_tsf}/tof_scan_row_starts.mmappet/shape.txt"
        " && test -f {ms2_tsf}/stats.json"
    )
    return ms2_tsf


@command
def find_ms1_argmaxes(
    ms1: Ms1Events,
    config: ScaleEstimationConfig,
    tool: FindArgmaxesModule,
    threads: Resource[int] = CORES,
):
    """The most intense MS1 events that are local intensity maxima, snapped to their peak
    centres (scale estimation, stage 1)."""
    argmaxes, stats = outputs(ArgmaxSample, ArgmaxSieveStats)
    shell("venvs/common/bin/ms1_find_argmaxes {ms1} {config} {argmaxes} {stats}")
    return argmaxes, stats


@command
def extract_ms1_sample_tensors(
    ms1: Ms1Events,
    argmaxes: ArgmaxSample,
    config: ScaleEstimationConfig,
    tool: ExtractSampleTensorsModule,
    threads: Resource[int] = CORES,
):
    """A small (tof, urt, scan) intensity tensor around each argmax (scale estimation, stage
    2)."""
    tensors = outputs(SampleTensors)
    shell(
        "venvs/common/bin/ms1_extract_sample_tensors {ms1} {argmaxes} {config} {tensors}"
    )
    return tensors


@command
def fit_ms1_scale_estimates(
    argmaxes: ArgmaxSample,
    stats: ArgmaxSieveStats,
    tensors: SampleTensors,
    config: ScaleEstimationConfig,
    tool: FitScaleEstimatesModule,
    threads: Resource[int] = CORES,
):
    """Peak width along tof, urt and scan as a function of position (`loc2scale`), fitted on
    the tensors' marginals (scale estimation, stage 3)."""
    scales = outputs(ScaleEstimates)
    shell(
        "venvs/common/bin/ms1_fit_scale_estimates {argmaxes} {stats} {tensors} {config} {scales}"
    )
    return scales


@command
def count_candidate_neighbors(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    config: PrecursorNeighborCountConfig,
    tool: CountNeighborsModule,
    threads: Resource[int] = CORES,
):
    """Each MS1 event's neighbour count and TIC within its scaled box: candidate selection's
    method-independent first pass."""
    features = outputs(RawCandidateFeatures)
    shell(
        "venvs/common/bin/ms1_count_candidate_neighbors {ms1} {scale_estimates} {config} {features}"
    )
    return features


@command
def score_candidates(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    features: RawCandidateFeatures,
    config: PrecursorCandidateSelectionConfig,
    tool: ScoreCandidatesModule,
    threads: Resource[int] = CORES,
):
    """Precursor candidates: events with enough neighbours, scored by the configured method
    and kept where they are the local maximum of that score."""
    clusters = outputs(RawPrecursorClusters)
    shell(
        "venvs/common/bin/ms1_score_candidates {ms1} {scale_estimates} {features} {config} {clusters}"
    )
    return clusters


@command
def annotate_precursor_clusters(
    tdf: BrukerD,
    ms1: Ms1Events,
    candidates: RawPrecursorClusters,
    scale_estimates: ScaleEstimates,
    config: PrecursorAnnotationConfig,
    tool: AnnotateCandidatesModule,
    threads: Resource[int] = CORES,
):
    """Candidates resolved over overlapping boxes, re-centred, and given m/z, RT, 1/K0,
    window group and peak metrics."""
    annotated = outputs(AnnotatedPrecursorClusters)
    shell(
        "venvs/common/bin/ms1_annotate_candidates {tdf} {ms1} {candidates} {scale_estimates} {config} {annotated}"
    )
    return annotated


@command
def decharge_precursor_clusters(
    tdf: BrukerD,
    ms1: Ms1Events,
    annotated: AnnotatedPrecursorClusters,
    config: PostprocessingConfig,
    tool: DechargeCandidatesModule,
    threads: Resource[int] = CORES,
):
    """The annotated precursors with charges from their isotope envelopes and a
    `precursor_idx`."""
    clusters = outputs(PostprocessedPrecursorClusters)
    shell(
        "venvs/common/bin/ms1_decharge_candidates {tdf} {ms1} {annotated} {config} {clusters}"
    )
    return clusters


@command
def transmit_precursors_into_fragment_space(
    tdf: BrukerD,
    clusters: PostprocessedPrecursorClusters,
    config: PrecursorTransmissionConfig,
    tool: TransmitPrecursorsModule,
    threads: Resource[int] = CORES,
):
    """The MS2 (frame, scan) cells the quadrupole transmits each precursor into, most
    probable first, and the precursors transmitted at all."""
    transpec, precursors = outputs(TransmittedMs1Events, TransmittedPrecursorClusters)
    shell(
        "venvs/common/bin/transmit_precursors {tdf} {clusters} {config} {transpec}"
        " --output-precursors {precursors} --verbose"
        " && test -f {transpec}/schema.txt"
    )
    return transpec, precursors


@command
def filter_first_precursors(
    precursors: TransmittedPrecursorClusters,
    tool: FilterMmappetModule,
    filter: str,
    threads: Resource[int] = CORES,
):
    """The transmitted precursors passing the job's first `filter` expression."""
    filtered = outputs(FirstFilterPrecursors)
    shell(
        "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
    )
    return filtered


@command
def run_mkpmsms_binary(
    mkpmsms: MkpmsmsBinary,
    ms2: Ms2Events,
    transprec: TransmittedMs1Events,
    filter_mm: FirstFilterPrecursors,
    config: PseudomsmsConfig,
    threads: Resource[int] = CORES,
):
    """Pseudo-MS/MS spectra: each precursor's fragments picked by mkpmsms from its
    transmitted MS2 cells, with a `dataindex` of each precursor's slice, in thread
    completion order."""
    pmsms = outputs(Pmsms)
    shell(
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
        " && venvs/common/bin/plot_ms2peakpicking_stats {pmsms}"
    )
    return pmsms


@command
def cut_and_index_precursors(
    filter_mm: FirstFilterPrecursors,
    pmsms: Pmsms,
    tool: CutAndIndexPrecursorsModule,
    threads: Resource[int] = CORES,
):
    """The precursors that got fragments, with their pmsms slice (`fragment_spectrum_start`,
    `fragment_event_cnt`), sorted by `precursor_idx`. Every later precursor table keeps this
    order, so it does not follow mkpmsms' write order."""
    precursors = outputs(Ms2IndexedPrecursors)
    shell(
        "venvs/common/bin/cut_and_index_precursors {filter_mm} {pmsms}/dataindex.mmappet {precursors}"
    )
    return precursors


@command
def filter_pre_sage_precursors(
    precursors: Ms2IndexedPrecursors,
    tool: FilterMmappetModule,
    filter: str,
    threads: Resource[int] = CORES,
):
    """The indexed precursors passing the job's pre-SAGE `filter` expression."""
    filtered = outputs(PreSageFilteredPrecursors)
    shell(
        "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
    )
    return filtered


@command
def write_tof2mz_table(tdf: BrukerD, tool: WriteTof2mzTableModule):
    """The `.d`'s TOF-to-m/z calibration as a dense table, through which SAGE and the
    exports read fragment m/z."""
    table = outputs(Tof2MzTable)
    shell("venvs/common/bin/python -m timstofu.cli.write_tof2mz_table {tdf} {table}")
    return table


@command
def top_cell_precursors(
    precursors: RecalibrationPrecursors,
    transmitted: TransmittedMs1Events,
    ms2: Ms2Events,
    tool: TopCellPrecursorsModule,
):
    """The recalibration precursors pointed at their most probable raw MS2 spectrum, for
    SAGE to search `ms2_events` in place."""
    top_cell = outputs(TopCellPrecursors)
    shell(
        "venvs/common/bin/python -m timstofu.cli.top_cell_precursors"
        " {precursors} {transmitted} {ms2} {top_cell}"
    )
    return top_cell


@command
def build_precursor_grid_index(
    precursors: PreSageFilteredPrecursors,
    tdf: BrukerD,
    config: PrecursorNeighborsConfig,
    tool: BuildPrecursorGridIndexModule,
):
    """Each precursor's (frame, scan, m/z) box, and a 3D grid index over the boxes."""
    grid = outputs(PrecursorGridIndex)
    shell(
        "venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid} --config {config}"
    )
    return grid


@command
def compute_precursor_neighbors(
    grid_index: PrecursorGridIndex,
    tdf: BrukerD,
    config: PrecursorNeighborsConfig,
    precursor_neighbors_csr: PrecursorNeighborsCsrBinary,
    tool: PrecursorNeighborsParamsModule,
    threads: Resource[int] = CORES,
):
    """Each precursor's overlapping neighbours, as CSR keyed by `precursor_idx`."""
    csr = outputs(PrecursorNeighborsCsr)
    shell(
        "{precursor_neighbors_csr}"
        " --boxes-input {grid_index}/boxes.mmappet --index-input {grid_index} --output {csr}"
        " $(venvs/common/bin/precursor-neighbors-params {config} {tdf})"
        " --n-threads {threads}"
    )
    return csr


@command
def tof_score_filter(
    pmsms: Pmsms,
    neighbors_csr: PrecursorNeighborsCsr,
    tof_filter: TofFilterBinary,
    threads: Resource[int] = CORES,
):
    """For every pmsms fragment, the highest mkpmsms score at the same TOF (±1) among its
    precursor's neighbours; NaN where no neighbour has one."""
    score = outputs(
        NeighborScore
    )  # no config arg -- confirmed vestigial in the Snakemake rule
    shell(
        "{tof_filter} --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
        " --out-path {score} --n-threads {threads}"
    )
    return score


@command
def materialize_tof_filtered_pmsms(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    neighbor_score: NeighborScore,
    tool: ScoreBasedPmsmsFilterModule,
    score_margin: int | float,
    threads: Resource[int] = CORES,
):
    """pmsms without the fragments a neighbour outscores by more than `score_margin`, and
    the precursors pointed at their new slices; precursors left without fragments are
    dropped."""
    pmsms_out, precursors_out = outputs(TofFilteredPmsms, TofFilteredPrecursors)
    shell(
        "venvs/common/bin/python -m timstofu.cli.score_based_pmsms_filter"
        " {pmsms} {precursors} {neighbor_score} {pmsms_out} {precursors_out}"
        " --threads {threads} --score-margin {score_margin}"
    )
    return pmsms_out, precursors_out


@command
def select_recalibration_precursors(
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationPrecursorSelectionConfig,
    tool: SelectRecalibrationPrecursorsModule,
    threads: Resource[int] = CORES,
):
    """The precursors the first (recalibration) search uses, chosen after sorting by
    `precursor_idx` so the choice does not depend on row order."""
    selected = outputs(RecalibrationPrecursors)
    shell(
        "venvs/common/bin/python -m timstofu.cli.select_recalibration_precursors"
        " {precursors} {config} {selected}"
    )
    return selected


@command
def fit_fragment_mz(
    sage_results_tsv: SageResultsTsv,
    matched_fragments: SageMatchedFragments,
    tof2mz: Tof2MzTable,
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationConfig,
    tool: RecalibratePmsmsMzModule,
    fdr: int | float,
):
    """Fragment m/z correction fitted on the first search's confident matches: a
    recalibrated tof2mz table, each precursor's `fragment_shift_ppm`, the fragment
    tolerance and a plot."""
    fragment_shift_ppm, recalibrated_tof2mz, tolerance, plot = outputs(
        FragmentShiftPpm, RecalibratedTof2MzTable, FragmentTolerance, FragmentMzFitPlot
    )
    shell(
        "venvs/common/bin/recalibrate-pmsms-mz {sage_results_tsv} {matched_fragments}"
        " {tof2mz} {precursors} {fragment_shift_ppm} {recalibrated_tof2mz} {tolerance} {plot}"
        " --config {config} --fdr {fdr}"
    )
    return fragment_shift_ppm, recalibrated_tof2mz, tolerance, plot


@command
def fit_precursor_mz(
    sage_results_tsv: SageResultsTsv,
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationConfig,
    tool: RecalibratePrecursorsModule,
    fdr: int | float,
    threads: Resource[int] = CORES,
):
    """Precursor m/z corrected on the first search's confident PSMs, with the precursor
    tolerance, a plot and the fitted model."""
    mz, tolerance, plot, model = outputs(
        CorrectedMz, PrecursorTolerance, PrecursorMzFitPlot, SerializedMzModel
    )
    shell(
        "venvs/common/bin/recalibrate-precursors {sage_results_tsv} {precursors}"
        " {mz} {tolerance} {plot} {model} --config {config} --fdr {fdr}"
    )
    return mz, tolerance, plot, model


@command
def plot_recalibrated_ppm(
    initial_sage_results_tsv: SageResultsTsv,
    sage_results_tsv: SageResultsTsv,
    initial_matched_fragments: SageMatchedFragments,
    matched_fragments: SageMatchedFragments,
    precursor_tolerance: PrecursorTolerance,
    fragment_tolerance: FragmentTolerance,
    tool: PlotRecalibratedPpmModule,
    fdr: int | float,
):
    """Precursor and fragment ppm errors before and after recalibration, with the
    tolerances."""
    plot = outputs(RecalibratedPpmPlot)
    shell(
        "venvs/common/bin/plot-recalibrated-ppm {initial_sage_results_tsv} {sage_results_tsv}"
        " {initial_matched_fragments} {matched_fragments} {precursor_tolerance} {fragment_tolerance}"
        " {plot} --fdr {fdr}"
    )
    return plot


@command
def dump_compact_peptides(
    fasta: Fasta,
    dump_peptides_config: DumpPeptidesConfig,
    dump_peptides_binary: DumpPeptidesBinary,
    threads: Resource[int] = CORES,
):
    """The FASTA's digest as SAGE's compact peptide table: the row space of the SPECLIBs and
    of the sliced search."""
    compact_peptides = outputs(CompactPeptidesTable)
    shell(
        "{dump_peptides_binary} -f {fasta} -c {dump_peptides_config} --compact {compact_peptides}"
    )
    return compact_peptides


def _json_text(obj) -> str:
    """A JSON file's text: sorted keys, indented, ending with a newline."""
    return json.dumps(obj, sort_keys=True, indent=2) + "\n"


@command
def check_speclib_fragments(
    compact_peptides: CompactPeptidesTable,
    fragments: SpeclibFragments,
    tool: FeaturePredictionCliModule,
    min_charge: int,
    max_charge: int,
):
    """The SPECLIB's fragments part, linked once `check-speclib` has matched it against the
    peptide table and charges; searches read this one."""
    checked = outputs(SpeclibFragments)
    shell(
        "ln -s $(realpath {fragments}) {checked}"
        " && venvs/featureprediction/bin/feature-prediction-check-speclib {compact_peptides} {checked}"
        " --part fragments --min-charge {min_charge} --max-charge {max_charge}"
    )
    return checked


@command
def check_speclib_rt(
    compact_peptides: CompactPeptidesTable,
    rt: SpeclibRt,
    tool: FeaturePredictionCliModule,
    min_charge: int,
    max_charge: int,
):
    """The SPECLIB's retention_time part, linked once `check-speclib` has matched it against
    the peptide table and charges; the RT fit reads this one."""
    checked = outputs(SpeclibRt)
    shell(
        "ln -s $(realpath {rt}) {checked}"
        " && venvs/featureprediction/bin/feature-prediction-check-speclib {compact_peptides} {checked}"
        " --part retention_time --min-charge {min_charge} --max-charge {max_charge}"
    )
    return checked


@command
def check_speclib_iim(
    compact_peptides: CompactPeptidesTable,
    iim: SpeclibIim,
    tool: FeaturePredictionCliModule,
    min_charge: int,
    max_charge: int,
):
    """The SPECLIB's inverse_ion_mobility part, linked once `check-speclib` has matched it
    against the peptide table and charges; the IIM fit and SAGE read this one."""
    checked = outputs(SpeclibIim)
    shell(
        "ln -s $(realpath {iim}) {checked}"
        " && venvs/featureprediction/bin/feature-prediction-check-speclib {compact_peptides} {checked}"
        " --part inverse_ion_mobility --min-charge {min_charge} --max-charge {max_charge}"
    )
    return checked


def _quoted(path: Path) -> str:
    """`path` shell-quoted, as literal text of a command template."""
    return shlex.quote(str(path)).replace("{", "{{").replace("}", "}}")


def _fdr_levels_record_flags(tool_command: str, *parents: Node | None) -> str:
    """`--node-dir`/`--command`/`--parent` flags of `sage-fdr-levels`/`mokapot-fdr-levels`:
    what the run database records about this node. `tool_command` is template text without
    placeholders."""
    return f" --node-dir {{workdir}} --command {shlex.quote(tool_command)}" + "".join(
        f" --parent {_quoted(parent.path)}" for parent in parents if parent is not None
    )


@command
def run_sage(
    pmsms: Pmsms | Ms2Events,
    precursors: PreSageFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
    sage_binary: SageBinary,
    tof2mz: Tof2MzTable,
    compact_peptides: CompactPeptidesTable,
    fdr_levels_tool: FdrLevelsModule,
    predicted_rt: PredictedRt | None = None,
    predicted_iim: SpeclibIim | None = None,
    speclib_fragments: SpeclibFragments | None = None,
    threads: Resource[int] = CORES,
):
    """Run Sage. The predictions are optional (mixed Node/`None` inputs) -- omitted for
    pass-1 and mode 1/2's plain search; `speclib_fragments` gives the `ms2_*` scoring
    columns, feature-only. `sage_config` must
    already carry `rt_tol_sec`+`rt_sigma_sec` and/or `mobility_tol`+`iim_sigma` for
    whichever of RT/IIM is passed here (`update_sage_config`), same
    requirement `Input::build()` enforces Rust-side. See
    plans/better_sage_filtering.md's B.4-B.6 and
    plans/rt_iim_independent_dimensions.md.

    Sage reads fragment m/z as `tof2mz[tof]` from `pmsms`' `tof` column. `pmsms`
    is either a pmsms, or the raw `Ms2Events` store, read in place as
    `<events.ms2>/data.mmappet` with `precursors` addressing its spectra (e.g.
    `top_cell_precursors`). For the recalibrated final search `tof2mz` is a
    RecalibratedTof2MzTable and `precursors` carry `fragment_shift_ppm`
    (RecalibratedPrecursors, from update_precursors), which Sage applies
    whenever the column is present. See plans/fragment_mz_corrected_tof2mz_table.md.

    `fdr_levels` holds the confident spectra, peptidoforms and ions at 1% FDR
    (`sage-fdr-levels`), which also records the run in `$ION_MAIDEN_REGRESSION_DB`.
    """
    results_json, results_pin, results_tsv, matched_fragments, fdr_levels = outputs(
        SageResultsJson, SageResultsPin, SageResultsTsv, SageMatchedFragments, SageFdrLevels
    )

    workdir = results_json.path.parent
    spectra = pmsms.path / "data.mmappet" if isinstance(pmsms, Ms2Events) else pmsms.path
    flags = (
        f" --tof2mz {_quoted(tof2mz.path)}"
        f" --compact-peptides {_quoted(compact_peptides.path)}"
    )
    if predicted_rt is not None:
        flags += f" --predicted-rt {_quoted(predicted_rt.path)}"
    if predicted_iim is not None:
        flags += f" --predicted-iim {_quoted(predicted_iim.path / 'iim.mmappet')}"
    if speclib_fragments is not None:
        fragments = speclib_fragments.path
        flags += (
            f" --predicted-fragment-intensity-index {_quoted(fragments / 'index.mmappet')}"
            f" --predicted-fragment-intensity-cache {_quoted(fragments / 'arrays.mmappet')}"
        )
    sage = (
        f"{_quoted(sage_binary.path)} -f {_quoted(fasta.path)}"
        f" --annotate-matches --write-pin --output_directory {_quoted(workdir)}"
        f" --pmsms {_quoted(spectra)} --precursors {_quoted(precursors.path)}{flags}"
        f" {_quoted(sage_config.path)}"
    )
    shell(
        f"{{sage_binary}} --version && {sage}"
        " && test -f {results_json} && test -f {results_pin}"
        " && test -f {results_tsv} && test -f {matched_fragments}"
        " && venvs/common/bin/sage-fdr-levels {results_tsv} {fdr_levels}"
        + _fdr_levels_record_flags(
            sage, pmsms, precursors, fasta, sage_config, sage_binary, tof2mz,
            compact_peptides, predicted_rt, predicted_iim, speclib_fragments,
        )
    )

    return results_json, results_pin, results_tsv, matched_fragments, fdr_levels


@command
def fit_precursor_rt(
    sage_results_tsv: SageResultsTsv,
    speclib_rt: SpeclibRt,
    precursors: PreSageFilteredPrecursors,
    tool: FeaturePredictionCliModule,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    rt_fit: str,
    fdr: int | float,
    threads: Resource[int] = CORES,
):
    """The SPECLIB's HI calibrated to this run's RT (`predicted_rt`, every peptide row),
    and the precursors' RT corrected on top of it (`rt`), both fitted on the first
    search's confident PSMs."""
    predicted_rt, rt, rt_tolerance, rt_model, plot = outputs(
        PredictedRt, CorrectedRt, RtTolerance, PrecursorCorrectionRtModel, PrecursorCorrectionRtFitPlot
    )
    shell(
        "venvs/featureprediction/bin/feature-prediction-correct-precursors-rt"
        " {sage_results_tsv} {speclib_rt} {precursors}"
        " {predicted_rt} {rt} {rt_tolerance} {rt_model} {plot}"
        " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
        " --tolerance-method {tolerance_method} --rt-fit {rt_fit} --fdr {fdr}"
    )
    return predicted_rt, rt, rt_tolerance, rt_model, plot


@command
def fit_precursor_iim(
    sage_results_tsv: SageResultsTsv,
    speclib_iim: SpeclibIim,
    precursors: PreSageFilteredPrecursors,
    tool: FeaturePredictionCliModule,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
    threads: Resource[int] = CORES,
):
    """The precursors' 1/K0 corrected against the SPECLIB's IIM on the first search's
    confident PSMs, with the mobility tolerance, the fitted models and a plot."""
    inv_ion_mobility, mobility_tolerance, iim_models, plot = outputs(
        CorrectedIim, MobilityTolerance, PrecursorCorrectionIimModels, PrecursorCorrectionIimFitPlot
    )
    shell(
        "venvs/featureprediction/bin/feature-prediction-correct-precursors-iim"
        " {sage_results_tsv} {speclib_iim} {precursors}"
        " {inv_ion_mobility} {mobility_tolerance} {iim_models} {plot}"
        " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
        " --tolerance-method {tolerance_method} --fdr {fdr}"
    )
    return inv_ion_mobility, mobility_tolerance, iim_models, plot


@command
def update_sage_config(
    sage_config: SageConfig,
    tool: ConfigSetModule,
    precursor_tolerance: PrecursorTolerance | None = None,
    fragment_tolerance: FragmentTolerance | None = None,
    rt_tolerance: RtTolerance | None = None,
    mobility_tolerance: MobilityTolerance | None = None,
):
    """Patches every given tolerance into `sage_config` with `config_set`."""
    recalibrated_sage_config = outputs(SageConfig)

    patches = [
        ("precursor_tol.ppm", precursor_tolerance, "ppm"),
        ("fragment_tol.ppm", fragment_tolerance, "ppm"),
        ("rt_tol_sec", rt_tolerance, "rt_tol_sec"),
        ("rt_sigma_sec", rt_tolerance, "rt_sigma_sec"),
        ("mobility_tol", mobility_tolerance, "mobility_tol"),
        ("iim_sigma", mobility_tolerance, "iim_sigma"),
    ]
    current = "{sage_config}"
    steps = []
    for target, source, field in patches:
        if source is None:
            continue
        patched = f"{{workdir}}/{target}.json"
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {patched}"
            f" --target {target} --source {_quoted(source.path)} --source-field {field}"
        )
        current = patched
    steps.append(f"cp {current} {{recalibrated_sage_config}}")
    shell(" && ".join(steps))

    return recalibrated_sage_config


@command
def update_precursors(
    precursors: PreSageFilteredPrecursors,
    tool: UpdatePrecursorsModule,
    fragment_shift_ppm: FragmentShiftPpm | None = None,
    mz: CorrectedMz | None = None,
    rt: CorrectedRt | None = None,
    inv_ion_mobility: CorrectedIim | None = None,
):
    """`precursors` with each given correction column (fragment shift, m/z, RT, 1/K0) set
    in."""
    recalibrated_precursors = outputs(RecalibratedPrecursors)

    columns = (fragment_shift_ppm, mz, rt, inv_ion_mobility)
    column_paths = "".join(f" {_quoted(column.path)}" for column in columns if column is not None)
    shell("venvs/common/bin/update-precursors {precursors} {recalibrated_precursors}" + column_paths)

    return recalibrated_precursors


@command
def mokapot(
    pin: Pin,
    sage_results_tsv: SageResultsTsv,
    pin_adapter: MokapotPinAdapterModule,
    fdr_levels_tool: FdrLevelsModule,
    plugin_module: MokapotPluginModule | None = None,
    train_fdr: float = 0.05,
    test_fdr: float = 0.01,
    folds: int = 3,
    train_top_rank: bool = False,
    plugin: str = "",
    rt_source: str | None = None,
    iim_source: str | None = None,
    seed: int = 1,
    xgboost_bagging_n_estimators: int = 0,
    xgboost_bagging_max_samples: int = 30_000,
    xgboost_bagging_seed: int = 0,
    threads: Resource[int] = CORES,
):
    """Builds its command in Python, not as a literal template -- lets `--plugin`
    be added conditionally without necroflow's string-template placeholders
    needing to express a conditional substring. Same reasoning for
    `--mode`/`--rt-source`/`--iim-source`: passing `rt_source`/`iim_source`
    selects `scripts/mokapot_pin_adapter.py --mode sage` (its leakage-safe
    feature registry, see `plans/mokapot_leakage_safe_pin.md`); omitting
    both leaves the adapter at its default `--mode passthrough` (original,
    unchanged, FileName-drop-only behavior).

    Calls the real `mokapot` CLI binary, not a Python-API wrapper --
    `git/mokapot` (our fork, pinned to a `v0.10.0`-based branch) patches
    `read_percolator` to use DuckDB instead of its original hand-rolled
    line-by-line parser and fixes `--seed` to actually seed `brew()`'s
    rng (previously a no-op: only the legacy `np.random.seed` global was
    touched). Both fixed at the source instead of worked around
    externally -- see `mokapot_integration.md`.

    `--max_workers`/`--xgboost_n_jobs`/`--lightgbm_n_jobs` are derived
    from the executing machine's core count, not job config -- mokapot's
    own CV folds are embarrassingly parallel (`mokapot.brew`) and each
    fold's fit is independently threadable, so both are free wall-clock
    wins with no effect on results beyond mokapot's own stochasticity
    (now seedable, see `seed` below). Computed here rather than passed as
    `mokapot()` kwargs because they're pure execution tuning, not part of
    the model; the provenance hash covers declared rule config, not the
    realized command string, so this can't make node identity
    core-count-dependent across machines.

    `xgboost_bagging_*` flags only emitted when `plugin == "xgboost"` and
    `xgboost_bagging_n_estimators > 0` -- 0 (both call sites' Python-level
    default) means the plugin builds a plain `XGBClassifier` exactly as
    before, byte-for-byte, so existing jobs that don't set these are
    unaffected. See `mokapot_integration.md` for what bagging mode is and
    why.

    `seed` is always passed (mokapot's own CLI default is 1; matching it
    here means omitting `[mokapot].seed` changes nothing) -- job-config
    seeding for every random-state site in this part of the pipeline:
    `brew()`'s own rng (mokapot core, fixed above) plus each plugin's
    `random_state` (`--xgboost_seed`/`--lightgbm_seed`, reusing the same
    integer). Verified end-to-end: two real runs with the same seed
    produced byte-identical output files.

    `--keep_decoys` also writes `mokapot.decoy.psms.txt` and
    `mokapot.decoy.peptides.txt` into the workdir, next to the target
    tables. The last step, `mokapot-fdr-levels`, counts confident spectra,
    peptidoforms and ions from them (ions take their charge from
    `sage_results_tsv`) and records the run in `$ION_MAIDEN_REGRESSION_DB`.
    """
    used_pin, peptides, psms, fdr_levels = outputs(
        MokapotUsedPin, MokapotPeptides, MokapotPsms, MokapotFdrLevels
    )

    workdir = psms.path.parent
    plugin_flag = f" --plugin {plugin}" if plugin else ""
    train_flag = " --train_top_rank" if train_top_rank else ""

    # The threads this call requests, not the machine's core
    # count: mokapot fans out to `max_workers` fold processes, each running
    # its own `--*_n_jobs`-wide plugin fit, so reading the machine here would
    # oversubscribe whatever else the scheduler is running alongside it.
    n_cores = max(1, int(threads))
    max_workers = min(folds, n_cores)
    per_fit_jobs = max(1, n_cores // max_workers)

    plugin_flags = ""
    if plugin == "xgboost":
        plugin_flags = f" --xgboost_n_jobs {per_fit_jobs} --xgboost_seed {seed}"
        if xgboost_bagging_n_estimators > 0:
            plugin_flags += (
                f" --xgboost_bagging_n_estimators {xgboost_bagging_n_estimators}"
                f" --xgboost_bagging_max_samples {xgboost_bagging_max_samples}"
                f" --xgboost_bagging_seed {xgboost_bagging_seed}"
            )
    elif plugin == "lightgbm":
        plugin_flags = f" --lightgbm_n_jobs {per_fit_jobs} --lightgbm_seed {seed}"

    adapter_flags = ""
    if rt_source is not None or iim_source is not None:
        adapter_flags = (
            f" --mode sage --rt-source {shlex.quote(rt_source or 'none')}"
            f" --iim-source {shlex.quote(iim_source or 'none')}"
        )
    mokapot = (
        "venvs/mokapot/bin/python scripts/mokapot_pin_adapter.py"
        f" -i {_quoted(pin.path)} -o {_quoted(used_pin.path)}{adapter_flags}"
        f" && venvs/mokapot/bin/mokapot {_quoted(used_pin.path)} --dest_dir {_quoted(workdir)}"
        f" --train_fdr {train_fdr} --test_fdr {test_fdr}"
        f" --folds {folds} --max_workers {max_workers} --seed {seed} --keep_decoys"
        f"{train_flag}{plugin_flag}{plugin_flags}"
    )
    shell(
        f"{mokapot} && test -f {{peptides}} && test -f {{psms}}"
        " && test -f {workdir}/mokapot.decoy.psms.txt"
        " && venvs/common/bin/mokapot-fdr-levels {psms} {peptides}"
        " {workdir}/mokapot.decoy.psms.txt {sage_results_tsv} {fdr_levels}"
        + _fdr_levels_record_flags(mokapot, pin, sage_results_tsv)
    )

    return used_pin, peptides, psms, fdr_levels


@command
def sage_summarize(
    sage_results_tsv: SageResultsTsv,
    sage_summarize_module: SageSummarizeModule,
    tool: SageSummarizeRawModule,
    fdr: int | float,
):
    """A summary of SAGE's results at `fdr`."""
    summary = outputs(SageSummary)
    shell(
        "venvs/common/bin/sage-summarize-raw {sage_results_tsv} {summary} --fdr {fdr}"
    )
    return summary


@command
def filter_sage_results(
    sage_results_tsv: SageResultsTsv, tool: SageFilterModule, fdr: int | float
):
    """SAGE's PSMs at `fdr`, as parquet."""
    confident_psms = outputs(ConfidentPsmsParquet)
    shell(
        "venvs/common/bin/sage-filter {sage_results_tsv} {confident_psms} --fdr {fdr}"
    )
    return confident_psms


@command
def score_comparison(
    pmsms: Pmsms,
    confident_psms: ConfidentPsmsParquet,
    matched_fragments: SageMatchedFragments,
    config: PseudomsmsConfig,
    tool: SageScoreMapperModule,
):
    """mkpmsms scores of the peaks SAGE matched for confident PSMs (its
    `fragment_pmsms_row`) against all other pmsms peaks."""
    plots = outputs(ScoreComparisonPlots)
    shell(
        "venvs/common/bin/sage_score_mapper {pmsms} {confident_psms} {matched_fragments}"
        " --config {config} -o {plots}"
    )
    return plots


@command
def convert_search_pmsms_to_mzml(
    pmsms: Pmsms,
    tof2mz: Tof2MzTable,
    precursors: PreSageFilteredPrecursors,
    pmsms2mzml: Pmsms2mzmlBinary,
    threads: Resource[int] = CORES,
):
    """The searched pmsms as mzML (numpress, zlib), with a spectrum id map."""
    mzml, idmap = outputs(TofFilteredMzml, TofFilteredMzmlIdmap)
    shell(
        "{pmsms2mzml} {pmsms} {precursors} {workdir}"
        " --threads {threads} --numpress --zlib-level 9 --tof2mz {tof2mz}"
        " && test -f {mzml} && test -f {idmap}/schema.txt"
    )
    return mzml, idmap


# `chattr +m`: no btrfs compression for the MGF (a no-op elsewhere). Compressing
# ~24 GiB of text in writeback throttles the writer: F9477 84 s -> 70 s without it,
# at the cost of the file's full size on disk.
@command
def convert_search_pmsms_to_mgf(
    pmsms: Pmsms,
    tof2mz: Tof2MzTable,
    precursors: PreSageFilteredPrecursors,
    tool: MgfWritersModule,
    config_path: str,
    threads: Resource[int] = CORES,
):
    """The searched pmsms as MGF, one entry per listed charge."""
    mgf = outputs(TofFilteredMgf)
    shell(
        "(chattr +m {workdir} 2>/dev/null || true)"
        " && venvs/common/bin/msms2mgf_multicharge {pmsms} {precursors} {config_path} {mgf}"
        " --threads_cnt {threads} --tof2mz {tof2mz}"
        " && test -f {mgf}"
    )
    return mgf


@command
def compress_with_pigz(path: Node, threads: Resource[int] = CORES):
    """`path` as a pigz-compressed tarball."""
    archive = outputs(TarGz)
    shell(
        "tar -cf - -C $(dirname {path}) $(basename {path}) | pigz -p {threads} > {archive}"
        " && test -f {archive}"
    )
    return archive


@workflow
def ionmaiden_pipeline(P: Pipeline, config: dict) -> None:
    """Main IonMaiden pipeline. A job with `[fragpipe]` needs
    `fragpipe_pipeline.py:ionmaiden_fragpipe_pipeline` instead."""
    if "fragpipe" in config:
        raise ValueError(
            "[fragpipe] runs only under git/ionmaidentools/src/ionmaidentools/"
            "fragpipe_pipeline.py:ionmaiden_fragpipe_pipeline -- point the job's "
            ".pipeline there"
        )
    build_ionmaiden_pipeline(P, config)


def build_ionmaiden_pipeline(P: Pipeline, config: dict) -> None:
    """The main IonMaiden pipeline's rules, called inside an active workflow
    (`ionmaiden_pipeline`, `fragpipe_pipeline.ionmaiden_fragpipe_pipeline`)."""
    cfg = DotDict.Recursive(config)
    P.pipeline_config = text_file(PipelineConfig, tomlkit.dumps(config))

    # Acquisition
    P.tdf = symlink(BrukerD, cfg.tdf_path)
    P.fasta = symlink(Fasta, cfg.fasta_path)
    P.mkpmsms_binary = symlink(MkpmsmsBinary, "git/ionmaidenmetal/build/mkpmsms")
    # Symlinked so edits to the installed-editable module invalidate sage_summarize.
    P.sage_summarize_module = symlink(
        SageSummarizeModule, "git/searchops/src/searchops/sage.py"
    )
    P.sage_binary = symlink(SageBinary, "git/sage/target/release/sage")
    P.dump_peptides_binary = symlink(
        DumpPeptidesBinary, "git/sage/target/release/dump_peptides"
    )

    # Raw Extraction
    tdf2ms = symlink(Tdf2msBinary, "git/ionmaidenmetal/build/tdf2ms")
    P.ms1_events = tdf2ms1(P.tdf, tdf2ms=tdf2ms)
    P.ms2_events = tdf2ms2(P.tdf, tdf2ms=tdf2ms)
    P.ms2_tfs_events = tdf2ms2_tfs(P.tdf, tdf2ms=tdf2ms)
    P.ms2_tsf_events = tdf2ms2_tsf(P.tdf, tdf2ms=tdf2ms)

    # MS1 Scale Calibration
    P.scale_estimation_config = text_file(
        ScaleEstimationConfig, tomlkit.dumps(cfg.scale_estimation)
    )
    P.argmaxes, P.argmax_sieve_stats = find_ms1_argmaxes(
        P.ms1_events, P.scale_estimation_config,
        tool=tool_module(FindArgmaxesModule),
    )
    P.sample_tensors = extract_ms1_sample_tensors(
        P.ms1_events, P.argmaxes, P.scale_estimation_config,
        tool=tool_module(ExtractSampleTensorsModule),
    )
    P.scale_estimates = fit_ms1_scale_estimates(
        P.argmaxes,
        P.argmax_sieve_stats,
        P.sample_tensors,
        P.scale_estimation_config,
        tool=tool_module(FitScaleEstimatesModule),
    )

    # Precursor Selection
    P.precursor_neighbor_count_config = text_file(
        PrecursorNeighborCountConfig, tomlkit.dumps(cfg.precursor_neighbor_count)
    )
    P.precursor_candidate_selection_config = text_file(
        PrecursorCandidateSelectionConfig, tomlkit.dumps(cfg.precursor_candidate_selection)
    )
    P.raw_candidate_features = count_candidate_neighbors(
        P.ms1_events,
        P.scale_estimates,
        P.precursor_neighbor_count_config,
        tool=tool_module(CountNeighborsModule),
    )
    P.raw_precursor_clusters = score_candidates(
        P.ms1_events,
        P.scale_estimates,
        P.raw_candidate_features,
        P.precursor_candidate_selection_config,
        tool=tool_module(ScoreCandidatesModule),
    )

    # Precursor Postprocessing
    P.precursor_annotation_config = text_file(
        PrecursorAnnotationConfig, tomlkit.dumps(cfg.precursor_annotation)
    )
    P.postprocessing_config = text_file(
        PostprocessingConfig, tomlkit.dumps(cfg.postprocessing_of_precursors)
    )
    P.annotated_precursor_clusters = annotate_precursor_clusters(
        P.tdf,
        P.ms1_events,
        P.raw_precursor_clusters,
        P.scale_estimates,
        P.precursor_annotation_config,
        tool=tool_module(AnnotateCandidatesModule),
    )
    P.postprocessed_precursor_clusters = decharge_precursor_clusters(
        P.tdf,
        P.ms1_events,
        P.annotated_precursor_clusters,
        P.postprocessing_config,
        tool=tool_module(DechargeCandidatesModule),
    )

    # Precursor Transmission
    P.precursor_transmission_config = text_file(
        PrecursorTransmissionConfig, tomlkit.dumps(cfg.precursor_transmission)
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = (
        transmit_precursors_into_fragment_space(
            P.tdf,
            P.postprocessed_precursor_clusters,
            P.precursor_transmission_config,
            tool=tool_module(TransmitPrecursorsModule),
        )
    )

    P.first_filter_precursors = filter_first_precursors(
        P.transmitted_precursor_clusters,
        filter=cfg.precursor_filters.mkpmsms.get("filter", ""),
        tool=tool_module(FilterMmappetModule),
    )

    # Pseudo-MS/MS Assembly
    P.pseudomsms_config = text_file(PseudomsmsConfig, tomlkit.dumps(cfg.pseudomsms))
    P.pmsms = run_mkpmsms_binary(
        P.mkpmsms_binary,
        P.ms2_events,
        P.transmitted_ms1events,
        P.first_filter_precursors,
        P.pseudomsms_config,
    )

    # Precursor Indexing
    P.ms2indexed_precursors = cut_and_index_precursors(
        P.first_filter_precursors, P.pmsms,
        tool=tool_module(CutAndIndexPrecursorsModule),
    )

    P.pre_sage_filtered_precursors = filter_pre_sage_precursors(
        P.ms2indexed_precursors,
        filter=cfg.precursor_filters.pre_sage.get("filter", ""),
        tool=tool_module(FilterMmappetModule),
    )

    # Score competition: each precursor's ToF scores against its neighbors'
    if "score_competition" in cfg:
        neighbors = {
            k: v for k, v in cfg.score_competition.items() if k != "score_margin"
        }
        P.precursor_neighbors_config = text_file(
            PrecursorNeighborsConfig, tomlkit.dumps(neighbors)
        )
        P.precursor_grid_index = build_precursor_grid_index(
            P.pre_sage_filtered_precursors,
            P.tdf,
            P.precursor_neighbors_config,
            tool=tool_module(BuildPrecursorGridIndexModule),
        )
        P.precursor_neighbors_csr = compute_precursor_neighbors(
            P.precursor_grid_index,
            P.tdf,
            P.precursor_neighbors_config,
            precursor_neighbors_csr=symlink(
                PrecursorNeighborsCsrBinary, "git/ionmaidenmetal/build/precursor_neighbors_csr"
            ),
            tool=tool_module(PrecursorNeighborsParamsModule),
        )
        P.neighbor_score = tof_score_filter(
            P.pmsms,
            P.precursor_neighbors_csr,
            tof_filter=symlink(TofFilterBinary, "git/ionmaidenmetal/build/tof_filter"),
        )
        P.search_pmsms, P.search_precursors = materialize_tof_filtered_pmsms(
            P.pmsms,
            P.pre_sage_filtered_precursors,
            P.neighbor_score,
            score_margin=cfg.score_competition.score_margin,
            tool=tool_module(ScoreBasedPmsmsFilterModule),
        )
    else:
        P.search_pmsms = P.pmsms
        P.search_precursors = P.pre_sage_filtered_precursors

    P.tof2mz_table = write_tof2mz_table(P.tdf, tool=tool_module(WriteTof2mzTableModule))
    current_tof2mz_table = P.tof2mz_table
    current_precursors = P.search_precursors

    # Search

    if "sage" in cfg:
        P.sage_config = text_file(SageConfig, _json_text(cfg.sage))
        current_sage_config = P.sage_config

        database = cfg.sage.database
        dump_peptides_database = {
            "decoy_tag": database["decoy_tag"],
            "enzyme": database["enzyme"],
            "generate_decoys": database["generate_decoys"],
            "max_variable_mods": database["max_variable_mods"],
            "peptide_max_mass": database["peptide_max_mass"],
            "peptide_min_mass": database["peptide_min_mass"],
            "static_mods": database["static_mods"],
            "variable_mods": database["variable_mods"],
        }
        P.dump_peptides_config = text_file(
            DumpPeptidesConfig, _json_text(dump_peptides_database)
        )
        P.compact_peptides = dump_compact_peptides(
            P.fasta, P.dump_peptides_config, P.dump_peptides_binary
        )

        speclib = cfg.get("speclib", {})
        min_charge, max_charge = cfg.sage.get("precursor_charge", (2, 4))
        charges = dict(min_charge=min_charge, max_charge=max_charge)

        def speclib_path(part: str) -> str:
            return f"{speclib.get('path')}/{part}/{speclib[part]}"

        speclib_fragments = (
            check_speclib_fragments(
                P.compact_peptides,
                symlink(SpeclibFragments, speclib_path("fragments")),
                **charges,
                tool=tool_module(FeaturePredictionCliModule),
            )
            if "fragments" in speclib
            else None
        )

        predictions = {}  # run_sage's predicted_rt / predicted_iim
        if "recalibration" in cfg:
            P.recalibration_precursor_selection_config = text_file(
                RecalibrationPrecursorSelectionConfig,
                tomlkit.dumps(cfg.recalibration_precursor_selection),
            )
            P.recalibration_precursors = select_recalibration_precursors(
                P.search_precursors,
                P.recalibration_precursor_selection_config,
                tool=tool_module(SelectRecalibrationPrecursorsModule),
            )
            if "recalibration_top_cell" in cfg:
                # The recalibration search reads each sampled precursor's most
                # probable raw MS2 spectrum in place; the fits below still apply to
                # the main pmsms and precursors.
                P.recalibration_top_cell_precursors = top_cell_precursors(
                    P.recalibration_precursors, P.transmitted_ms1events, P.ms2_events,
                    tool=tool_module(TopCellPrecursorsModule),
                )
                recalibration_spectra = P.ms2_events
                recalibration_search_precursors = P.recalibration_top_cell_precursors
            else:
                recalibration_spectra = P.search_pmsms
                recalibration_search_precursors = P.recalibration_precursors
            # The fits read rank 1 only, and extra ranks would change the SAGE
            # q-values they select on.
            P.recalibration_sage_config = text_file(
                SageConfig, _json_text({**cfg.sage, "report_psms": 1})
            )
            (
                P.filtered_sage_results_json,
                P.filtered_sage_results_pin,
                P.filtered_sage_results_tsv,
                P.filtered_sage_matched_fragments,
                P.filtered_sage_fdr_levels,
            ) = run_sage(
                recalibration_spectra,
                recalibration_search_precursors,
                P.fasta,
                P.recalibration_sage_config,
                P.sage_binary,
                tof2mz=P.tof2mz_table,
                compact_peptides=P.compact_peptides,
                fdr_levels_tool=tool_module(FdrLevelsModule),
            )
            corrections = {}  # update_precursors' columns
            tolerances = {}  # update_sage_config's tolerances
            if "mz" in cfg.recalibration:
                P.recalibration_config = text_file(
                    RecalibrationConfig, tomlkit.dumps(cfg.recalibration.mz)
                )
                (
                    P.fragment_shift_ppm,
                    P.recalibrated_tof2mz_table,
                    P.fragment_mz_search_tolerance,
                    P.fragment_mz_recalibration_fit_plot,
                ) = fit_fragment_mz(
                    P.filtered_sage_results_tsv,
                    P.filtered_sage_matched_fragments,
                    P.tof2mz_table,
                    P.search_precursors,
                    P.recalibration_config,
                    fdr=cfg.sage_summarize.fdr,
                    tool=tool_module(RecalibratePmsmsMzModule),
                )
                (
                    P.corrected_mz,
                    P.precursor_mz_search_tolerance,
                    P.precursor_mz_recalibration_fit_plot,
                    P.precursro_mz_recalibration_model_serialization,
                ) = fit_precursor_mz(
                    P.filtered_sage_results_tsv,
                    P.search_precursors,
                    P.recalibration_config,
                    fdr=cfg.sage_summarize.fdr,
                    tool=tool_module(RecalibratePrecursorsModule),
                )
                corrections.update(fragment_shift_ppm=P.fragment_shift_ppm, mz=P.corrected_mz)
                tolerances.update(
                    precursor_tolerance=P.precursor_mz_search_tolerance,
                    fragment_tolerance=P.fragment_mz_search_tolerance,
                )
                current_tof2mz_table = P.recalibrated_tof2mz_table

            if "rt" in cfg.recalibration:
                speclib_rt = check_speclib_rt(
                    P.compact_peptides,
                    symlink(SpeclibRt, speclib_path("retention_time")),
                    **charges,
                    tool=tool_module(FeaturePredictionCliModule),
                )
                rt_tolerance_lo, rt_tolerance_hi = cfg.recalibration.rt[
                    "tolerance_percentiles"
                ]
                rt_tolerance_method = cfg.recalibration.rt.get(
                    "tolerance_method", "theoretic"
                )
                (
                    P.predicted_rt,
                    P.corrected_rt,
                    tolerances["rt_tolerance"],
                    P.precursor_correction_rt_model,
                    P.precursor_correction_rt_fit_plot,
                ) = fit_precursor_rt(
                    P.filtered_sage_results_tsv,
                    speclib_rt,
                    P.search_precursors,
                    tolerance_lo=rt_tolerance_lo,
                    tolerance_hi=rt_tolerance_hi,
                    tolerance_method=rt_tolerance_method,
                    rt_fit=cfg.recalibration.rt.get("fit", "robust"),
                    fdr=cfg.sage_summarize.fdr,
                    tool=tool_module(FeaturePredictionCliModule),
                )
                corrections["rt"] = P.corrected_rt
                predictions["predicted_rt"] = P.predicted_rt

            if "iim" in cfg.recalibration:
                speclib_iim = check_speclib_iim(
                    P.compact_peptides,
                    symlink(SpeclibIim, speclib_path("inverse_ion_mobility")),
                    **charges,
                    tool=tool_module(FeaturePredictionCliModule),
                )
                iim_tolerance_lo, iim_tolerance_hi = cfg.recalibration.iim[
                    "tolerance_percentiles"
                ]
                iim_tolerance_method = cfg.recalibration.iim.get(
                    "tolerance_method", "theoretic"
                )
                (
                    P.corrected_inv_ion_mobility,
                    tolerances["mobility_tolerance"],
                    P.precursor_correction_iim_models,
                    P.precursor_correction_iim_fit_plot,
                ) = fit_precursor_iim(
                    P.filtered_sage_results_tsv,
                    speclib_iim,
                    P.search_precursors,
                    tolerance_lo=iim_tolerance_lo,
                    tolerance_hi=iim_tolerance_hi,
                    tolerance_method=iim_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
                    tool=tool_module(FeaturePredictionCliModule),
                )
                corrections["inv_ion_mobility"] = P.corrected_inv_ion_mobility
                predictions["predicted_iim"] = speclib_iim

            P.recalibrated_precursors = update_precursors(
                P.search_precursors, **corrections, tool=tool_module(UpdatePrecursorsModule)
            )
            P.recalibrated_sage_config = update_sage_config(
                P.sage_config, **tolerances, tool=tool_module(ConfigSetModule)
            )
            current_precursors = P.recalibrated_precursors
            current_sage_config = P.recalibrated_sage_config

        (
            P.sage_results_json,
            P.sage_results_pin,
            P.sage_results_tsv,
            P.sage_matched_fragments,
            P.sage_fdr_levels,
        ) = run_sage(
            P.search_pmsms,
            current_precursors,
            P.fasta,
            current_sage_config,
            P.sage_binary,
            tof2mz=current_tof2mz_table,
            compact_peptides=P.compact_peptides,
            speclib_fragments=speclib_fragments,
            **predictions,
            fdr_levels_tool=tool_module(FdrLevelsModule),
        )

        if "mz" in cfg.get("recalibration", {}):
            P.recalibrated_ppm_plot = plot_recalibrated_ppm(
                P.filtered_sage_results_tsv,
                P.sage_results_tsv,
                P.filtered_sage_matched_fragments,
                P.sage_matched_fragments,
                P.precursor_mz_search_tolerance,
                P.fragment_mz_search_tolerance,
                fdr=cfg.sage_summarize.fdr,
                tool=tool_module(PlotRecalibratedPpmModule),
            )

        P.confident_psms = filter_sage_results(
            P.sage_results_tsv, fdr=cfg.sage_summarize.fdr,
            tool=tool_module(SageFilterModule),
        )
        P.score_comparison = score_comparison(
            P.search_pmsms,
            P.confident_psms,
            P.sage_matched_fragments,
            P.pseudomsms_config,
            tool=tool_module(SageScoreMapperModule),
        )

        if "mokapot" in cfg:
            plugin = cfg.mokapot.get("plugin", "")
            mokapot_plugin_module = tool_module(MOKAPOT_PLUGIN_MODULES[plugin]) if plugin else None
            (
                P.mokapot_used_pin,
                P.mokapot_peptides,
                P.mokapot_psms,
                P.mokapot_fdr_levels,
            ) = mokapot(
                P.sage_results_pin,
                P.sage_results_tsv,
                rt_source="external" if "predicted_rt" in predictions else "none",
                iim_source="external" if "predicted_iim" in predictions else "none",
                **cfg.mokapot,
                pin_adapter=tool_module(MokapotPinAdapterModule),
                fdr_levels_tool=tool_module(FdrLevelsModule),
                plugin_module=mokapot_plugin_module,
            )

        # FDR Summary
        P.sage_summary = sage_summarize(
            P.sage_results_tsv, P.sage_summarize_module, fdr=cfg.sage_summarize.fdr,
            tool=tool_module(SageSummarizeRawModule),
        )

    # Exports -- current_tof2mz_table/current_precursors are the recalibrated
    # table and the mz/rt/iim-corrected precursors (carrying `fragment_shift_ppm`)
    # when recalibration ran (mode 2: mz only; mode 3: mz+RT+IIM), otherwise the
    # raw table and search_precursors (mode 1). Both must come from the same mode
    # SAGE1 actually searched against, or the exported headers and the SAGE
    # results disagree on what a peak's mz/rt/iim actually was -- see
    # plans/better_sage_filtering.md's B.6. The writers read fragment m/z as SAGE
    # does, from `tof` (plans/exports_from_tof2mz_table.md).
    P.search_mzml, P.search_mzml_idmap = convert_search_pmsms_to_mzml(
        P.search_pmsms,
        current_tof2mz_table,
        current_precursors,
        pmsms2mzml=symlink(Pmsms2mzmlBinary, "git/pmsms2mzml/pmsms2mzml"),
    )
    if mgf_config_path := cfg.get("mgf", {}).get("config_path"):
        P.search_mgf = convert_search_pmsms_to_mgf(
            P.search_pmsms,
            current_tof2mz_table,
            current_precursors,
            config_path=mgf_config_path,
            tool=tool_module(MgfWritersModule),
        )
