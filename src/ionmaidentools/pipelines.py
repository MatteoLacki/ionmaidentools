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
    CommandArgs,
    NodeType,
    Pipeline,
    command,
    output,
    symlink_file,
    text_file,
    workflow,
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


class SageFdrLevels(NodeType):
    filename = "results.sage.fdr_levels.tsv"


class DumpPeptidesBinary(NodeType):
    filename = "dump_peptides"


class DumpPeptidesConfig(NodeType):
    """The digestion settings from `cfg.sage.database` that `dump_peptides`
    reads (fasta, enzyme, mods, mass bounds, decoys)."""

    filename = "dump_peptides_config.json"


class CompactPeptidesTable(MmappetDataset):
    """`dump_peptides --compact`'s peptide table: the rows as mmappet columns plus
    `mass_chunks.bin`, `wide_masses.bin` and `meta.json` sidecars, read in place by
    `sage --compact-peptides` (sliced search only)."""

    filename = "compact_peptides.mmappet"


class MokapotUsedPin(NodeType):
    filename = "used.pin"


class MokapotPeptides(NodeType):
    filename = "mokapot.peptides.txt"


class MokapotPsms(NodeType):
    filename = "mokapot.psms.txt"


class MokapotFdrLevels(NodeType):
    filename = "mokapot.fdr_levels.tsv"


class ConfidentPsmsParquet(NodeType):
    filename = "confident_psms.parquet"


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


class SageSummary(NodeType):
    filename = "results.sage.summary.tsv"


class RecalibrationPrecursorSelectionConfig(NodeType):
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


class RecalibrationConfig(NodeType):
    filename = "recalibration_config.toml"


class FragmentTolerance(NodeType):
    """Fragment-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    fit_fragment_mz from its own fitted model's residuals."""

    filename = "fragment_tolerance.json"


class PrecursorTolerance(NodeType):
    """Precursor-residual tolerance JSON (`{"ppm": [lo, hi]}`), produced by
    fit_precursor_mz from its own fitted model's residuals."""

    filename = "precursor_tolerance.json"


class SerializedMzModel(NodeType):
    """A fitted searchops.models.MzCorrectionModel, dumped via its own save()
    (inspection only -- fit_fragment_mz/fit_precursor_mz already
    apply the model themselves, nothing downstream reads this back in)."""

    filename = "mz_model"


class FragmentShiftPpm(MmappetDataset):
    """Each search precursor's `fragment_shift_ppm` (`bias + f_rt(rt)`, from
    fit_fragment_mz): SAGE divides its fragment m/z by `1 + fragment_shift_ppm·1e-6`."""

    filename = "fragment_shift_ppm.mmappet"


class CorrectedMz(MmappetDataset):
    filename = "corrected_mz.mmappet"


class CorrectedRt(MmappetDataset):
    filename = "corrected_rt.mmappet"


class CorrectedIim(MmappetDataset):
    filename = "corrected_inv_ion_mobility.mmappet"


class RecalibratedPrecursors(TofFilteredPrecursors):
    """search_precursors with update_precursors' columns set: each a one-column
    dataset above, row-aligned and evaluated on the raw search_precursors."""

    filename = "recalibrated_precursors.mmappet"


class RecalibratedPpmPlot(Png):
    filename = "recalibrated_ppm.png"


class PredictedRt(MmappetDataset):
    """The SPECLIB's Chronologer HI calibrated for this run: `rt` (float64 minutes),
    one row per peptide row, SAGE's `--predicted-rt`."""

    filename = "predicted_rt.mmappet"


class SpeclibFragments(NodeType):
    """A SPECLIB fragments part: `index.mmappet` + `arrays.mmappet`."""

    filename = "speclib_fragments"


class SpeclibRt(NodeType):
    """A SPECLIB retention_time part: `rt.mmappet`, one predicted HI per peptide row."""

    filename = "speclib_rt"


class SpeclibIim(NodeType):
    """A SPECLIB inverse_ion_mobility part: `iim.mmappet`, SAGE's `--predicted-iim`."""

    filename = "speclib_iim"


class RtTolerance(NodeType):
    """rt_tol_sec tolerance JSON (`ValueTolSpline`-shaped, wrapped as
    `{"rt_tol_sec": {...}}` -- see `feature_prediction.tolerance
    .write_value_tol_spline_json`). Produced by both `predict_rt`
    (pre-correction) and `fit_precursor_rt` (post-correction,
    tighter, the one actually fed to SAGE1 when RT correction runs)."""

    filename = "rt_tolerance.json"


class MobilityTolerance(NodeType):
    """Same as `RtTolerance`, for `mobility_tol`/IIM. Produced by both
    `predict_iim` and `fit_precursor_iim`."""

    filename = "mobility_tolerance.json"


class RtFitPlot(Png):
    filename = "rt_fit.png"


class IimFitPlot(Png):
    filename = "iim_fit.png"


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
    " --threads {threads} --paced-writeback-mib 1024 --implicit-tof-urt --overwrite"
    " && grep -qx 'uint64 x' {ms1}/tof_row_starts.mmappet/schema.txt"
    " && test -f {ms1}/tof_urt_diff_index.mmappet/schema.txt"
    " && grep -qx 'uint32 scan' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
    " && grep -qx 'uint32 intensity' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt"
    " && ! grep -Eq ' (tof|urt)$' {ms1}/tof_urt_scan_ordered_data.mmappet/schema.txt",
    threads=CORES,
)
def tdf2ms1(tdf: BrukerD):
    ms1 = output(Ms1Events)
    return ms1


@command(
    "git/ionmaidenmetal/build/tdf2ms ms2 {tdf} {ms2} --overwrite"
    " && test -f {ms2}/frame_row_starts.mmappet/0.bin"
    " && test -f {ms2}/frame_scan_diff.mmappet/0.bin"
)
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


@command(
    "venvs/common/bin/ms1_find_argmaxes {ms1} {config} {argmaxes} {stats}",
    threads=CORES,
)
def find_ms1_argmaxes(ms1: Ms1Events, config: ScaleEstimationConfig):
    argmaxes = output(ArgmaxSample)
    stats = output(ArgmaxSieveStats)
    return argmaxes, stats


@command(
    "venvs/common/bin/ms1_extract_sample_tensors {ms1} {argmaxes} {config} {tensors}",
    threads=CORES,
)
def extract_ms1_sample_tensors(
    ms1: Ms1Events, argmaxes: ArgmaxSample, config: ScaleEstimationConfig
):
    tensors = output(SampleTensors)
    return tensors


@command(
    "venvs/common/bin/ms1_fit_scale_estimates {argmaxes} {stats} {tensors} {config} {scales}",
    threads=CORES,
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
    "venvs/common/bin/ms1_count_candidate_neighbors {ms1} {scale_estimates} {config} {features}",
    threads=CORES,
)
def count_candidate_neighbors(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    config: PrecursorNeighborCountConfig,
):
    features = output(RawCandidateFeatures)
    return features


@command(
    "venvs/common/bin/ms1_score_candidates {ms1} {scale_estimates} {features} {config} {clusters}",
    threads=CORES,
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
    "venvs/common/bin/ms1_annotate_candidates {tdf} {ms1} {candidates} {scale_estimates} {config} {annotated}",
    threads=CORES,
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
    "venvs/common/bin/ms1_decharge_candidates {tdf} {ms1} {annotated} {config} {clusters}",
    threads=CORES,
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
    " && test -f {transpec}/schema.txt",
    threads=CORES,
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
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}",
    threads=CORES,
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
    "venvs/common/bin/cut_and_index_precursors {filter_mm} {pmsms}/dataindex.mmappet {precursors}",
    threads=CORES,
)
def cut_and_index_precursors(filter_mm: FirstFilterPrecursors, pmsms: Pmsms):
    precursors = output(Ms2IndexedPrecursors)
    return precursors


@command(
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}",
    threads=CORES,
)
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, filter: str):
    filtered = output(PreSageFilteredPrecursors)
    return filtered


@command("venvs/common/bin/python -m timstofu.cli.write_tof2mz_table {tdf} {table}")
def write_tof2mz_table(tdf: BrukerD):
    table = output(Tof2MzTable)
    return table


@command(
    "venvs/common/bin/python -m timstofu.cli.top_cell_precursors"
    " {precursors} {transmitted} {ms2} {top_cell}"
)
def top_cell_precursors(
    precursors: RecalibrationPrecursors,
    transmitted: TransmittedMs1Events,
    ms2: Ms2Events,
):
    top_cell = output(TopCellPrecursors)
    return top_cell


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
    " {precursors} {config} {selected}",
    threads=CORES,
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
    "venvs/common/bin/recalibrate-pmsms-mz {sage_results_tsv} {matched_fragments}"
    " {tof2mz} {precursors} {fragment_shift_ppm} {recalibrated_tof2mz} {tolerance} {plot}"
    " --config {config} --fdr {fdr}",
)
def fit_fragment_mz(
    sage_results_tsv: SageResultsTsv,
    matched_fragments: SageMatchedFragments,
    tof2mz: Tof2MzTable,
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationConfig,
    fdr: int | float,
):
    fragment_shift_ppm = output(FragmentShiftPpm)
    recalibrated_tof2mz = output(RecalibratedTof2MzTable)
    tolerance = output(FragmentTolerance)
    plot = output(Png)
    return fragment_shift_ppm, recalibrated_tof2mz, tolerance, plot


@command(
    "venvs/common/bin/recalibrate-precursors {sage_results_tsv} {precursors}"
    " {mz} {tolerance} {plot} {model} --config {config} --fdr {fdr}",
    threads=CORES,
)
def fit_precursor_mz(
    sage_results_tsv: SageResultsTsv,
    precursors: PreSageFilteredPrecursors,
    config: RecalibrationConfig,
    fdr: int | float,
):
    mz = output(CorrectedMz)
    tolerance = output(PrecursorTolerance)
    plot = output(Png)
    model = output(SerializedMzModel)
    return mz, tolerance, plot, model


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


@text_file
def write_sage_config(text: str):
    config = output(SageConfig)
    return config


@text_file
def write_dump_peptides_config(text: str):
    dump_peptides_config = output(DumpPeptidesConfig)
    return dump_peptides_config


@command(
    "{dump_peptides_binary} -f {fasta} -c {dump_peptides_config} --compact {compact_peptides}",
    threads=CORES,
)
def dump_compact_peptides(
    fasta: Fasta,
    dump_peptides_config: DumpPeptidesConfig,
    dump_peptides_binary: DumpPeptidesBinary,
):
    compact_peptides = output(CompactPeptidesTable)
    return compact_peptides


def _json_text(obj) -> str:
    """A JSON file's text: sorted keys, indented, ending with a newline."""
    return json.dumps(obj, sort_keys=True, indent=2) + "\n"


# A SPECLIB part reflinked from the job's store into its own node, then checked against the
# peptide table and the search's charges.
_SOURCE_SPECLIB_PART = (
    "cp -r --reflink=always {store}/%(part)s/{name} {%(output)s}"
    " && venvs/featureprediction/bin/feature-prediction-check-speclib {compact_peptides} {%(output)s}"
    " --part %(part)s --min-charge {min_charge} --max-charge {max_charge}"
)


@command(_SOURCE_SPECLIB_PART % {"part": "fragments", "output": "fragments"})
def source_speclib_fragments(
    compact_peptides: CompactPeptidesTable,
    store: str,
    name: str,
    min_charge: int,
    max_charge: int,
):
    fragments = output(SpeclibFragments)
    return fragments


@command(_SOURCE_SPECLIB_PART % {"part": "retention_time", "output": "rt"})
def source_speclib_rt(
    compact_peptides: CompactPeptidesTable,
    store: str,
    name: str,
    min_charge: int,
    max_charge: int,
):
    rt = output(SpeclibRt)
    return rt


@command(_SOURCE_SPECLIB_PART % {"part": "inverse_ion_mobility", "output": "iim"})
def source_speclib_iim(
    compact_peptides: CompactPeptidesTable,
    store: str,
    name: str,
    min_charge: int,
    max_charge: int,
):
    iim = output(SpeclibIim)
    return iim


def _fdr_levels_record_flags(args: CommandArgs, tool_command: str) -> str:
    """`--node-dir`/`--command`/`--parent` flags of `sage-fdr-levels`/`mokapot-fdr-levels`:
    what the run database records about this node."""
    parents = []
    for value in args.inputs.values():
        if isinstance(value, tuple):
            parents.extend(value)
        elif value is not None:
            parents.append(value)
    return (
        f" --node-dir {shlex.quote(str(args.workdir))} --command {shlex.quote(tool_command)}"
        + "".join(f" --parent {shlex.quote(str(p))}" for p in parents)
    )


def _run_sage_command(args: CommandArgs) -> str:
    """Sliced SAGE search; each prediction input adds its flag when given."""
    sage_binary = shlex.quote(str(args.inputs.sage_binary))
    fasta = shlex.quote(str(args.inputs.fasta))
    reads_raw_ms2 = Path(args.inputs.pmsms).name == Ms2Events.filename
    if reads_raw_ms2:
        pmsms = shlex.quote(str(Path(args.inputs.pmsms) / "data.mmappet"))
    else:
        pmsms = shlex.quote(str(args.inputs.pmsms))
    precursors = shlex.quote(str(args.inputs.precursors))
    sage_config = shlex.quote(str(args.inputs.sage_config))
    workdir = shlex.quote(str(args.workdir))
    results_json = shlex.quote(str(args.outputs.results_json))
    results_pin = shlex.quote(str(args.outputs.results_pin))
    results_tsv = shlex.quote(str(args.outputs.results_tsv))
    matched_fragments = shlex.quote(str(args.outputs.matched_fragments))
    fdr_levels = shlex.quote(str(args.outputs.fdr_levels))

    flags = (
        f" --tof2mz {shlex.quote(str(args.inputs.tof2mz))}"
        f" --compact-peptides {shlex.quote(str(args.inputs.compact_peptides))}"
    )
    if args.inputs.predicted_rt is not None:
        flags += f" --predicted-rt {shlex.quote(str(args.inputs.predicted_rt))}"
    if args.inputs.predicted_iim is not None:
        iim = Path(args.inputs.predicted_iim) / "iim.mmappet"
        flags += f" --predicted-iim {shlex.quote(str(iim))}"
    if args.inputs.speclib_fragments is not None:
        fragments = Path(args.inputs.speclib_fragments)
        flags += (
            f" --predicted-fragment-intensity-index {shlex.quote(str(fragments / 'index.mmappet'))}"
            f" --predicted-fragment-intensity-cache {shlex.quote(str(fragments / 'arrays.mmappet'))}"
        )

    sage = (
        f"{sage_binary} -f {fasta}"
        f" --annotate-matches --write-pin --output_directory {workdir}"
        f" --pmsms {pmsms} --precursors {precursors}{flags} {sage_config}"
    )
    return (
        f"{sage_binary} --version && {sage}"
        f" && test -f {results_json} && test -f {results_pin}"
        f" && test -f {results_tsv} && test -f {matched_fragments}"
        f" && venvs/common/bin/sage-fdr-levels {results_tsv} {fdr_levels}"
        f"{_fdr_levels_record_flags(args, sage)}"
    )


@command(_run_sage_command, threads=CORES)
def run_sage(
    pmsms: Pmsms | Ms2Events,
    precursors: PreSageFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
    sage_binary: SageBinary,
    tof2mz: Tof2MzTable,
    compact_peptides: CompactPeptidesTable,
    predicted_rt: PredictedRt | None = None,
    predicted_iim: SpeclibIim | None = None,
    speclib_fragments: SpeclibFragments | None = None,
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
    results_json = output(SageResultsJson)
    results_pin = output(SageResultsPin)
    results_tsv = output(SageResultsTsv)
    matched_fragments = output(SageMatchedFragments)
    fdr_levels = output(SageFdrLevels)
    return results_json, results_pin, results_tsv, matched_fragments, fdr_levels


@command(
    "venvs/featureprediction/bin/feature-prediction-generate-rt"
    " {speclib_rt} {sage_results_tsv} {predicted_rt} {rt_tolerance} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}",
    threads=CORES,
)
def predict_rt(
    speclib_rt: SpeclibRt,
    sage_results_tsv: SageResultsTsv,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    """The SPECLIB's Chronologer HI calibrated to this run's RT on the first
    search's confident PSMs (found by their `peptide_row`), for every peptide row."""
    predicted_rt = output(PredictedRt)
    rt_tolerance = output(RtTolerance)
    plot = output(RtFitPlot)
    return predicted_rt, rt_tolerance, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-generate-iim"
    " {speclib_iim} {sage_results_tsv} {mobility_tolerance} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}",
    threads=CORES,
)
def predict_iim(
    speclib_iim: SpeclibIim,
    sage_results_tsv: SageResultsTsv,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    """The tolerance of the SPECLIB's IM2Deep 1/K0 against the first search's
    confident PSMs; SAGE reads that IIM from the SPECLIB as it is."""
    mobility_tolerance = output(MobilityTolerance)
    plot = output(IimFitPlot)
    return mobility_tolerance, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-correct-precursors-rt"
    " {sage_results_tsv} {predicted_rt} {precursors}"
    " {rt} {rt_tolerance} {rt_model} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}",
    threads=CORES,
)
def fit_precursor_rt(
    sage_results_tsv: SageResultsTsv,
    predicted_rt: PredictedRt,
    precursors: PreSageFilteredPrecursors,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    rt = output(CorrectedRt)
    rt_tolerance = output(RtTolerance)
    rt_model = output(PrecursorCorrectionRtModel)
    plot = output(PrecursorCorrectionRtFitPlot)
    return rt, rt_tolerance, rt_model, plot


@command(
    "venvs/featureprediction/bin/feature-prediction-correct-precursors-iim"
    " {sage_results_tsv} {speclib_iim} {precursors}"
    " {inv_ion_mobility} {mobility_tolerance} {iim_models} {plot}"
    " --tolerance-lo {tolerance_lo} --tolerance-hi {tolerance_hi}"
    " --tolerance-method {tolerance_method} --fdr {fdr}",
    threads=CORES,
)
def fit_precursor_iim(
    sage_results_tsv: SageResultsTsv,
    speclib_iim: SpeclibIim,
    precursors: PreSageFilteredPrecursors,
    tolerance_lo: int | float,
    tolerance_hi: int | float,
    tolerance_method: str,
    fdr: int | float,
):
    inv_ion_mobility = output(CorrectedIim)
    mobility_tolerance = output(MobilityTolerance)
    iim_models = output(PrecursorCorrectionIimModels)
    plot = output(PrecursorCorrectionIimFitPlot)
    return inv_ion_mobility, mobility_tolerance, iim_models, plot


def _update_sage_config_command(args: CommandArgs) -> str:
    """Patches every given tolerance into `sage_config` with `config_set`."""
    inputs = args.inputs
    patches = [
        ("precursor_tol.ppm", inputs.precursor_tolerance, "ppm"),
        ("fragment_tol.ppm", inputs.fragment_tolerance, "ppm"),
        ("rt_tol_sec", inputs.rt_tolerance, "rt_tol_sec"),
        ("rt_sigma_sec", inputs.rt_tolerance, "rt_sigma_sec"),
        ("mobility_tol", inputs.mobility_tolerance, "mobility_tol"),
        ("iim_sigma", inputs.mobility_tolerance, "iim_sigma"),
    ]
    current = shlex.quote(str(inputs.sage_config))
    steps = []
    for target, source, field in patches:
        if source is None:
            continue
        patched = shlex.quote(str(args.workdir / f"{target}.json"))
        steps.append(
            f".venv/bin/python -m necroflow.tools.config_set {current} {patched}"
            f" --target {target} --source {shlex.quote(str(source))}"
            f" --source-field {field}"
        )
        current = patched
    recalibrated_sage_config = shlex.quote(str(args.outputs.recalibrated_sage_config))
    steps.append(f"cp {current} {recalibrated_sage_config}")
    return " && ".join(steps)


@command(_update_sage_config_command)
def update_sage_config(
    sage_config: SageConfig,
    precursor_tolerance: PrecursorTolerance | None = None,
    fragment_tolerance: FragmentTolerance | None = None,
    rt_tolerance: RtTolerance | None = None,
    mobility_tolerance: MobilityTolerance | None = None,
):
    recalibrated_sage_config = output(SageConfig)
    return recalibrated_sage_config


def _update_precursors_command(args: CommandArgs) -> str:
    inputs = args.inputs
    columns = (inputs.fragment_shift_ppm, inputs.mz, inputs.rt, inputs.inv_ion_mobility)
    return shlex.join(
        [
            "venvs/common/bin/update-precursors",
            str(inputs.precursors),
            str(args.outputs.recalibrated_precursors),
            *(str(column) for column in columns if column is not None),
        ]
    )


@command(_update_precursors_command)
def update_precursors(
    precursors: PreSageFilteredPrecursors,
    fragment_shift_ppm: FragmentShiftPpm | None = None,
    mz: CorrectedMz | None = None,
    rt: CorrectedRt | None = None,
    inv_ion_mobility: CorrectedIim | None = None,
):
    recalibrated_precursors = output(RecalibratedPrecursors)
    return recalibrated_precursors


def _mokapot_command(args: CommandArgs) -> str:
    """Python command callback, not a static template -- lets `--plugin`
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
    pin = shlex.quote(str(args.inputs.pin))
    used_pin = shlex.quote(str(args.outputs.used_pin))
    peptides = shlex.quote(str(args.outputs.peptides))
    psms = shlex.quote(str(args.outputs.psms))
    fdr_levels = shlex.quote(str(args.outputs.fdr_levels))
    decoy_psms = shlex.quote(str(Path(args.workdir) / "mokapot.decoy.psms.txt"))
    sage_results_tsv = shlex.quote(str(args.inputs.sage_results_tsv))
    workdir = shlex.quote(str(args.workdir))
    plugin_flag = f" --plugin {args.config.plugin}" if args.config.plugin else ""

    # The threads this call was actually granted, not the machine's core
    # count: mokapot fans out to `max_workers` fold processes, each running
    # its own `--*_n_jobs`-wide plugin fit, so reading the machine here would
    # oversubscribe whatever else the scheduler is running alongside it.
    n_cores = max(1, int(args.constraints.threads))
    max_workers = min(args.config.folds, n_cores)
    per_fit_jobs = max(1, n_cores // max_workers)

    plugin_flags = ""
    if args.config.plugin == "xgboost":
        plugin_flags = (
            f" --xgboost_n_jobs {per_fit_jobs} --xgboost_seed {args.config.seed}"
        )
        if args.config.xgboost_bagging_n_estimators > 0:
            plugin_flags += (
                f" --xgboost_bagging_n_estimators {args.config.xgboost_bagging_n_estimators}"
                f" --xgboost_bagging_max_samples {args.config.xgboost_bagging_max_samples}"
                f" --xgboost_bagging_seed {args.config.xgboost_bagging_seed}"
            )
    elif args.config.plugin == "lightgbm":
        plugin_flags = (
            f" --lightgbm_n_jobs {per_fit_jobs} --lightgbm_seed {args.config.seed}"
        )

    adapter_flags = ""
    if args.config.rt_source is not None or args.config.iim_source is not None:
        rt_source = shlex.quote(args.config.rt_source or "none")
        iim_source = shlex.quote(args.config.iim_source or "none")
        adapter_flags = (
            f" --mode sage --rt-source {rt_source} --iim-source {iim_source}"
        )
    mokapot = (
        f"venvs/mokapot/bin/python scripts/mokapot_pin_adapter.py -i {pin} -o {used_pin}"
        f"{adapter_flags}"
        f" && venvs/mokapot/bin/mokapot {used_pin} --dest_dir {workdir}"
        f" --train_fdr {args.config.train_fdr} --test_fdr {args.config.test_fdr}"
        f" --folds {args.config.folds} --max_workers {max_workers} --seed {args.config.seed} --keep_decoys"
        f"{plugin_flag}{plugin_flags}"
    )
    return (
        f"{mokapot} && test -f {peptides} && test -f {psms} && test -f {decoy_psms}"
        f" && venvs/common/bin/mokapot-fdr-levels {psms} {peptides} {decoy_psms}"
        f" {sage_results_tsv} {fdr_levels}{_fdr_levels_record_flags(args, mokapot)}"
    )


@command(_mokapot_command, threads=CORES)
def mokapot(
    pin: Pin,
    sage_results_tsv: SageResultsTsv,
    train_fdr: float = 0.05,
    test_fdr: float = 0.01,
    folds: int = 3,
    plugin: str | None = None,
    rt_source: str | None = None,
    iim_source: str | None = None,
    seed: int = 1,
    xgboost_bagging_n_estimators: int = 0,
    xgboost_bagging_max_samples: int = 30_000,
    xgboost_bagging_seed: int = 0,
):
    used_pin = output(MokapotUsedPin)
    peptides = output(MokapotPeptides)
    psms = output(MokapotPsms)
    fdr_levels = output(MokapotFdrLevels)
    return used_pin, peptides, psms, fdr_levels


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
    "venvs/common/bin/sage_score_mapper {pmsms} {confident_psms} {matched_fragments}"
    " --config {config} -o {plots}"
)
def score_comparison(
    pmsms: Pmsms,
    confident_psms: ConfidentPsmsParquet,
    matched_fragments: SageMatchedFragments,
    config: PseudomsmsConfig,
):
    """mkpmsms scores of the peaks SAGE matched for confident PSMs (its
    `fragment_pmsms_row`) against all other pmsms peaks."""
    plots = output(ScoreComparisonPlots)
    return plots


@command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --threads {threads} --numpress --zlib-level 9 --tof2mz {tof2mz}"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_search_pmsms_to_mzml(
    pmsms: Pmsms,
    tof2mz: Tof2MzTable,
    precursors: PreSageFilteredPrecursors,
):
    mzml = output(TofFilteredMzml)
    idmap = output(TofFilteredMzmlIdmap)
    return mzml, idmap


# `chattr +m`: no btrfs compression for the MGF (a no-op elsewhere). Compressing
# ~24 GiB of text in writeback throttles the writer: F9477 84 s -> 70 s without it,
# at the cost of the file's full size on disk.
@command(
    "(chattr +m {workdir} 2>/dev/null || true)"
    " && venvs/common/bin/msms2mgf_multicharge {pmsms} {precursors} {config_path} {mgf}"
    " --threads_cnt {threads} --tof2mz {tof2mz}"
    " && test -f {mgf}",
    threads=CORES,
)
def convert_search_pmsms_to_mgf(
    pmsms: Pmsms,
    tof2mz: Tof2MzTable,
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
    P.pipeline_config = write_pipeline_config(text=tomlkit.dumps(config))

    # Acquisition
    P.tdf = source_bruker_d(path=cfg.tdf_path)
    P.fasta = source_fasta(path=cfg.fasta_path)
    P.mkpmsms_binary = source_mkpmsms_binary(path="git/ionmaidenmetal/build/mkpmsms")
    P.sage_summarize_module = source_sage_summarize_module(
        path="git/searchops/src/searchops/sage.py"
    )
    P.sage_binary = source_sage_binary(path="git/sage/target/release/sage")
    P.dump_peptides_binary = source_dump_peptides_binary(
        path="git/sage/target/release/dump_peptides"
    )

    # Raw Extraction
    P.ms1_events = tdf2ms1(P.tdf)
    P.ms2_events = tdf2ms2(P.tdf)
    P.ms2_tfs_events = tdf2ms2_tfs(P.tdf)
    P.ms2_tsf_events = tdf2ms2_tsf(P.tdf)

    # MS1 Scale Calibration
    P.scale_estimation_config = write_scale_estimation_config(
        text=tomlkit.dumps(cfg.scale_estimation)
    )
    P.argmaxes, P.argmax_sieve_stats = find_ms1_argmaxes(
        P.ms1_events, P.scale_estimation_config
    )
    P.sample_tensors = extract_ms1_sample_tensors(
        P.ms1_events, P.argmaxes, P.scale_estimation_config
    )
    P.scale_estimates = fit_ms1_scale_estimates(
        P.argmaxes,
        P.argmax_sieve_stats,
        P.sample_tensors,
        P.scale_estimation_config,
    )

    # Precursor Selection
    P.precursor_neighbor_count_config = write_precursor_neighbor_count_config(
        text=tomlkit.dumps(cfg.precursor_neighbor_count)
    )
    P.precursor_candidate_selection_config = write_precursor_candidate_selection_config(
        text=tomlkit.dumps(cfg.precursor_candidate_selection)
    )
    P.raw_candidate_features = count_candidate_neighbors(
        P.ms1_events,
        P.scale_estimates,
        P.precursor_neighbor_count_config,
    )
    P.raw_precursor_clusters = score_candidates(
        P.ms1_events,
        P.scale_estimates,
        P.raw_candidate_features,
        P.precursor_candidate_selection_config,
    )

    # Precursor Postprocessing
    P.precursor_annotation_config = write_precursor_annotation_config(
        text=tomlkit.dumps(cfg.precursor_annotation)
    )
    P.postprocessing_config = write_postprocessing_config(
        text=tomlkit.dumps(cfg.postprocessing_of_precursors)
    )
    P.annotated_precursor_clusters = annotate_precursor_clusters(
        P.tdf,
        P.ms1_events,
        P.raw_precursor_clusters,
        P.scale_estimates,
        P.precursor_annotation_config,
    )
    P.postprocessed_precursor_clusters = decharge_precursor_clusters(
        P.tdf,
        P.ms1_events,
        P.annotated_precursor_clusters,
        P.postprocessing_config,
    )

    # Precursor Transmission
    P.precursor_transmission_config = write_precursor_transmission_config(
        text=tomlkit.dumps(cfg.precursor_transmission)
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = (
        transmit_precursors_into_fragment_space(
            P.tdf,
            P.postprocessed_precursor_clusters,
            P.precursor_transmission_config,
        )
    )

    P.first_filter_precursors = filter_first_precursors(
        P.transmitted_precursor_clusters,
        filter=cfg.precursor_filters.mkpmsms.get("filter", ""),
    )

    # Pseudo-MS/MS Assembly
    P.pseudomsms_config = write_pseudomsms_config(text=tomlkit.dumps(cfg.pseudomsms))
    P.pmsms = run_mkpmsms_binary(
        P.mkpmsms_binary,
        P.ms2_events,
        P.transmitted_ms1events,
        P.first_filter_precursors,
        P.pseudomsms_config,
    )

    # Precursor Indexing
    P.ms2indexed_precursors = cut_and_index_precursors(
        P.first_filter_precursors, P.pmsms
    )

    P.pre_sage_filtered_precursors = filter_pre_sage_precursors(
        P.ms2indexed_precursors,
        filter=cfg.precursor_filters.pre_sage.get("filter", ""),
    )

    # Score competition: each precursor's ToF scores against its neighbors'
    if "score_competition" in cfg:
        neighbors = {
            k: v for k, v in cfg.score_competition.items() if k != "score_margin"
        }
        P.precursor_neighbors_config = write_precursor_neighbors_config(
            text=tomlkit.dumps(neighbors)
        )
        P.precursor_grid_index = build_precursor_grid_index(
            P.pre_sage_filtered_precursors,
            P.tdf,
            P.precursor_neighbors_config,
        )
        P.precursor_neighbors_csr = compute_precursor_neighbors(
            P.precursor_grid_index,
            P.tdf,
            P.precursor_neighbors_config,
        )
        P.neighbor_score = tof_score_filter(P.pmsms, P.precursor_neighbors_csr)
        P.search_pmsms, P.search_precursors = materialize_tof_filtered_pmsms(
            P.pmsms,
            P.pre_sage_filtered_precursors,
            P.neighbor_score,
            score_margin=cfg.score_competition.score_margin,
        )
    else:
        P.search_pmsms = P.pmsms
        P.search_precursors = P.pre_sage_filtered_precursors

    P.tof2mz_table = write_tof2mz_table(P.tdf)
    current_tof2mz_table = P.tof2mz_table
    current_precursors = P.search_precursors

    # Search

    if "sage" in cfg:
        P.sage_config = write_sage_config(text=_json_text(cfg.sage))

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
        P.dump_peptides_config = write_dump_peptides_config(
            text=_json_text(dump_peptides_database)
        )
        P.compact_peptides = dump_compact_peptides(
            P.fasta, P.dump_peptides_config, P.dump_peptides_binary
        )

        speclib = cfg.get("speclib", {})
        min_charge, max_charge = cfg.sage.get("precursor_charge", (2, 4))
        part = dict(
            store=speclib.get("path"), min_charge=min_charge, max_charge=max_charge
        )

        speclib_fragments = (
            source_speclib_fragments(
                P.compact_peptides,
                name=speclib["fragments"],
                **part,
            )
            if "fragments" in speclib
            else None
        )

        if "recalibration" in cfg:
            P.recalibration_precursor_selection_config = (
                write_recalibration_precursor_selection_config(
                    text=tomlkit.dumps(cfg.recalibration_precursor_selection)
                )
            )
            P.recalibration_precursors = select_recalibration_precursors(
                P.search_precursors,
                P.recalibration_precursor_selection_config,
            )
            if "recalibration_top_cell" in cfg:
                # The recalibration search reads each sampled precursor's most
                # probable raw MS2 spectrum in place; the fits below still apply to
                # the main pmsms and precursors.
                P.recalibration_top_cell_precursors = top_cell_precursors(
                    P.recalibration_precursors, P.transmitted_ms1events, P.ms2_events
                )
                recalibration_spectra = P.ms2_events
                recalibration_search_precursors = P.recalibration_top_cell_precursors
            else:
                recalibration_spectra = P.search_pmsms
                recalibration_search_precursors = P.recalibration_precursors
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
                P.sage_config,
                P.sage_binary,
                tof2mz=P.tof2mz_table,
                compact_peptides=P.compact_peptides,
            )
            corrections = {}  # update_precursors' columns
            tolerances = {}  # update_sage_config's tolerances
            if "mz" in cfg.recalibration:
                P.recalibration_config = write_recalibration_config(
                    text=tomlkit.dumps(cfg.recalibration.mz)
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
                )
                corrections.update(fragment_shift_ppm=P.fragment_shift_ppm, mz=P.corrected_mz)
                tolerances.update(
                    precursor_tolerance=P.precursor_mz_search_tolerance,
                    fragment_tolerance=P.fragment_mz_search_tolerance,
                )
                current_tof2mz_table = P.recalibrated_tof2mz_table

            # tolerance_percentiles/tolerance_method: required explicitly in
            # `[recalibration.rt]`/`[recalibration.iim]` (separate tables,
            # 2026-08-25) -- never inherited from `cfg.recalibration`'s own
            # (mz-scoped) `[recalibration.mz]` value; each of the three
            # dimensions (mz, rt, iim) gets its own percentiles/method. See
            # `git/featureprediction`'s `tolerance.select_tolerance`.
            if "rt" in cfg.recalibration:
                speclib_rt = source_speclib_rt(
                    P.compact_peptides, name=speclib["retention_time"], **part
                )
                rt_tolerance_lo, rt_tolerance_hi = cfg.recalibration.rt[
                    "tolerance_percentiles"
                ]
                rt_tolerance_method = cfg.recalibration.rt.get(
                    "tolerance_method", "theoretic"
                )
                P.predicted_rt, P.rt_tolerance, P.rt_fit_plot = predict_rt(
                    speclib_rt,
                    P.filtered_sage_results_tsv,
                    tolerance_lo=rt_tolerance_lo,
                    tolerance_hi=rt_tolerance_hi,
                    tolerance_method=rt_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
                )
                (
                    P.corrected_rt,
                    tolerances["rt_tolerance"],
                    P.precursor_correction_rt_model,
                    P.precursor_correction_rt_fit_plot,
                ) = fit_precursor_rt(
                    P.filtered_sage_results_tsv,
                    P.predicted_rt,
                    P.search_precursors,
                    tolerance_lo=rt_tolerance_lo,
                    tolerance_hi=rt_tolerance_hi,
                    tolerance_method=rt_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
                )
                corrections["rt"] = P.corrected_rt
            else:
                P.predicted_rt = None

            if "iim" in cfg.recalibration:
                speclib_iim = source_speclib_iim(
                    P.compact_peptides, name=speclib["inverse_ion_mobility"], **part
                )
                iim_tolerance_lo, iim_tolerance_hi = cfg.recalibration.iim[
                    "tolerance_percentiles"
                ]
                iim_tolerance_method = cfg.recalibration.iim.get(
                    "tolerance_method", "theoretic"
                )
                P.mobility_tolerance, P.iim_fit_plot = predict_iim(
                    speclib_iim,
                    P.filtered_sage_results_tsv,
                    tolerance_lo=iim_tolerance_lo,
                    tolerance_hi=iim_tolerance_hi,
                    tolerance_method=iim_tolerance_method,
                    fdr=cfg.sage_summarize.fdr,
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
                )
                corrections["inv_ion_mobility"] = P.corrected_inv_ion_mobility
            else:
                speclib_iim = None

            P.recalibrated_precursors = update_precursors(P.search_precursors, **corrections)
            current_precursors = P.recalibrated_precursors
            P.recalibrated_sage_config = update_sage_config(P.sage_config, **tolerances)
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
                P.recalibrated_sage_config,
                P.sage_binary,
                tof2mz=current_tof2mz_table,
                predicted_rt=P.predicted_rt,
                predicted_iim=speclib_iim,
                compact_peptides=P.compact_peptides,
                speclib_fragments=speclib_fragments,
            )

            if "mz" in cfg.recalibration:
                P.recalibrated_ppm_plot = plot_recalibrated_ppm(
                    P.filtered_sage_results_tsv,
                    P.sage_results_tsv,
                    P.filtered_sage_matched_fragments,
                    P.sage_matched_fragments,
                    P.precursor_mz_search_tolerance,
                    P.fragment_mz_search_tolerance,
                    fdr=cfg.sage_summarize.fdr,
                )
        else:
            (
                P.sage_results_json,
                P.sage_results_pin,
                P.sage_results_tsv,
                P.sage_matched_fragments,
                P.sage_fdr_levels,
            ) = run_sage(
                P.search_pmsms,
                P.search_precursors,
                P.fasta,
                P.sage_config,
                P.sage_binary,
                tof2mz=P.tof2mz_table,
                compact_peptides=P.compact_peptides,
                speclib_fragments=speclib_fragments,
            )

        P.confident_psms = filter_sage_results(
            P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
        )
        P.score_comparison = score_comparison(
            P.search_pmsms,
            P.confident_psms,
            P.sage_matched_fragments,
            P.pseudomsms_config,
        )

        mokapot_cfg = cfg.get("mokapot", {})
        P.mokapot_used_pin, P.mokapot_peptides, P.mokapot_psms, P.mokapot_fdr_levels = (
            mokapot(
                P.sage_results_pin,
                P.sage_results_tsv,
                folds=mokapot_cfg.get("folds", 3),
                # "" (not None) when unset -- necroflow's dependency-provenance
                # recording can't serialize a bare `None` scalar kwarg value to
                # TOML; "" is falsy in `_mokapot_command`'s own `if
                # args.config.plugin` check, so the effect (no --plugin flag,
                # mokapot's own default model) is identical.
                plugin=mokapot_cfg.get("plugin", ""),
                rt_source=(
                    "external"
                    if getattr(P, "predicted_rt", None) is not None
                    else "none"
                ),
                iim_source=(
                    "external" if "iim" in cfg.get("recalibration", {}) else "none"
                ),
                seed=mokapot_cfg.get("seed", 1),
                xgboost_bagging_n_estimators=mokapot_cfg.get(
                    "xgboost_bagging_n_estimators", 0
                ),
                xgboost_bagging_max_samples=mokapot_cfg.get(
                    "xgboost_bagging_max_samples", 30_000
                ),
                xgboost_bagging_seed=mokapot_cfg.get("xgboost_bagging_seed", 0),
            )
        )

        # FDR Summary
        P.sage_summary = sage_summarize(
            P.sage_results_tsv, P.sage_summarize_module, fdr=cfg.sage_summarize.fdr
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
    )
    if mgf_config_path := cfg.get("mgf", {}).get("config_path"):
        P.search_mgf = convert_search_pmsms_to_mgf(
            P.search_pmsms,
            current_tof2mz_table,
            current_precursors,
            config_path=mgf_config_path,
        )
