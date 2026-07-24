"""necroflow migration of necromerge2's Snakemake `short_test` chain.

First-pass scope only: tof-filtered branch, through `sage_summarize`. Plain-MGF/mzML,
`studies.smk`, and the `configs.smk` config generator are not ported.

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
Every rule invokes a real installed CLI directly (no custom Python wrapper CLIs). Sage's
own JSON config file is generated separately by `write_sage_config` straight from
`cfg.sage`, which already matches Sage's config schema 1:1. `run_sage` passes the pmsms,
tof2mz, and precursors mmappet datasets straight to Sage's `--pmsms`/`--tof2mz`/
`--precursors` flags (`software/sage/devel_fixed`, a patched fork -- upstream Sage only
takes mzML/MGF/TDF paths and reads `precursors.parquet` from a fixed-filename directory;
this fork reads all three inputs as explicit paths, precursors as `.mmappet` directly, no
staging directory or format conversion needed). Sage's four fixed-name output files
(`results.json`, `results.sage.pin`, `results.sage.tsv`, `matched_fragments.sage.tsv`)
are each their own typed necroflow output -- written straight into the rule call's shared
`{workdir}`, not wrapped in an opaque directory NodeType, since `SageResultsJson`/
`SageResultsPin`/`SageResultsTsv`/`SageMatchedFragments`'s filenames already match what
Sage itself writes there.

TODO(regression-db): the old Snakemake `sage_summarize`/`short_test` rules recorded
results into a SQLite regression DB and did an interactive baseline comparison. Neither
is ported here -- `sage_summary` is the pipeline's terminal output.

Optional m/z recalibration: when a job config has a `[recalibration]` section,
`ionmaiden_pipeline` runs Sage twice. The first (calibration) pass searches a selected subset
of `tof_filtered_precursors` (`select_recalibration_precursors`, default top-K most
intense, configurable via `[recalibration_precursor_selection]`) against the original
`tof2mz`. Its confident, top-ranked, FDR-filtered PSMs are used to fit a single
ppm-error-vs-(precursor)-m/z correction (`searchops.recalibration.fit_correction`),
applied in two separate places: `recalibrate_mz` applies it to `tof2mz` (the
ToF-bin-indexed lookup array Sage uses only for *fragment* m/z -- precursor m/z
reaches Sage as an already-materialized column and is never looked up through
`tof2mz`, so correcting `tof2mz` alone does not touch precursors, despite what an
earlier version of this docstring claimed), and `recalibrate_precursor_mz` applies
it directly to `tof_filtered_precursors`'s own `mz` column, producing
`recalibrated_precursors`. `recalibrate_mz` also derives shared `precursor_tol`/
`fragment_tol` bounds from a user-specified percentile cut of the precursor residual
error distribution (`update_sage_config`, via two chained `necroflow.tools.config_set`
calls) -- fragment_tol reuses the precursor window rather than being computed from
Sage's own `fragment_ppm`, which is an absolute-value quantity with no usable sign
(see `sage_rescoring.md`). The second, final pass re-runs Sage with the corrected
`tof2mz`, the corrected `recalibrated_precursors`, and the narrowed config. Jobs
without `[recalibration]` keep today's single-pass behaviour untouched.

Optional FragPipe comparison run: when a job config has a `[fragpipe]` section,
`ionmaiden_pipeline` additionally runs FragPipe on the same `tof_filtered_pmsms`/
`tof_filtered_precursors`/`tof2mz` Sage already searched, reproducing the old Snakemake
`search_test` side-by-side comparison. `convert_tof_filtered_to_mzml` (`git/pmsms2mzml`)
turns those into an mzML + idmap; `cfg.fragpipe.workflow_path` points at a FragPipe
`.workflow` file, symlinked in as-is (unlike Sage's JSON config, this file is long,
hand-maintained, and rarely changes -- treated as a data file, not something the
pipeline dumps or patches; the user is responsible for its `database.db-path=` line
already pointing at the right fasta). `run_fragpipe` shells out to a fixed, pre-installed
`software/fragpipe/fragpipe-24.0` install (same fixed-path convention as Sage's own
`software/sage/devel_fixed`). The old rule's log-scraping is the only result handling
ported (`extract_fragpipe_log`/`summarize_fragpipe`) -- no run_info.json, no
regression-DB recording, matching `sage_summarize`'s own current scope. Jobs without
`[fragpipe]` are unaffected.

Unlike Sage, FragPipe/Philosopher has no generate-decoys-at-search-time option -- its
database must already contain them, verified empirically against a real FragPipe 24.0
install (`philosopher database --custom fasta --nodecoys` fails outright with "Workspace
not found" without a `workspace --init` first; without decoys at all, MSFragger runs fine
but Percolator/ProteinProphet error out). `generate_fragpipe_decoy_fasta` covers this via
FragPipe's own bundled Philosopher tool (`database --custom {fasta} --prefix rev_`,
matching every bundled `.workflow`'s default `database.decoy-tag=rev_`). It's a
standalone output (`P.fragpipe_decoy_fasta`), deliberately *not* wired into
`run_fragpipe` -- keeping it a dependency would mean patching the `.workflow` file's
`database.db-path=` at run time, reopening the "plain data file, never patched" design
above. For now, point that line at `generate_fragpipe_decoy_fasta`'s output by hand.
"""

from __future__ import annotations

import json
import os
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


class SageConfig(NodeType):
    filename = "sage_config.json"


class SageResultsJson(NodeType):
    filename = "results.json"


class SageResultsPin(NodeType):
    filename = "results.sage.pin"


class SageResultsTsv(NodeType):
    filename = "results.sage.tsv"


class SageMatchedFragments(NodeType):
    filename = "matched_fragments.sage.tsv"


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


class SyntheticPmsms(MmappetDataset):
    """Fragment ions for peptides simulated from a FASTA (Koina-predicted, no
    real acquisition) -- see scripts/simulate_peptides_to_pmsms.py. mz is
    baked in directly (unlike TofFilteredPmsms's tof-index form), so
    downstream converters run without --tof2mz."""

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


class RecalibrationTolerance(NodeType):
    filename = "recalibration_tolerance.json"


class RecalibratedPrecursors(TofFilteredPrecursors):
    """tof_filtered_precursors with `mz` corrected by the same ppm-error fit
    recalibrate_mz derives for tof2mz -- accepted wherever TofFilteredPrecursors
    is, e.g. run_sage's final-pass `precursors` input."""

    filename = "recalibrated_precursors.mmappet"


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


@R.command("venvs/common/bin/ms1_find_argmaxes {ms1} {config} {argmaxes} {stats}")
def find_ms1_argmaxes(ms1: Ms1Events, config: ScaleEstimationConfig):
    return ArgmaxSample[argmaxes], ArgmaxSieveStats[stats]


@R.command(
    "venvs/common/bin/ms1_extract_sample_tensors {ms1} {argmaxes} {config} {tensors}"
)
def extract_ms1_sample_tensors(
    ms1: Ms1Events, argmaxes: ArgmaxSample, config: ScaleEstimationConfig
):
    return SampleTensors[tensors]


@R.command(
    "venvs/common/bin/ms1_fit_scale_estimates {argmaxes} {stats} {tensors} {config} {scales}"
)
def fit_ms1_scale_estimates(
    argmaxes: ArgmaxSample,
    stats: ArgmaxSieveStats,
    tensors: SampleTensors,
    config: ScaleEstimationConfig,
):
    return ScaleEstimates[scales]


R.text_file(
    "write_precursor_candidate_selection_config", PrecursorCandidateSelectionConfig
)


@R.command(
    "venvs/common/bin/ms1_select_candidates {ms1} {scale_estimates} {config} {clusters}"
)
def select_precursor_candidates(
    ms1: Ms1Events,
    scale_estimates: ScaleEstimates,
    config: PrecursorCandidateSelectionConfig,
):
    return RawPrecursorClusters[clusters]


R.text_file("write_postprocessing_config", PostprocessingConfig)


@R.command(
    "venvs/common/bin/ms1_postprocess_candidates {tdf} {ms1} {candidates} {scale_estimates} {config} {clusters}"
)
def postprocess_precursor_candidates(
    tdf: BrukerD,
    ms1: Ms1Events,
    candidates: RawPrecursorClusters,
    scale_estimates: ScaleEstimates,
    config: PostprocessingConfig,
):
    return PostprocessedPrecursorClusters[clusters]


R.text_file("write_precursor_transmission_config", PrecursorTransmissionConfig)


@R.command(
    "venvs/common/bin/transmit_precursors {tdf} {clusters} {config} {transpec}"
    " --output-precursors {precursors} --verbose"
    " && test -f {transpec}/schema.txt"
)
def transmit_precursors_into_fragment_space(
    tdf: BrukerD,
    clusters: PostprocessedPrecursorClusters,
    config: PrecursorTransmissionConfig,
):
    return TransmittedMs1Events[transpec], TransmittedPrecursorClusters[precursors]


@R.command(
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
)
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
    ms2: Ms2Events,
    transprec: TransmittedMs1Events,
    filter_mm: FirstFilterPrecursors,
    config: PseudomsmsConfig,
):
    return Pmsms[pmsms]


@R.command(
    "venvs/common/bin/cut_and_index_precursors {filter_mm} {pmsms}/dataindex.mmappet {precursors}"
)
def cut_and_index_precursors(filter_mm: FirstFilterPrecursors, pmsms: Pmsms):
    return Ms2IndexedPrecursors[precursors]


@R.command(
    "venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}"
)
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, filter: str):
    return PreSageFilteredPrecursors[filtered]


R.text_file("write_precursor_neighbors_config", PrecursorNeighborsConfig)


@R.command(
    "venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid} --config {config}"
)
def build_precursor_grid_index(
    precursors: PreSageFilteredPrecursors,
    tdf: BrukerD,
    config: PrecursorNeighborsConfig,
):
    return PrecursorGridIndex[grid]


@R.command(
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
    return PrecursorNeighborsCsr[csr]


@R.command(
    "git/ionmaidenmetal/build/tof_filter --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
    " --out-path {score} --n-threads {threads}",
    threads=CORES,
)
def tof_score_filter(pmsms: Pmsms, neighbors_csr: PrecursorNeighborsCsr):
    return NeighborScore[
        score
    ]  # no config arg -- confirmed vestigial in the Snakemake rule


@R.command(
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
    return TofFilteredPmsms[pmsms_out], TofFilteredPrecursors[precursors_out]


R.text_file(
    "write_recalibration_precursor_selection_config",
    RecalibrationPrecursorSelectionConfig,
)


@R.command(
    "venvs/common/bin/python -m timstofu.cli.select_recalibration_precursors"
    " {precursors} {config} {selected}"
)
def select_recalibration_precursors(
    precursors: TofFilteredPrecursors,
    config: RecalibrationPrecursorSelectionConfig,
):
    return RecalibrationPrecursors[selected]


R.text_file("write_recalibration_config", RecalibrationConfig)


@R.command(
    "venvs/common/bin/recalibrate-mz {sage_results_tsv} {tof2mz}"
    " {recalibrated_tof2mz} {tolerance} --config {config} --fdr {fdr}"
)
def recalibrate_mz(
    sage_results_tsv: SageResultsTsv,
    tof2mz: Tof2Mz,
    config: RecalibrationConfig,
    fdr: int | float,
):
    return Tof2Mz[recalibrated_tof2mz], RecalibrationTolerance[tolerance]


@R.command(
    "venvs/common/bin/recalibrate-precursor-mz {sage_results_tsv} {precursors}"
    " {recalibrated_precursors} --config {config} --fdr {fdr}"
)
def recalibrate_precursor_mz(
    sage_results_tsv: SageResultsTsv,
    precursors: TofFilteredPrecursors,
    config: RecalibrationConfig,
    fdr: int | float,
):
    return RecalibratedPrecursors[recalibrated_precursors]


@R.command(
    ".venv/bin/python -m necroflow.tools.config_set"
    " {sage_config} {workdir}/precursor_tol_updated.json"
    " --target precursor_tol --source {tolerance} --source-field precursor_tol"
    " && .venv/bin/python -m necroflow.tools.config_set"
    " {workdir}/precursor_tol_updated.json {recalibrated_sage_config}"
    " --target fragment_tol --source {tolerance} --source-field fragment_tol"
)
def update_sage_config(sage_config: SageConfig, tolerance: RecalibrationTolerance):
    return SageConfig[recalibrated_sage_config]


R.text_file("write_sage_config", SageConfig)


@R.command(
    "software/sage/devel_fixed/sage --version"
    " && software/sage/devel_fixed/sage -f {fasta} --annotate-matches --write-pin"
    " --output_directory {workdir} --pmsms {pmsms} --tof2mz {tof2mz} --precursors {precursors}"
    " {sage_config}"
    " && test -f {results_json} && test -f {results_pin}"
    " && test -f {results_tsv} && test -f {matched_fragments}"
)
def run_sage(
    pmsms: TofFilteredPmsms,
    tof2mz: Tof2Mz,
    precursors: TofFilteredPrecursors,
    fasta: Fasta,
    sage_config: SageConfig,
):
    return (
        SageResultsJson[results_json],
        SageResultsPin[results_pin],
        SageResultsTsv[results_tsv],
        SageMatchedFragments[matched_fragments],
    )


@R.command(
    "venvs/mokapot/bin/python scripts/mokapot_pin_adapter.py -i {sage_results_pin} -o {used_pin}"
    " && venvs/mokapot/bin/mokapot {used_pin} --dest_dir {workdir}"
    " --train_fdr 0.05 --test_fdr 0.01"
    " && test -f {peptides} && test -f {psms}"
)
def mokapot(sage_results_pin: SageResultsPin):
    return MokapotUsedPin[used_pin], MokapotPeptides[peptides], MokapotPsms[psms]


@R.command(
    "venvs/common/bin/sage-summarize-raw {sage_results_tsv} {summary} --fdr {fdr}"
)
def sage_summarize(sage_results_tsv: SageResultsTsv, fdr: int | float):
    return SageSummary[summary]


@R.command("venvs/common/bin/sage-filter {sage_results_tsv} {confident_psms} --fdr {fdr}")
def filter_sage_results(sage_results_tsv: SageResultsTsv, fdr: int | float):
    return ConfidentPsmsParquet[confident_psms]


@R.command(
    "venvs/common/bin/sage-pmsms-mapper {confident_psms} {matched_fragments}"
    " {precursors} {pmsms} {mapped} --tof2mz {tof2mz}"
)
def sage_map_to_pmsms(
    confident_psms: ConfidentPsmsParquet,
    matched_fragments: SageMatchedFragments,
    precursors: TofFilteredPrecursors,
    pmsms: TofFilteredPmsms,
    tof2mz: Tof2Mz,
):
    return SagePmsmsMapping[mapped]


@R.command(
    "venvs/common/bin/sage_score_mapper {precursors} {pmsms}"
    " {mapping}/precursors.parquet {mapping}/mapping.parquet"
    " --config {config} -o {plots}"
)
def score_comparison(
    precursors: TofFilteredPrecursors,
    pmsms: TofFilteredPmsms,
    mapping: SagePmsmsMapping,
    config: PseudomsmsConfig,
):
    return ScoreComparisonPlots[plots]


@R.command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --tof2mz {tof2mz} --threads {threads} --numpress --zlib-level 9"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_tof_filtered_to_mzml(
    pmsms: TofFilteredPmsms,
    precursors: TofFilteredPrecursors,
    tof2mz: Tof2Mz,
):
    return TofFilteredMzml[mzml], TofFilteredMzmlIdmap[idmap]


@R.command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --tof2mz {tof2mz} --threads {threads} --numpress --zlib-level 9"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_pmsms_to_mzml(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    tof2mz: Tof2Mz,
):
    return Mzml[mzml], MzmlIdmap[idmap]


@R.command(
    "venvs/common/bin/msms2mgf_multicharge {pmsms} {precursors} configs/mgf/default.toml {mgf}"
    " --tof2mz_path {tof2mz} --threads_cnt {threads}"
    " && test -f {mgf}",
    threads=CORES,
)
def convert_pmsms_to_mgf(
    pmsms: Pmsms,
    precursors: PreSageFilteredPrecursors,
    tof2mz: Tof2Mz,
):
    return Mgf[mgf]


@R.command(
    "venvs/common/bin/msms2mgf_multicharge {pmsms} {precursors} configs/mgf/default.toml {mgf}"
    " --tof2mz_path {tof2mz} --threads_cnt {threads}"
    " && test -f {mgf}",
    threads=CORES,
)
def convert_tof_filtered_to_mgf(
    pmsms: TofFilteredPmsms,
    precursors: TofFilteredPrecursors,
    tof2mz: Tof2Mz,
):
    return TofFilteredMgf[mgf]


@R.command(
    "tar -cf - -C $(dirname {path}) $(basename {path}) | pigz -p {threads} > {archive}"
    " && test -f {archive}",
    threads=CORES,
)
def compress_with_pigz(path: NodeType):
    return TarGz[archive]


@R.command("ln -s $(realpath {path}) {workflow}")
def source_fragpipe_workflow(path: str):
    return FragpipeWorkflow[workflow]


@R.command(
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
    return FragpipeDecoyFasta[decoy_fasta]


@R.command(
    "venvs/common/bin/python scripts/simulate_peptides_to_pmsms.py"
    " {fasta} {pmsms} {precursors}"
    " --charges {charges} --max-peptides-per-protein {max_peptides_per_protein}"
    " --seed {seed}"
    " && test -f {pmsms}/schema.txt && test -f {pmsms}/dataindex.mmappet/schema.txt"
    " && test -f {precursors}/schema.txt"
)
def simulate_peptides_to_pmsms(
    fasta: Fasta, charges: str, max_peptides_per_protein: int, seed: int,
):
    return SyntheticPmsms[pmsms], SyntheticPrecursors[precursors]


@R.command(
    "git/pmsms2mzml/pmsms2mzml {pmsms} {precursors} {workdir}"
    " --threads {threads} --numpress --zlib-level 9"
    " && test -f {mzml} && test -f {idmap}/schema.txt",
    threads=CORES,
)
def convert_synthetic_pmsms_to_mzml(pmsms: SyntheticPmsms, precursors: SyntheticPrecursors):
    return SyntheticMzml[mzml], SyntheticMzmlIdmap[idmap]


@R.command(
    "venvs/common/bin/msms2mgf {pmsms} {precursors} configs/mgf/default.toml {mgf}"
    " && test -f {mgf}"
)
def convert_synthetic_pmsms_to_mgf(pmsms: SyntheticPmsms, precursors: SyntheticPrecursors):
    return SyntheticMgf[mgf]


@R.command('printf "%s\\tA\\t1\\tDDA" "$(realpath {mzml})" > {manifest}')
def write_fragpipe_manifest(mzml: TofFilteredMzml):
    return FragpipeManifest[manifest]


@R.command(
    "software/fragpipe/fragpipe-24.0/bin/fragpipe --headless"
    " --workflow {workflow} --manifest {manifest} --workdir {dir}"
    " --config-tools-folder software/fragpipe/fragpipe-24.0/tools"
    " --threads {threads} --ram {ram}"
    ' && test -n "$(ls {dir}/log_*.txt 2>/dev/null)"',
    threads=CORES,
)
def run_fragpipe(manifest: FragpipeManifest, workflow: FragpipeWorkflow, ram: int = 0):
    return FragpipeResultsDir[dir]


@R.command("cp $(ls {dir}/log_*.txt | head -n 1) {log}")
def extract_fragpipe_log(dir: FragpipeResultsDir):
    return FragpipeLog[log]


@R.command(
    "grep -m 1 temp {log} > {summary}"
    " && grep -n 'MASS CALIBRATION' {log} -A 8 >> {summary}"
    " && grep 'Final report numbers after FDR filtering, and post-processing' {log} >> {summary}"
)
def summarize_fragpipe(log: FragpipeLog):
    return FragpipeSummary[summary]


def ionmaiden_pipeline(config: dict) -> Pipeline:
    """tof-filtered Sage search chain reproducing the old short_test Snakemake target."""
    cfg = DotDict.Recursive(config)
    P = Pipeline()

    P.section("Acquisition")
    P.tdf = R.source_bruker_d(path=cfg.tdf_path)
    P.fasta = R.source_fasta(path=cfg.fasta_path)

    P.section("Raw Extraction")
    P.ms1_events = R.tdf2ms1(P.tdf)
    P.ms2_events = R.tdf2ms2(P.tdf)
    P.tof2mz = R.tdf2tof2mz(P.tdf, P.ms2_events)

    P.section("MS1 Scale Calibration")
    P.scale_estimation_config = R.write_scale_estimation_config(
        text=tomlkit.dumps(cfg.scale_estimation)
    )
    P.argmaxes, P.argmax_sieve_stats = R.find_ms1_argmaxes(
        P.ms1_events, P.scale_estimation_config
    )
    P.sample_tensors = R.extract_ms1_sample_tensors(
        P.ms1_events, P.argmaxes, P.scale_estimation_config
    )
    P.scale_estimates = R.fit_ms1_scale_estimates(
        P.argmaxes,
        P.argmax_sieve_stats,
        P.sample_tensors,
        P.scale_estimation_config,
    )

    P.section("Precursor Selection")
    P.precursor_candidate_selection_config = (
        R.write_precursor_candidate_selection_config(
            text=tomlkit.dumps(cfg.precursor_candidate_selection)
        )
    )
    P.raw_precursor_clusters = R.select_precursor_candidates(
        P.ms1_events,
        P.scale_estimates,
        P.precursor_candidate_selection_config,
    )

    P.section("Precursor Postprocessing")
    P.postprocessing_config = R.write_postprocessing_config(
        text=tomlkit.dumps(cfg.postprocessing_of_precursors)
    )
    P.postprocessed_precursor_clusters = R.postprocess_precursor_candidates(
        P.tdf,
        P.ms1_events,
        P.raw_precursor_clusters,
        P.scale_estimates,
        P.postprocessing_config,
    )

    P.section("Precursor Transmission")
    P.precursor_transmission_config = R.write_precursor_transmission_config(
        text=tomlkit.dumps(cfg.precursor_transmission)
    )
    P.transmitted_ms1events, P.transmitted_precursor_clusters = (
        R.transmit_precursors_into_fragment_space(
            P.tdf,
            P.postprocessed_precursor_clusters,
            P.precursor_transmission_config,
        )
    )

    P.first_filter_precursors = R.filter_first_precursors(
        P.transmitted_precursor_clusters,
        filter=cfg.precursor_filters.mkpmsms.get("filter", ""),
    )

    P.section("Pseudo-MS/MS Assembly")
    P.pseudomsms_config = R.write_pseudomsms_config(text=tomlkit.dumps(cfg.pseudomsms))
    P.pmsms = R.run_mkpmsms_binary(
        P.ms2_events,
        P.transmitted_ms1events,
        P.first_filter_precursors,
        P.pseudomsms_config,
    )

    P.section("Precursor Indexing")
    P.ms2indexed_precursors = R.cut_and_index_precursors(
        P.first_filter_precursors, P.pmsms
    )

    P.pre_sage_filtered_precursors = R.filter_pre_sage_precursors(
        P.ms2indexed_precursors,
        filter=cfg.precursor_filters.pre_sage.get("filter", ""),
    )

    P.section("Neighbor Graph")
    P.precursor_neighbors_config = R.write_precursor_neighbors_config(
        text=tomlkit.dumps(cfg.precursor_neighbors)
    )
    P.precursor_grid_index = R.build_precursor_grid_index(
        P.pre_sage_filtered_precursors,
        P.tdf,
        P.precursor_neighbors_config,
    )
    P.precursor_neighbors_csr = R.compute_precursor_neighbors(
        P.precursor_grid_index,
        P.tdf,
        P.precursor_neighbors_config,
    )

    P.section("ToF Score Filtering")
    P.neighbor_score = R.tof_score_filter(P.pmsms, P.precursor_neighbors_csr)

    P.tof_filtered_pmsms, P.tof_filtered_precursors = R.materialize_tof_filtered_pmsms(
        P.pmsms,
        P.pre_sage_filtered_precursors,
        P.neighbor_score,
        score_margin=cfg.tof_score_filter.score_margin,
    )

    P.plain_mzml, P.plain_mzml_idmap = R.convert_pmsms_to_mzml(
        P.pmsms, P.pre_sage_filtered_precursors, P.tof2mz,
    )
    P.plain_mgf = R.convert_pmsms_to_mgf(
        P.pmsms, P.pre_sage_filtered_precursors, P.tof2mz,
    )
    P.tof_filtered_mzml, P.tof_filtered_mzml_idmap = R.convert_tof_filtered_to_mzml(
        P.tof_filtered_pmsms, P.tof_filtered_precursors, P.tof2mz,
    )
    P.tof_filtered_mgf = R.convert_tof_filtered_to_mgf(
        P.tof_filtered_pmsms, P.tof_filtered_precursors, P.tof2mz,
    )

    P.section("Search")

    if "sage" in cfg:
        P.sage_config = R.write_sage_config(
            text=json.dumps(cfg.sage, sort_keys=True, indent=2) + "\n"
        )

        if "recalibration" in cfg:
            P.recalibration_precursor_selection_config = (
                R.write_recalibration_precursor_selection_config(
                    text=tomlkit.dumps(cfg.recalibration_precursor_selection)
                )
            )
            P.recalibration_precursors = R.select_recalibration_precursors(
                P.tof_filtered_precursors,
                P.recalibration_precursor_selection_config,
            )
            (
                P.filtered_sage_results_json,
                P.filtered_sage_results_pin,
                P.filtered_sage_results_tsv,
                P.filtered_sage_matched_fragments,
            ) = R.run_sage(
                P.tof_filtered_pmsms,
                P.tof2mz,
                P.recalibration_precursors,
                P.fasta,
                P.sage_config,
            )
            P.recalibration_config = R.write_recalibration_config(
                text=tomlkit.dumps(cfg.recalibration)
            )
            P.recalibrated_tof2mz, P.recalibration_tolerance = R.recalibrate_mz(
                P.filtered_sage_results_tsv,
                P.tof2mz,
                P.recalibration_config,
                fdr=cfg.sage_summarize.fdr,
            )
            P.recalibrated_precursors = R.recalibrate_precursor_mz(
                P.filtered_sage_results_tsv,
                P.tof_filtered_precursors,
                P.recalibration_config,
                fdr=cfg.sage_summarize.fdr,
            )
            P.recalibrated_sage_config = R.update_sage_config(
                P.sage_config, P.recalibration_tolerance
            )
            (
                P.sage_results_json,
                P.sage_results_pin,
                P.sage_results_tsv,
                P.sage_matched_fragments,
            ) = R.run_sage(
                P.tof_filtered_pmsms,
                P.recalibrated_tof2mz,
                P.recalibrated_precursors,
                P.fasta,
                P.recalibrated_sage_config,
            )
            P.confident_psms = R.filter_sage_results(
                P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
            )
            P.sage_pmsms_mapping = R.sage_map_to_pmsms(
                P.confident_psms,
                P.sage_matched_fragments,
                P.tof_filtered_precursors,
                P.tof_filtered_pmsms,
                P.recalibrated_tof2mz,
            )
            P.score_comparison = R.score_comparison(
                P.tof_filtered_precursors,
                P.tof_filtered_pmsms,
                P.sage_pmsms_mapping,
                P.pseudomsms_config,
            )
        else:
            (
                P.sage_results_json,
                P.sage_results_pin,
                P.sage_results_tsv,
                P.sage_matched_fragments,
            ) = R.run_sage(
                P.tof_filtered_pmsms,
                P.tof2mz,
                P.tof_filtered_precursors,
                P.fasta,
                P.sage_config,
            )
            P.confident_psms = R.filter_sage_results(
                P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
            )
            P.sage_pmsms_mapping = R.sage_map_to_pmsms(
                P.confident_psms,
                P.sage_matched_fragments,
                P.tof_filtered_precursors,
                P.tof_filtered_pmsms,
                P.tof2mz,
            )
            P.score_comparison = R.score_comparison(
                P.tof_filtered_precursors,
                P.tof_filtered_pmsms,
                P.sage_pmsms_mapping,
                P.pseudomsms_config,
            )

        P.mokapot_used_pin, P.mokapot_peptides, P.mokapot_psms = R.mokapot(
            P.sage_results_pin
        )

        P.section("FDR Summary")
        P.sage_summary = R.sage_summarize(
            P.sage_results_tsv, fdr=cfg.sage_summarize.fdr
        )

    if "fragpipe" in cfg:
        P.fragpipe_workflow = R.source_fragpipe_workflow(
            path=cfg.fragpipe.workflow_path
        )
        # Standalone: not an input to run_fragpipe. Point database.db-path= at
        # this file's output yourself -- see FragpipeDecoyFasta's docstring.
        P.fragpipe_decoy_fasta = R.generate_fragpipe_decoy_fasta(P.fasta)
        P.fragpipe_manifest = R.write_fragpipe_manifest(P.tof_filtered_mzml)
        P.fragpipe_results_dir = R.run_fragpipe(
            P.fragpipe_manifest, P.fragpipe_workflow, ram=cfg.fragpipe.get("ram", 0)
        )
        P.fragpipe_log = R.extract_fragpipe_log(P.fragpipe_results_dir)
        P.fragpipe_summary = R.summarize_fragpipe(P.fragpipe_log)

    return P


def fragpipe_synthetic_pipeline(config: dict) -> Pipeline:
    """FragPipe smoke test on simulated data -- no Bruker .d input, no Sage.

    1. Simulate tryptic peptides from cfg.fasta_path (a small multi-protein FASTA)
       via Koina-predicted fragment ions, written as a pmsms dataset
       (`simulate_peptides_to_pmsms`).
    2. Convert that pmsms to a spectrum file: mzML (`convert_synthetic_pmsms_to_mzml`,
       via `git/pmsms2mzml`) or MGF (`convert_synthetic_pmsms_to_mgf`, via
       `venvs/common/bin/msms2mgf` + `configs/mgf/default.toml`), picked by
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
    P = Pipeline()

    P.section("Peptide Simulation")
    P.fasta = R.source_fasta(path=cfg.fasta_path)
    P.synthetic_pmsms, P.synthetic_precursors = R.simulate_peptides_to_pmsms(
        P.fasta,
        charges=cfg.get("simulation", {}).get("charges", "2,3"),
        max_peptides_per_protein=cfg.get("simulation", {}).get(
            "max_peptides_per_protein", 12
        ),
        seed=cfg.get("simulation", {}).get("seed", 20260715),
    )

    P.section("Spectrum File Creation")
    output_format = cfg.get("output_format", "mzml")
    if output_format == "mzml":
        P.synthetic_mzml, P.synthetic_mzml_idmap = R.convert_synthetic_pmsms_to_mzml(
            P.synthetic_pmsms, P.synthetic_precursors,
        )

        P.section("FragPipe Search")
        P.fragpipe_workflow = R.source_fragpipe_workflow(path=cfg.fragpipe.workflow_path)
        P.fragpipe_decoy_fasta = R.generate_fragpipe_decoy_fasta(P.fasta)
        P.fragpipe_manifest = R.write_fragpipe_manifest(P.synthetic_mzml)
        P.fragpipe_results_dir = R.run_fragpipe(P.fragpipe_manifest, P.fragpipe_workflow)
        P.fragpipe_log = R.extract_fragpipe_log(P.fragpipe_results_dir)
        P.fragpipe_summary = R.summarize_fragpipe(P.fragpipe_log)
    elif output_format == "mgf":
        P.synthetic_mgf = R.convert_synthetic_pmsms_to_mgf(
            P.synthetic_pmsms, P.synthetic_precursors,
        )
    else:
        raise ValueError(f"cfg.output_format must be 'mzml' or 'mgf', got {output_format!r}")

    return P
