"""necroflow migration of necromerge2's Snakemake `short_test` chain.

First-pass scope only: tof-filtered branch, through `sage_summarize`. Plain-MGF/mzML,
FragPipe, `studies.smk`, and the `configs.smk` config generator are not ported.

Config fields travel as explicit named kwargs, one per fixed config field, straight from
the job TOML's own parsed dict -- no config file, no blob. The external tools this
pipeline shells out to (timstofu/quadops/boxing) were patched to accept these fields as
flags directly, with their old config-file argument made optional (and still
positionally compatible with the old Snakemake pipeline). The two tools whose CLI we own
outright (necromerge2-run-mkpmsms, necromerge2-run-sage) already took flags with no file
at all, except for Sage's own JSON file, which necromerge2-run-sage still builds
internally since the real `sage` binary requires one.

TODO(regression-db): the old Snakemake `sage_summarize`/`short_test` rules recorded
results into a SQLite regression DB and did an interactive baseline comparison. Neither
is ported here -- `sage_run_info`/`sage_summary` are the pipeline's terminal outputs.
"""
from __future__ import annotations

import json
import shlex

from necroflow import NodeType, Pipeline, Rules

R = Rules()


def _plain(x):
    """Recursively convert tomlkit's dict-like Table/Container/Array wrappers into
    plain dict/list/scalar values, so json.dumps can handle them directly."""
    if hasattr(x, "items"):
        return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_plain(v) for v in x]
    return x


def _json_kwarg(value) -> str:
    """Canonicalize+encode a genuinely open-ended (variable-key-set) config sub-dict as
    a shell-safe JSON string kwarg. Used only for fields with no fixed schema (shape
    depends on which algorithm/method is selected, or an arbitrary user-defined key set).
    Quoting happens here (in the factory), since necroflow does not shell-quote plain
    string kwargs on its own (only Path-typed substitutions get shlex.quote -- confirmed
    in necroflow/dag.py's resolve_command)."""
    return shlex.quote(json.dumps(_plain(value), sort_keys=True))


# --- source node types ---
class BrukerD(NodeType):
    filename = "input.d"


class Fasta(NodeType):
    filename = "fasta.fasta"


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


@R.command(
    "venvs/common/bin/ms1_estimate_scales {ms1} {scales} --dataset_name {dataset}"
    " --argmax-candidates-cnt {argmax_candidates_cnt}"
    " --max-iter-smoothed-argmax-search {max_iter_smoothed_argmax_search}"
    " --times-top-intensity-higher {times_top_intensity_higher}"
    " --hill-descent-buffer {hill_descent_buffer}"
    " --radii-tof {radii_tof} --radii-urt {radii_urt} --radii-scan {radii_scan}"
    " --smoothing-radii-tof {smoothing_radii_tof} --smoothing-radii-urt {smoothing_radii_urt}"
    " --smoothing-radii-scan {smoothing_radii_scan}"
    " --hq-pearson-tof {hq_pearson_tof} --hq-pearson-urt {hq_pearson_urt} --hq-pearson-scan {hq_pearson_scan}"
)
def estimate_precursor_scales(
    ms1: Ms1Events, dataset: str,
    argmax_candidates_cnt: int, max_iter_smoothed_argmax_search: int,
    times_top_intensity_higher: float, hill_descent_buffer: int,
    radii_tof: int, radii_urt: int, radii_scan: int,
    smoothing_radii_tof: int, smoothing_radii_urt: int, smoothing_radii_scan: int,
    hq_pearson_tof: float, hq_pearson_urt: float, hq_pearson_scan: float,
):
    return ScaleEstimates[scales]


@R.command(
    "venvs/common/bin/ms1_select_candidates {ms1} {scale_estimates} {clusters} --dataset_name {dataset}"
    " --tof-radius {tof_radius} --minimal-number-of-events-in-cluster {minimal_number_of_events_in_cluster}"
    " --method {method} --method-settings-json {method_settings_json}"
)
def select_precursor_candidates(
    ms1: Ms1Events, scale_estimates: ScaleEstimates, dataset: str,
    tof_radius: int, minimal_number_of_events_in_cluster: int, method: str, method_settings_json: str,
):
    return RawPrecursorClusters[clusters]


@R.command(
    "venvs/common/bin/ms1_postprocess_candidates {tdf} {ms1} {candidates} {scale_estimates} {clusters}"
    " --dataset_name {dataset}"
    " --tof-radius {tof_radius}"
    " --precursor-connected-component-policy {precursor_connected_component_policy}"
    " --location-reestimator {location_reestimator}"
    " --estimate-frame-location-method {estimate_frame_location_method}"
    " --sigmas-multiplier {sigmas_multiplier} --overlap-geometry {overlap_geometry}"
    " --decharging-method {decharging_method}"
    " --decharging-precursor-no-charge-policy {decharging_precursor_no_charge_policy}"
    " --decharging-mass-diff {decharging_mass_diff}"
    " --decharging-charges-to-look-for {decharging_charges_to_look_for}"
    " --decharging-max-jumps-right {decharging_max_jumps_right}"
    " --decharging-min-ratio-of-previous-to-next-isotope-intensity {decharging_min_ratio_of_previous_to_next_isotope_intensity}"
    " --decharging-force-same-tof-diffs {decharging_force_same_tof_diffs}"
    " --decharging-verbose {decharging_verbose}"
    " --decharging-sigmas-multiplier {decharging_sigmas_multiplier}"
)
def postprocess_precursor_candidates(
    tdf: BrukerD, ms1: Ms1Events, candidates: RawPrecursorClusters, scale_estimates: ScaleEstimates, dataset: str,
    tof_radius: int, precursor_connected_component_policy: str, location_reestimator: str,
    estimate_frame_location_method: str, sigmas_multiplier: float, overlap_geometry: str,
    decharging_method: str, decharging_precursor_no_charge_policy: str,
    decharging_mass_diff: float, decharging_charges_to_look_for: str,
    decharging_max_jumps_right: int, decharging_min_ratio_of_previous_to_next_isotope_intensity: float,
    decharging_force_same_tof_diffs: bool, decharging_verbose: bool, decharging_sigmas_multiplier: float,
):
    return PostprocessedPrecursorClusters[clusters]


@R.command(
    "venvs/common/bin/transmit_precursors {tdf} {clusters} {transpec}"
    " --output-precursors {precursors} --verbose"
    " --sigmas-multiplier {sigmas_multiplier} --transmission-geometry {transmission_geometry}"
    " --delta-left {delta_left} --delta-right {delta_right}"
    " && test -f {transpec}/schema.txt"
)
def transmit_precursors_into_fragment_space(
    tdf: BrukerD, clusters: PostprocessedPrecursorClusters,
    sigmas_multiplier: float, transmission_geometry: str, delta_left: float, delta_right: float,
):
    return TransmittedMs1Events[transpec], TransmittedPrecursorClusters[precursors]


@R.command("venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}")
def filter_first_precursors(precursors: TransmittedPrecursorClusters, filter: str):
    return FirstFilterPrecursors[filtered]


@R.command(
    ".venv/bin/necromerge2-stage-dir {staged}"
    " --link events.ms2={ms2}"
    " --link transmitted_precursors.mmappet={transprec}"
    " --link precursor_filter.mmappet={filter_mm}"
)
def stage_mkpmsms_input(ms2: Ms2Events, transprec: TransmittedMs1Events, filter_mm: FirstFilterPrecursors):
    return MkpmsmsInputStaged[staged]


@R.command(
    ".venv/bin/necromerge2-run-mkpmsms {staged} {pmsms} {precursors}"
    " --tofs-extraction-method {tofs_extraction_method}"
    " --tofs-extraction-params-json {tofs_extraction_params_json}"
    " --threads {threads}",
    threads=1,
)
def run_mkpmsms(staged: MkpmsmsInputStaged, tofs_extraction_method: str, tofs_extraction_params_json: str, threads: int):
    return Pmsms[pmsms], Ms2IndexedPrecursors[precursors]


@R.command("venvs/common/bin/filter_mmappet {precursors} {filtered} --verbose --filter {filter}")
def filter_pre_sage_precursors(precursors: Ms2IndexedPrecursors, filter: str):
    return PreSageFilteredPrecursors[filtered]


@R.command(
    "venvs/common/bin/build-precursor-grid-index {precursors} {tdf} {grid}"
    " --frame-mult {frame_mult} --scan-mult {scan_mult}"
    " --frame-inner-mult {frame_inner_mult} --scan-inner-mult {scan_inner_mult}"
    " --mz-inner-radius-da {mz_inner_radius_da} --top-k {top_k} --geometry {geometry}"
)
def build_precursor_grid_index(
    precursors: PreSageFilteredPrecursors, tdf: BrukerD,
    frame_mult: float, scan_mult: float, frame_inner_mult: float, scan_inner_mult: float,
    mz_inner_radius_da: float, top_k: int, geometry: str,
):
    return PrecursorGridIndex[grid]


@R.command(
    "git/ionmaidenmetal/build/precursor_neighbors_csr"
    " --boxes-input {grid_index}/boxes.mmappet --index-input {grid_index} --output {csr}"
    " $(venvs/common/bin/precursor-neighbors-params {tdf}"
    " --frame-mult {frame_mult} --scan-mult {scan_mult}"
    " --frame-inner-mult {frame_inner_mult} --scan-inner-mult {scan_inner_mult}"
    " --mz-inner-radius-da {mz_inner_radius_da} --top-k {top_k} --geometry {geometry})"
    " --n-threads {threads}",
    threads=1,
)
def compute_precursor_neighbors(
    grid_index: PrecursorGridIndex, tdf: BrukerD, threads: int,
    frame_mult: float, scan_mult: float, frame_inner_mult: float, scan_inner_mult: float,
    mz_inner_radius_da: float, top_k: int, geometry: str,
):
    return PrecursorNeighborsCsr[csr]


@R.command(
    "git/ionmaidenmetal/build/tof_filter --pmsms-path {pmsms} --neighbors-csr-path {neighbors_csr}"
    " --out-path {score} --n-threads {threads}",
    threads=1,
)
def tof_score_filter(pmsms: Pmsms, neighbors_csr: PrecursorNeighborsCsr, threads: int):
    return NeighborScore[score]   # no config arg -- confirmed vestigial in the Snakemake rule


@R.command(
    "venvs/common/bin/python -m timstofu.cli.score_based_pmsms_filter"
    " {pmsms} {precursors} {neighbor_score} {pmsms_out} {precursors_out}"
    " --threads {threads} --score-margin {score_margin}",
    threads=1,
)
def materialize_tof_filtered_pmsms(
    pmsms: Pmsms, precursors: PreSageFilteredPrecursors, neighbor_score: NeighborScore, threads: int,
    score_margin: float,
):
    return TofFilteredPmsms[pmsms_out], TofFilteredPrecursors[precursors_out]


@R.command(".venv/bin/necromerge2-stage-sage-input {pmsms} {tof2mz} {precursors} {staged}")
def make_tof_filtered_sage_input(pmsms: TofFilteredPmsms, tof2mz: Tof2Mz, precursors: TofFilteredPrecursors):
    return SageInputStaged[staged]


@R.command(
    ".venv/bin/necromerge2-run-sage {spectra} {fasta} {outdir} {run_info}"
    " --deisotope {deisotope} --chimera {chimera} --wide-window {wide_window} --predict-rt {predict_rt}"
    " --min-peaks {min_peaks} --max-peaks {max_peaks} --min-matched-peaks {min_matched_peaks}"
    " --max-fragment-charge {max_fragment_charge} --ignore-precursor-charge {ignore_precursor_charge}"
    " --parallel {parallel} --report-psms {report_psms}"
    " --isotope-errors-lo {isotope_errors_lo} --isotope-errors-hi {isotope_errors_hi}"
    " --database-bucket-size {database_bucket_size}"
    " --database-fragment-min-mz {database_fragment_min_mz} --database-fragment-max-mz {database_fragment_max_mz}"
    " --database-peptide-min-mass {database_peptide_min_mass} --database-peptide-max-mass {database_peptide_max_mass}"
    " --database-min-ion-index {database_min_ion_index} --database-max-variable-mods {database_max_variable_mods}"
    " --database-decoy-tag {database_decoy_tag} --database-generate-decoys {database_generate_decoys}"
    " --database-enzyme-missed-cleavages {database_enzyme_missed_cleavages}"
    " --database-enzyme-cleave-at {database_enzyme_cleave_at} --database-enzyme-restrict {database_enzyme_restrict}"
    " --database-enzyme-min-len {database_enzyme_min_len} --database-enzyme-max-len {database_enzyme_max_len}"
    " --database-enzyme-c-terminal {database_enzyme_c_terminal}"
    " --precursor-tol-ppm-lo {precursor_tol_ppm_lo} --precursor-tol-ppm-hi {precursor_tol_ppm_hi}"
    " --fragment-tol-ppm-lo {fragment_tol_ppm_lo} --fragment-tol-ppm-hi {fragment_tol_ppm_hi}"
    " --static-mods-json {static_mods_json} --variable-mods-json {variable_mods_json}"
)
def run_sage(
    spectra: SageInputStaged, fasta: Fasta,
    deisotope: bool, chimera: bool, wide_window: bool, predict_rt: bool,
    min_peaks: int, max_peaks: int, min_matched_peaks: int, max_fragment_charge: int,
    ignore_precursor_charge: bool, parallel: bool, report_psms: int,
    isotope_errors_lo: int, isotope_errors_hi: int,
    database_bucket_size: int, database_fragment_min_mz: float, database_fragment_max_mz: float,
    database_peptide_min_mass: float, database_peptide_max_mass: float, database_min_ion_index: int,
    database_max_variable_mods: int, database_decoy_tag: str, database_generate_decoys: bool,
    database_enzyme_missed_cleavages: int, database_enzyme_cleave_at: str, database_enzyme_restrict: str,
    database_enzyme_min_len: int, database_enzyme_max_len: int, database_enzyme_c_terminal: bool,
    precursor_tol_ppm_lo: float, precursor_tol_ppm_hi: float,
    fragment_tol_ppm_lo: float, fragment_tol_ppm_hi: float,
    static_mods_json: str, variable_mods_json: str,
):
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

    P.ms1_events = R.tdf2ms1(P.tdf)
    P.ms2_events, P.tof2mz = R.tdf2ms2(P.tdf)

    se = cfg["scale_estimation"]
    P.scale_estimates = R.estimate_precursor_scales(
        P.ms1_events, dataset=dataset,
        argmax_candidates_cnt=int(se["argmax_candidates_cnt"]),
        max_iter_smoothed_argmax_search=int(se["max_iter_smoothed_argmax_search"]),
        times_top_intensity_higher=float(se["times_top_intensity_higher"]),
        hill_descent_buffer=int(se["hill_descent_buffer"]),
        radii_tof=int(se["radii"]["tof"]), radii_urt=int(se["radii"]["urt"]), radii_scan=int(se["radii"]["scan"]),
        smoothing_radii_tof=int(se["smoothing_radii"]["tof"]), smoothing_radii_urt=int(se["smoothing_radii"]["urt"]),
        smoothing_radii_scan=int(se["smoothing_radii"]["scan"]),
        hq_pearson_tof=float(se["high_quality_min_pearson_corr"]["tof"]),
        hq_pearson_urt=float(se["high_quality_min_pearson_corr"]["urt"]),
        hq_pearson_scan=float(se["high_quality_min_pearson_corr"]["scan"]),
    )

    pcs = cfg["precursor_candidate_selection"]
    P.raw_precursor_clusters = R.select_precursor_candidates(
        P.ms1_events, P.scale_estimates, dataset=dataset,
        tof_radius=int(pcs["tof_radius"]),
        minimal_number_of_events_in_cluster=int(pcs["minimal_number_of_events_in_cluster"]),
        method=str(pcs["method"]),
        method_settings_json=_json_kwarg(pcs["method_settings"]),
    )

    pop = cfg["postprocessing_of_precursors"]
    dch = pop["decharging"]
    dchs = dch["method_settings"]
    P.postprocessed_precursor_clusters = R.postprocess_precursor_candidates(
        P.tdf, P.ms1_events, P.raw_precursor_clusters, P.scale_estimates, dataset=dataset,
        tof_radius=int(pop["tof_radius"]),
        precursor_connected_component_policy=str(pop["precursor_connected_component_policy"]),
        location_reestimator=str(pop["location_reestimator"]),
        estimate_frame_location_method=str(pop["estimate_frame_location_method"]),
        sigmas_multiplier=float(pop["sigmas_multiplier"]),
        overlap_geometry=str(pop["overlap_geometry"]),
        decharging_method=str(dch["method"]),
        decharging_precursor_no_charge_policy=str(dch["precursor_no_charge_policy"]),
        decharging_mass_diff=float(dchs["mass_diff"]),
        decharging_charges_to_look_for=",".join(str(c) for c in dchs["charges_to_look_for"]),
        decharging_max_jumps_right=int(dchs["max_jumps_right"]),
        decharging_min_ratio_of_previous_to_next_isotope_intensity=float(dchs["min_ratio_of_previous_to_next_isotope_intensity"]),
        decharging_force_same_tof_diffs=bool(dchs["force_same_tof_diffs"]),
        decharging_verbose=bool(dchs["verbose"]),
        decharging_sigmas_multiplier=float(dchs["sigmas_multiplier"]),
    )

    pt = cfg["precursor_transmission"]
    P.transmitted_ms1events, P.transmitted_precursor_clusters = R.transmit_precursors_into_fragment_space(
        P.tdf, P.postprocessed_precursor_clusters,
        sigmas_multiplier=float(pt["sigmas_multiplier"]), transmission_geometry=str(pt["transmission_geometry"]),
        delta_left=float(pt["delta_left"]), delta_right=float(pt["delta_right"]),
    )

    P.first_filter_precursors = R.filter_first_precursors(
        P.transmitted_precursor_clusters, filter=str(cfg["precursor_filters"]["mkpmsms"].get("filter", ""))
    )
    P.mkpmsms_input = R.stage_mkpmsms_input(P.ms2_events, P.transmitted_ms1events, P.first_filter_precursors)

    pmm = cfg["pseudomsms"]
    P.pmsms, P.ms2indexed_precursors = R.run_mkpmsms(
        P.mkpmsms_input, threads=threads,
        tofs_extraction_method=str(pmm["tofs_extraction_method"]),
        tofs_extraction_params_json=_json_kwarg(pmm["tofs_extraction_params"]),
    )

    P.pre_sage_filtered_precursors = R.filter_pre_sage_precursors(
        P.ms2indexed_precursors, filter=str(cfg["precursor_filters"]["pre_sage"].get("filter", ""))
    )

    pn = cfg["precursor_neighbors"]
    pn_kwargs = dict(
        frame_mult=float(pn["frame_mult"]), scan_mult=float(pn["scan_mult"]),
        frame_inner_mult=float(pn["frame_inner_mult"]), scan_inner_mult=float(pn["scan_inner_mult"]),
        mz_inner_radius_da=float(pn["mz_inner_radius_da"]), top_k=int(pn["top_k"]), geometry=str(pn["geometry"]),
    )
    P.precursor_grid_index = R.build_precursor_grid_index(P.pre_sage_filtered_precursors, P.tdf, **pn_kwargs)
    P.precursor_neighbors_csr = R.compute_precursor_neighbors(P.precursor_grid_index, P.tdf, threads=threads, **pn_kwargs)
    P.neighbor_score = R.tof_score_filter(P.pmsms, P.precursor_neighbors_csr, threads=threads)

    P.tof_filtered_pmsms, P.tof_filtered_precursors = R.materialize_tof_filtered_pmsms(
        P.pmsms, P.pre_sage_filtered_precursors, P.neighbor_score, threads=threads,
        score_margin=float(cfg["tof_score_filter"]["score_margin"]),
    )
    P.sage_input = R.make_tof_filtered_sage_input(P.tof_filtered_pmsms, P.tof2mz, P.tof_filtered_precursors)

    sg = cfg["sage"]
    db = sg["database"]
    enz = db["enzyme"]
    P.sage_raw_outdir, P.sage_run_info = R.run_sage(
        P.sage_input, P.fasta,
        deisotope=bool(sg["deisotope"]), chimera=bool(sg["chimera"]), wide_window=bool(sg["wide_window"]),
        predict_rt=bool(sg["predict_rt"]), min_peaks=int(sg["min_peaks"]), max_peaks=int(sg["max_peaks"]),
        min_matched_peaks=int(sg["min_matched_peaks"]), max_fragment_charge=int(sg["max_fragment_charge"]),
        ignore_precursor_charge=bool(sg["ignore_precursor_charge"]), parallel=bool(sg["parallel"]),
        report_psms=int(sg["report_psms"]),
        isotope_errors_lo=int(sg["isotope_errors"][0]), isotope_errors_hi=int(sg["isotope_errors"][1]),
        database_bucket_size=int(db["bucket_size"]),
        database_fragment_min_mz=float(db["fragment_min_mz"]), database_fragment_max_mz=float(db["fragment_max_mz"]),
        database_peptide_min_mass=float(db["peptide_min_mass"]), database_peptide_max_mass=float(db["peptide_max_mass"]),
        database_min_ion_index=int(db["min_ion_index"]), database_max_variable_mods=int(db["max_variable_mods"]),
        database_decoy_tag=str(db["decoy_tag"]), database_generate_decoys=bool(db["generate_decoys"]),
        database_enzyme_missed_cleavages=int(enz["missed_cleavages"]), database_enzyme_cleave_at=str(enz["cleave_at"]),
        database_enzyme_restrict=str(enz["restrict"]), database_enzyme_min_len=int(enz["min_len"]),
        database_enzyme_max_len=int(enz["max_len"]), database_enzyme_c_terminal=bool(enz["c_terminal"]),
        precursor_tol_ppm_lo=float(sg["precursor_tol"]["ppm"][0]), precursor_tol_ppm_hi=float(sg["precursor_tol"]["ppm"][1]),
        fragment_tol_ppm_lo=float(sg["fragment_tol"]["ppm"][0]), fragment_tol_ppm_hi=float(sg["fragment_tol"]["ppm"][1]),
        static_mods_json=_json_kwarg(db["static_mods"]), variable_mods_json=_json_kwarg(db["variable_mods"]),
    )
    P.results_json, P.sage_precursors, P.matched_fragments, P.sage_pin = R.extract_sage_results(P.sage_raw_outdir)
    P.sage_precursors_after_fdr = R.sage_filter(P.sage_precursors)
    P.sage_summary = R.sage_summarize(P.sage_precursors_after_fdr)
    return P
