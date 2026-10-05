# Recalibration: three independently selectable modes (B.6, 2026-08-20; gating simplified 2026-08-25)

`ionmaiden_pipeline`'s `if "sage" in cfg:` branch has three modes. Each of
mz/RT/IIM gates on whether its own config table exists — no separate
"is this feature on" flag, table presence *is* the flag (RT/IIM used to
share one `[recalibration.rt_iim]` umbrella table with a `dimensions`
sub-list; removed 2026-08-25, see below) — see
`plans/better_sage_filtering.md`'s B.6 for the full original design and
`plans/rt_iim_independent_dimensions.md` for the RT/IIM independence split:

1. **No recalibration** — no `[recalibration]` section at all. One SAGE
   pass, `run_sage` directly on `search_pmsms` (+ `tof2mz_table`)/`search_precursors`.
2. **mz recalibration alone** — `[recalibration]` present, neither
   `[recalibration.rt]` nor `[recalibration.iim]`. Existing two-pass mz
   correction (`recalibrate_pmsms_mz`/`recalibrate_precursors`/
   `update_sage_config`), unchanged since before B.6.
3. **mz + RT and/or IIM recalibration** — `[recalibration.rt]` and/or
   `[recalibration.iim]` also present (`"rt" in cfg.recalibration or "iim"
   in cfg.recalibration"`, not a `[recalibration.rt_iim]` gate anymore).
   Nests *inside* mode 2's branch (chains onto mode 2's already
   mz-corrected outputs, doesn't duplicate them): `predict_rt`/`predict_iim`
   (`git/featureprediction`, independently called per active dimension) run
   off the same `filtered_sage_results_tsv` anchors the mz fits use, then
   `correct_precursors_rt`/`correct_precursors_iim` chain onto
   `recalibrated_precursors`, then `update_sage_config_rt_iim` chains onto
   `update_sage_config`'s output to add `rt_tol_sec`+`rt_sigma_sec` and/or
   `mobility_tol`+`iim_sigma` (two `config_set` calls per active dimension,
   not one — see below) for whichever dimensions are active, then
   `run_sage` runs the final pass (see `run_sage_merge.md` — this used to
   be a distinct `run_sage_with_predicted` rule, since merged).

## Why table presence, not a separate `dimensions` list (2026-08-25)

`[recalibration.rt_iim].dimensions` used to be the sole way to select which
of RT/IIM were active, but by the time `tolerance_percentiles`/
`tolerance_method`/`min_charge`/`max_charge`/`server_url` had all moved into
their own `[recalibration.rt]`/`[recalibration.iim]` tables (the entries
just below), `dimensions` was carrying the exact same fact a second time —
and the name `[recalibration.rt_iim]` read as "both RT and IIM" while
actually meaning "the RT/IIM subsystem in general," which was a real
footgun reading a config cold. `dimensions` is now computed in the
workflow (`tuple(d for d in ("rt", "iim") if d in cfg.recalibration)`), not
read from config; `[recalibration.rt_iim]` no longer exists as a concept.
`NoPrediction`/the Python-callback-command pattern is unaffected by this —
that part remains necessary because necroflow's `@command` templates are
static and RT/IIM are independently optional (4 real combinations), even
though necroflow *does* support skipping a rule/node entirely via plain
`if/else` branching in the workflow (its `docs/rules.md`'s "Conditional
pipelines") — properly removing the sentinel would mean splitting
`update_sage_config_rt_iim` into up to 3 static-template rule variants
instead; not done, deliberately deferred as a larger, separate change.

`current_mz_pmsms`/`current_precursors` are three-way selections feeding
`convert_search_pmsms_to_mzml`/`convert_search_pmsms_to_mgf` — whichever
mode ran, exported MGF/mzML headers match exactly what SAGE1 actually
searched against. Before B.6, this variable didn't exist at all
(exports always used raw `search_precursors`, even in mode 2) — a
pre-existing gap this closes for mz too, not just RT/IIM. Named
`current_*`, not `final_*` (2026-09-07 rename): each is reassigned as
later recalibration/RT/IIM stages run, so it's only "final" once the
workflow function returns, not at its point of definition.

MGF export is optional: `convert_search_pmsms_to_mgf` is wired only when
`[mgf].config_path` is present. The path is a scalar rule input rather than a typed
artifact, so different paths produce different nodes; changing one config file in
place does not invalidate an existing node. `fragpipe_synthetic_pipeline` likewise
requires `[mgf].config_path` when `output_format = "mgf"`.

## `update_sage_config_rt_iim` writes `rt_sigma_sec`/`iim_sigma`, not just `rt_tol_sec`/`mobility_tol` (2026-08-26)

SAGE's `Input::build()` requires
`predicted_rt`+`rt_tol_sec`+`rt_sigma_sec` all-or-none (same for the IIM
trio, see the sage fork's `AI.md`) — `_update_sage_config_rt_iim_command`
previously only chained one `config_set` per active dimension (`rt_tol_sec`/
`mobility_tol`), so SAGE rejected the config the instant RT/IIM
recalibration mode ran end-to-end, on every job, forever (nothing before
this could have exercised mode 3 successfully). Now two `config_set` calls
per active dimension, both sourced from the same `rt_tolerance`/
`mobility_tolerance` artifact (`git/featureprediction`'s
`correct_precursors_rt`/`correct_precursors_iim` write `rt_sigma_sec`/
`iim_sigma` as sibling top-level keys in that same file — see that repo's
`AI.md`, which had the identical gap on its own side, fixed alongside this).
Found via a real end-to-end F9477 run while A/B-testing the sage fork's new
`combined_score` ranking, not by inspection — there was no test on either
side that would have caught this, since necroflow's `--requests` on smoke
configs used so far all resolved to nodes upstream of the final `run_sage`.

**Config derivation, not restatement** (see `git/featureprediction`'s
`AI.md` for the numbers): `[recalibration.iim]`'s `min_charge`/
`max_charge` default to `cfg.sage.get("precursor_charge", (2, 4))` (SAGE's
own compiled-in default) rather than crashing when `sage.precursor_charge`
is unset, which real job configs often do.

## `tolerance_percentiles`/`tolerance_method`: one table per dimension, not shared (2026-08-25)

`[recalibration.mz]`, `[recalibration.rt]`, and
`[recalibration.iim]` each carry their own `tolerance_percentiles`
(required explicitly, never inherited from a sibling dimension or from a
different table) and optional `tolerance_method` (`"theoretic"`, the
default — symmetric `median ± z*robust_sigma` — or `"empiric"`, plain
percentiles; see `git/featureprediction`'s `tolerance.select_tolerance`
and `git/searchops`'s `recalibration._select_tolerance`, separate
implementations of the same dispatch). Previously RT and IIM shared one
`[recalibration.rt_iim].tolerance_percentiles` value, and mz's own value
sat at `cfg.recalibration`'s root — a real F9477 finding
(2026-08-25: mass-error residuals are visibly right-skewed, RT residuals
much less so) showed treating all three the same was masking real
per-dimension differences. `[recalibration.rt_iim]` itself still exists as
the mode-3 trigger key (`if "rt_iim" in cfg.recalibration:`) and now only
carries `dimensions` (which of `rt`/`iim` are active) — `tolerance_lo`/
`tolerance_hi`/`min_charge`/`max_charge` are only read from
`[recalibration.rt]`/`[recalibration.iim]` when that dimension is actually
in `dimensions`, so a job enabling only one dimension doesn't need the
other table to exist at all. `predict_rt`/`predict_iim`/
`correct_precursors_rt`/`correct_precursors_iim` (the `@command` rules)
each gained a `tolerance_method: str` parameter, threaded to
`feature-prediction-*`'s `--tolerance-method` flag. mz needs no equivalent
Python-side change beyond the job-config nesting — `write_recalibration_config`
already serializes the whole `cfg.recalibration` dict verbatim as
`recalibration_config.toml`, so `searchops`'s own config-reading code is
what changed there, not this repo's.

**Note (2026-08-31): `tolerance_percentiles`'s `theoretic` method requires
`hi_pct` strictly between 50 and 100** — `[0, 100]` (a value several
committed job configs still use, sane only under the old `"empiric"`
default) now fails loudly in `searchops`/`git/featureprediction` rather
than silently producing an infinite tolerance window. See those repos'
`AI.md`s.

No back-compat shim for the old shared/root-level keys — existing job
configs were rewritten directly (5 of them set mz's `tolerance_percentiles`
using the current `fragment_model`/`precursor_model` schema:
`f9468_fragpipe.toml`, `f9477_gam_test.toml`, `b6699_gam_test.toml`,
`f9477_fragpipe.toml`, `short_test_recal.toml`). No committed job config
set `[recalibration.rt_iim]` at all as of this change — the RT/IIM mode-3
F9477 comparisons referenced elsewhere in this file were run via ad-hoc
CLI-level config overrides, not persisted job files, so there was nothing
to migrate for that key. Several other `jobs/*.toml` files
(`b6699_test_recal.toml`, `quick_test_10k_recal.toml`, `full_recal_grid.toml`,
`inspect_recalibrate.toml`, `quick_test_10k_recal_grid.toml`,
`f9468_test_recal.toml`, `short_test_recal_tol_grid.toml`,
`short_test_recal_shared_shape_grid.toml`, `short_test_recal_grid.toml`)
use an older, already-stale `recalibration` schema (flat `model = "..."`,
no `fragment_model`/`precursor_model` sub-tables) that predates the current
`searchops.recalibrate_pmsms_mz`/`recalibrate_precursors` split and would
already fail against current code independent of this change — left
untouched, not this change's problem to fix.

## `server_url` (2026-08-25)

`[recalibration.rt]`/`[recalibration.iim].server_url` — optional, either a
plain string or a TOML array (`server_url = ["ip0:8500",
"ip1:8500"]`, tried in that order, falling back on failure — see
`git/featureprediction`'s `AI.md` for the fallback design). `_server_url_arg`
(`pipelines.py`) normalizes either shape into the single comma-joined string
`predict_rt`/`predict_iim`'s `--server-url` flag expects (`feature-prediction-
generate-{rt,iim}` splits it back into a list), defaulting to
`_DEFAULT_KOINA_SERVER_URL` — a literal duplicate of `koina_client
.DEFAULT_SERVER_URL`'s value, since this repo shells out to a separate
`venvs/featureprediction` install and can't import that package directly
(same reasoning as `default_precursor_charge` duplicating SAGE's own
compiled-in `(2, 4)`, just above). Unlike `tolerance_percentiles`/
`tolerance_method`, `server_url` is *always* passed on the command line
(the `@command` template has no conditional-flag mechanism for static
templates), computed by Python first so the flag is never omitted, just
resolves to the same default `koina_client` itself would use when the job
config doesn't set it.

`ssl`/`timeout`/`retries`/`backoff_factor` (see `git/featureprediction`'s
`AI.md`) are **not** yet wired through `pipelines.py` — only `server_url`
was added at the pipeline-config level so far; those four are available at
the CLI/`feature_prediction.predict` config-dict level today, add them here
the same way if a job actually needs to set one.

## Mode-2/mode-3 branches flattened into one (2026-09-0x)

The `if "rt" in cfg.recalibration or "iim" in cfg.recalibration: <block>
else: <block>` split (what this file used to call "mode 2" vs "mode 3") no
longer exists as two separate code paths with two separate final `run_sage`
calls. RT and IIM are now independent sibling `if "rt" in dimensions:`/
`if "iim" in dimensions:` blocks at the same nesting level as the rest of
the `"recalibration" in cfg` branch, and there is exactly **one** final
`run_sage` call regardless of which of RT/IIM (zero, one, or both) are
active.

This was unlocked by `update_sage_config_rt_iim`'s own command builder
already supporting `dimensions=()` as a documented safe no-op (a plain `cp
{sage_config} {recalibrated_sage_config}`, per its own docstring, written
back when `dimensions` was assumed non-empty "in practice" but handled as a
safe fallback anyway) — so calling it unconditionally, even when neither
dimension is active, produces byte-identical config content to what
"mode 2" used to produce directly from `update_sage_config`'s output.
Verified on a real job: the realized `sage_config.json` for a mode-2 job
diffed byte-for-byte identical (modulo path strings) between the old
two-branch code and the new unified code, and real ion/PSM counts at 1%
FDR matched exactly (75,888 PSMs / 21,873 ions on the same F9477 job).

One real, one-time cost: every existing "mode 2" job's final `run_sage`
node gets a new hash the first time it runs after this change (its input,
`recalibrated_sage_config_rt_iim`, didn't exist as a concept for that mode
before), even though the config bytes are identical — a single cheap
`update_sage_config_rt_iim` node plus one `run_sage` rerun, not a
correctness concern.

The repeated `confident_psms → sage_pmsms_mapping → score_comparison`
trio (previously written out once per branch) is now one
`_finalize_confident_psms` closure inside `ionmaiden_pipeline`, called from
both the no-recalibration and recalibration code paths.

## Recalibration spectra from raw MS2 top cells: `[recalibration_top_cell]` (2026-09-29)

Optional table, gated by presence like `[fragment_intensity]`. When set, the
recalibration SAGE pass does not search the main pmsms. Each sampled precursor's
most probable (frame, scan) raw MS2 spectrum is searched in place:

- `write_tof2mz_table(tdf)` -> `Tof2MzTable`: the dense tof -> m/z table
  (float32 column `mz`; `timstofu.cli.write_tof2mz_table`). Since 2026-10 every
  `run_sage` call takes it (next section).
- `top_cell_precursors(recalibration_precursors, transmitted_ms1events, ms2_events)`
  -> `TopCellPrecursors`: the same sample, with `fragment_spectrum_start` /
  `fragment_event_cnt` pointing at the first footprint row's (frame, scan) slice of
  `events.ms2/data.mmappet` (`timstofu.cli.top_cell_precursors`; footprint rows are
  sorted by probability, most probable first).
- `run_sage(ms2_events, top_cell_precursors, ..., tof2mz=tof2mz_table)`: `run_sage`
  takes either a pmsms or the raw `Ms2Events` store (told apart by the input's
  `events.ms2` filename); for the latter it runs Sage with
  `--pmsms <events.ms2>/data.mmappet --tof2mz <table>`.

The fits (`recalibrate_pmsms_mz`, `recalibrate_precursors`, RT) are unchanged and
still apply to the main pmsms and `search_precursors`: they read only the
recalibration search's PSMs and matched fragments. Without the table the recalibration
`run_sage` gets the same nodes as before.

These are exactly the spectra mkpmsms' `top_probable_frame_scan` copies into a pmsms.
The mkpmsms route (`[recalibration_pseudomsms]`, commit 81ec5cb, replaced by this)
built them with a second mkpmsms run, cut-and-index and a materialized m/z column;
SAGE on the raw top cells reproduced that search exactly on F9477 (all 47,917 PSMs and
235,629 matched fragments identical), and the pipeline's two tables equal the ones
used there. Comparison with `f9477_best`:
`git/pipeline_analysis/docs/ai/recalibration_paths.md`.

## Fragment m/z correction applied inside SAGE (2026-10)

Plans: necromerge2 `plans/fragment_mz_correction_in_sage.md`, simplified by
`plans/fragment_mz_corrected_tof2mz_table.md`. Every `run_sage` reads fragment m/z
as `tof2mz[tof]` from the pmsms' `tof` column; no search reads a materialized `mz`
column (Sage's `mz`-column input was removed).

- `recalibrate_pmsms_mz(filtered results, filtered matched fragments, tof2mz_table,
  search_precursors)` fits `bias + f_mz(mz) + f_rt(rt)` and writes no pmsms.
  Outputs: `recalibrated_tof2mz_table` (`RecalibratedTof2MzTable`, float64,
  `table[t] / (1 + f_mz(table[t])·1e-6)`), `fragment_shifted_precursors`
  (`search_precursors` plus `fragment_shift_ppm = bias + f_rt(raw rt)`), tolerance
  and plot.
- `recalibrate_precursors` reads `fragment_shifted_precursors`, so the column rides
  through every later precursor rewrite (m/z, RT, IIM; each rewrites the whole
  table) into the final search's `current_precursors`.
- Final `run_sage(search_pmsms, current_precursors, ...,
  tof2mz=recalibrated_tof2mz_table)`: Sage reads each fragment m/z as
  `table[tof] / (1 + fragment_shift_ppm·1e-6)`, applying the shift whenever the
  precursors table has the column (git/sage `docs/ai/pmsms_input.md`). The
  recalibration search uses the raw `tof2mz_table` and precursors without the column.
- Exports and `sage_map_to_pmsms` read m/z the same way, from `search_pmsms` with
  `--tof2mz` (since 2026-10-05, `plans/exports_from_tof2mz_table.md`):
  `convert_search_pmsms_to_mzml`/`convert_search_pmsms_to_mgf` get
  `current_tof2mz_table` + `current_precursors` (raw table and `search_precursors` in
  mode 1, recalibrated table and the shifted precursors after recalibration), and
  `sage_map_to_pmsms` gets the table + `fragment_shifted_precursors`. Nothing in the
  pipeline materializes `mz`; the `materialize_pmsms_mz`/
  `materialize_recalibrated_pmsms_mz` rules and `MzPmsms`/`RecalibratedPmsms` are
  gone. timstofu's `materialize_pmsms_mz` CLI remains as a standalone tool.
  Verified on F9477 (2026-10-05) against materialize-then-write: MGF byte-identical,
  `--indexed` mzML byte-identical, pipeline (directory-mode) mzML identical per
  precursor across all 923,634 spectra (its write order is run-dependent either way),
  `sage_map_to_pmsms` tables identical. Times: mzML 47.4 s vs 8.8 s + 47.8 s, MGF
  88.5 s vs 8.8 s + 84.9 s, mapping 8.2 s vs 8.8 s + 8.4 s, with no 6.5 GB `mz` column.
- `convert_search_pmsms_to_mgf` runs `chattr +m {workdir}` first (a no-op off btrfs):
  compressing ~24 GiB of text with zstd during writeback throttled the writer. F9477
  MGF step 90.9 s -> 74.6 s, byte-identical output; the file then takes its full size
  on disk.

The table is float64 because a float32 one would round every m/z twice: measured on
F9477, that moves 25% of peaks by one float32 ulp against the single rounding, while
float64 moves essentially none.

Verified on `jobs/f9477_best.toml` (2026-10-05), paired with the `.mzcalib` route's
fresh run on the same first pass:

- Identical fit (`fragment_shift_ppm` bit-identical). 1,434,812 of 1,625,763,451
  peaks (0.088%) move by one float32 ulp, almost all from evaluating `f_mz` exactly
  instead of off the old 2,000-point grid (interpolation error up to 6e-4 ppm);
  dividing in two steps instead of adding ppm accounts for 40,117.
- Final SAGE: 2,038 of 2,245,145 matched m/z differ (one ulp, five by more), 6/7
  PSM rows differ; SAGE-level peptides/ions at 1% `peptide_q` 21,123 / 24,653 on
  both routes (3 peptides, 4 ions swapped). mokapot 27,478 -> 28,061 peptides,
  33,597 -> 34,290 ions: its input sensitivity amplifying those changes, not an
  effect of the route.
- Export: all 2,119,448 charge-1 matched peaks found bit-exact in
  `materialize_recalibrated_pmsms_mz`'s output, which now takes 8.4 s (one pass,
  no `.d`) instead of ~22 s.
- Steps: `recalibrate_pmsms_mz` 12.5 s (15.2 s before), final SAGE 79.5 s (79.7 s).

Why the shift is computed here and not by Sage from RT: in mode 3,
`correct_precursors_rt` overwrites `rt` with corrected RT before the final search,
while `f_rt` is fitted on raw RT. Evaluating it once, from raw RT, before any
rewrite avoids that. The column is also where a different per-precursor fragment
model would go without touching Sage (e.g. a 2-D `f(rt, 1/K0)`, or a free
per-precursor intercept).

First version (SAGE evaluated an `mz`-only `.mzcalib` per peak), verified on
`jobs/f9477_best.toml` (2026-10-02) against the materialized route's run of the same
job file from 2026-10-01:

- Final SAGE output: the same 358,374 PSMs (scannr, peptide, rank, charge) and the
  same 2,245,176 matched fragments; 15 matched m/z differ by exactly one float32 ulp,
  which moved `fragment_ppm`/`matched_intensity_pct`/`poisson`/`scored_candidates` on
  35 PSMs and the LDA score by ~1e-7. SAGE-level peptides/ions at 1% `peptide_q`:
  identical sets (21,124 / 24,656).
- Cause of those ulps, measured over all 1,625,763,451 peaks: the fit is identical,
  but `f_rt` used to be read off a 2,000-point linear grid at float32 RT and is now
  the P-spline evaluated exactly at float64 RT (difference ≤7.6e-6 ppm, median
  6.6e-7). That flips the float32 rounding of 20,720 peaks (1.3e-5) by one ulp.
  The changed summation order (`f_mz + (bias + f_rt)`) alone flips none.
- The old route's `.mzcalib` stored `bias = 0` while applying the real bias
  (−3.03 ppm on this run) in memory; the bias now lives in `fragment_shift_ppm`.
- mokapot: 27,679 -> 27,746 peptides, 33,674 -> 33,904 ions (+0.7%). Not a real
  effect: mokapot (seed 1) reproduces its output byte for byte on an identical pin,
  but applying only the 35 rows' tiny feature changes to the old pin already moves
  it by +237 peptides / +205 ions, and the PSM row order (SAGE sorts by its LDA
  score, which moved ~1e-7) shifts it further.
- Export bit-identity: all 2,119,447 charge-1 matched peaks' m/z (after SAGE's own
  `(mz - PROTON) + PROTON` float32 round trip) occur exactly in the same precursor's
  slice of `materialize_recalibrated_pmsms_mz`'s output; against the uncorrected
  `search_mz_pmsms` only 0.06% do.
- Time and disk: `materialize_pmsms_mz` (20.4 s, a 6.5 GB `mz` column) no longer
  runs, and `recalibrate_pmsms_mz` no longer writes the second 6.5 GB column; it
  still takes ~13-20 s, dominated by its diagnostic plot (~5 s), the fit (~4 s),
  imports and the `tof` range scan (~2.4 s). Final SAGE file IO is unchanged
  (15.8 -> 15.6 s) with the lookup and correction done per peak. The recalibration
  SAGE's file IO rose 6.4 -> 13.3 s in this incremental run because the `tof`
  column was cold in the page cache (the old route read an `mz` column
  `materialize_pmsms_mz` had just written); rerun warm it is 5.6 s.
