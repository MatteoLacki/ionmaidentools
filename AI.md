# AI.md — ionmaidentools

## What this repo is

`src/ionmaidentools/pipelines.py` is the necroflow pipeline factory for the ionmaiden
DIA timsTOF pipeline — the single module every job TOML's `.pipeline` key points at
(e.g. `git/ionmaidentools/src/ionmaidentools/pipelines.py:ionmaiden_pipeline`). It has
no other source files: one large module, `NodeType` subclasses + `@command`-decorated
rule functions + two top-level pipeline factory functions.

Package metadata (`pyproject.toml`): depends on `tomlkit` + `dictodot`, no CLI entry
points registered — this package is imported by necroflow's job runner (`./nf`, i.e.
`.venv/bin/necroflow`), not invoked directly.

## Necroflow rule conventions (as used here)

- A rule is a plain function decorated `@command("<shell template with {placeholders}>")`.
  Parameters annotated with a `NodeType` subclass are DAG-edge inputs (parent nodes);
  plain-typed parameters (`str`, `int`, `float`) are scalar config/CLI-flag inputs, not
  DAG edges.
- `output(SomeNodeType)` inside the function body declares a typed output node; its
  placeholder is substituted into the shell template. The function returns the output
  node(s).
- Pipeline factories (`ionmaiden_pipeline(P, config)`, `fragpipe_synthetic_pipeline(P, config)`)
  wire rules together: `P.some_label = some_rule(P, parent_node, ..., scalar_kwarg=...)`.
  The first positional arg is always the `Pipeline` object `P`.
- Both factories serialize their complete Necroflow-supplied config through
  `write_pipeline_config` and publish it as requestable `P.pipeline_config`
  (`pipeline_config.toml`). Necroflow strips metadata keys such as `.pipeline` and
  `.requests`, resolves `.extends`, and expands grids before calling the factory, so
  the file records the effective factory input.
- Every rule shells out to a real installed/patched CLI directly — no bespoke Python
  wrapper CLIs. External tools live at fixed, pre-installed paths (`software/sage/devel_fixed`,
  `software/fragpipe/fragpipe-24.0`), same convention on both sides.
- See root `CLAUDE.md`'s "Precursor table format" — all precursor tables in this pipeline
  are `.mmappet` directories (`MmappetDataset` subclasses here), never intermediate parquet.

## Two pipeline factories

- `ionmaiden_pipeline`: the real pipeline — Bruker `.d` → MS1 peak picking → quadrupole
  transmission → pseudo-MS/MS → SAGE search (+ optional recalibration pass, optional
  MGF export when the job config has an `[mgf]` section, optional FragPipe side-by-side
  comparison when it has a `[fragpipe]` section). `[mgf].config_path` is passed as a
  scalar input to `convert_search_pmsms_to_mgf`; its value therefore selects and
  fingerprints the MGF flavour. Jobs without that key do not add `search_mgf` to the
  pipeline.
- `fragpipe_synthetic_pipeline`: FragPipe smoke test on Koina-simulated peptides from a
  FASTA (no Bruker `.d` input, no Sage) — reuses `ionmaiden_pipeline`'s FragPipe rules
  verbatim (`source_fragpipe_workflow`/`generate_fragpipe_decoy_fasta`/
  `patch_fragpipe_workflow`/`write_fragpipe_manifest`/`run_fragpipe`/`extract_fragpipe_log`/
  `summarize_fragpipe`).

## C++ MS1 extraction

`ionmaiden_pipeline` uses `git/ionmaidenmetal/build/tdf2ms ms1` for `ms1_events`,
passing the Necroflow-allocated thread count and `--overwrite`. The rule validates both
mmappet split-index datasets and the event-dataset schema before completion. The in-memory
slice contract is unchanged from `d2ms1`; both indices now use the common mmappet-array
format, so all peak-picking and precursor consumers retain the same `(tof,urt,scan)` slices.

The converter decompresses each MS1 frame once, stable-count-sorts it to a `(tof,scan)` run, then
uses event-balanced TOF shards to merge runs into disjoint sequential mmap regions while building
the split index. Exact comparisons with `d2ms1` passed on F9477, F9468, and B6699; observed wall
speedups were 3.78x, 2.76x, and 3.43x. Python `d2ms1` remains installed as a reference/debug tool,
not the pipeline default.

## Compact TOF/scan/frame MS2 extraction benchmark target

`ionmaiden_pipeline` exposes `ms2_tsf_events` as an independently requestable target.
Its `tdf2ms2_tsf` rule runs `git/ionmaidenmetal/build/tdf2ms ms2-tsf` with the
Necroflow-allocated thread count and produces `events_ms2_tsf.mmappet`. The command
validates the payload and flat TOF/scan row-start mmappet schemas plus `stats.json`
before Necroflow marks the node complete.

This target deliberately has no downstream consumer yet. Existing pseudo-MS/MS stages
continue to use `ms2_events` and the established frame/scan layout; request
`ms2_tsf_events` explicitly to benchmark or inspect the compact `(tof,scan,frame)`
layout without changing production results. Necroflow records wall time and output size
automatically in the node `.rip/run.toml`.

## Query-oriented TOF/frame/scan MS2 extraction benchmark target

`ionmaiden_pipeline` also exposes `ms2_tfs_events` independently. Its
`tdf2ms2_tfs` rule runs `git/ionmaidenmetal/build/tdf2ms ms2-tfs` with the
Necroflow thread allocation and publishes `events_ms2_tfs.mmappet`. Completion
checks cover intensity data, packed-scan schema/shape, TOF index schema/shape,
frame index schema, and `stats.json`.

This node has no downstream consumer and does not replace `ms2_events` or
`ms2_tsf_events`. Request it to benchmark frame-selective box queries against the
`(tof,frame,scan)` layout. Scans use four little-endian uint10 values per five bytes;
TOF and frame headers locate the narrow ranges without expanding scan IDs. F9477
completed in 5.07 s and 1.554 GB, with all 160,689,740 events matching TSF exactly.

## Fragment-intensity prediction cache: this repo's first `mutable=True` rule (2026-08-31)

`predict_fragment_intensity` wraps `git/featureprediction`'s
`feature-prediction-generate-fragments` (see that repo's `AI.md` for the
Prosit/`compact_trt` predictor itself and the cache's own internal shape).
Independently requestable (same pattern as `ms2_tsf_events`/`ms2_tfs_events`
above) off `dumped_peptides` — never runs just because a job dumps peptides,
since a full human-proteome fill is a real, expensive, network-bound
operation (~43 minutes against the live Koina server measured on F9477).

**Its output (`FragmentIntensityCache`) is declared `mutable=True`** —
necroflow's own mechanism for "persistent single-output state whose external
byte changes should not invalidate consumers" (`git/necroflow`'s
`docs/rules.md`, "Mutable Rules"). This is the first use of that feature in
this repo. It exists because the underlying cache
(`mmappeteer.PredictionCache`) is deliberately append-only and growing
across runs — `feature_prediction.fragment_intensity.predict_fragment_intensity`
does `cache.lookup()` before ever calling Koina, so a rerun against the same
node only fills in genuinely missing keys. Without `mutable=True`, that
growth between runs would be indistinguishable from "this output changed,
re-run everything downstream" to a normal (immutable, content-addressed)
necroflow node — exactly the semantics this feature exists to opt out of.
Real effect verified: the cache was first populated (a real full-F9477
fill, ~2.8GB) by a standalone script, then moved by hand into this node's
exact hash directory. The first real `necroflow run` invocation of
`predict_fragment_intensity` against that pre-populated cache (necroflow's
own summary: `1 completed, 4 skipped (up-to-date)`) still had to run
(this was its first-ever necroflow-tracked invocation, so it can't be
skipped), but completed in 131.9s — dominated entirely by the cache's own
bulk SQLite lookup over 16.87M keys, confirming the underlying Python
function's lookup-before-append logic made zero Koina calls against
already-cached data.

No downstream consumer yet (SAGE doesn't read this cache) — that's the next
piece of `plans/sage_features_on_while_searching.md`, not done here. See
`plans/fragment_intensity_cache.md` for the cache's own design (predictor
choice, slot-numbering scheme, dtype choices).

### `export_fragment_intensity_for_sage`: ordinary rule consuming a mutable parent (2026-08-31)

Wraps `git/featureprediction`'s `feature-prediction-export-fragments-for-sage`
(DuckDB-based; see that repo's AI.md) — scopes the shared, ever-growing
`FragmentIntensityCache` down to one job's `dumped_peptides` x charge range,
as a single parquet file (`FragmentIntensityForSage`). An **ordinary**
(non-mutable) rule, even though one of its inputs
(`predict_fragment_intensity`'s output) is `mutable=True` — per
`docs/rules.md`'s "Mutable Rules", "if a mutable call executes during the
current run, every consumer replays", so necroflow itself guarantees this
always sees whatever cache state `predict_fragment_intensity` left behind in
the same invocation, with no extra bookkeeping needed here. Same
independently-requestable reasoning as its parent — no downstream SAGE
consumer yet, never runs just because `dumped_peptides` exists.

Verified end-to-end against the real full-F9477 setup (2026-08-31,
`jobs/f9477_gam_test.toml` with `.requests = ["fragment_intensity_for_sage"]`):
dry-run correctly showed only this one node would run (5 already up-to-date,
including the mutable cache itself — confirming the mutable-consumer
dependency doesn't force a spurious recompute of its parent). Real run:
116.6s, exit 0, 1.9GiB output parquet, `2,028,474` of `18,902,646`
(sequence, charge) pairs omitted (too-long peptides / not yet filled) —
identical warning counts to `git/featureprediction`'s own standalone-script
verification of the same export logic. Re-verified after that repo's
`ORDER BY ce.start` fix (see its AI.md) via `--invalidate
fragment_intensity_for_sage`: 62.4s, same row count, same warning counts,
output byte-for-byte re-checked against `PredictionCache.lookup()`.
**Re-verified again after the pointer-export redesign** (exports `sequence,
charge, start, end` instead of copying the sparse payload — see that
repo's AI.md): 18.0s, 479MiB (down from 1.9GiB), same row/warning counts,
output re-checked against `PredictionCache.lookup()` + direct
`arrays.mmappet` reads at the resolved ranges.

## Recalibration: three independently selectable modes (B.6, 2026-08-20; gating simplified 2026-08-25)

`ionmaiden_pipeline`'s `if "sage" in cfg:` branch has three modes. Each of
mz/RT/IIM gates on whether its own config table exists — no separate
"is this feature on" flag, table presence *is* the flag (RT/IIM used to
share one `[recalibration.rt_iim]` umbrella table with a `dimensions`
sub-list; removed 2026-08-25, see below) — see
`plans/better_sage_filtering.md`'s B.6 for the full original design and
`plans/rt_iim_independent_dimensions.md` for the RT/IIM independence split:

1. **No recalibration** — no `[recalibration]` section at all. One SAGE
   pass, `run_sage` directly on `search_mz_pmsms`/`search_precursors`.
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
   `run_sage_with_predicted` (a
   distinct rule from plain `run_sage`, adding `--predicted-rt`/
   `--predicted-iim` — necroflow `@command` templates are static, can't
   conditionally add a flag, so this uses a Python command callback
   instead, see `_run_sage_with_predicted_command`).

**Why table presence, not a separate `dimensions` list (2026-08-25):**
`[recalibration.rt_iim].dimensions` used to be the sole way to select which
of RT/IIM were active, but by the time `tolerance_percentiles`/
`tolerance_method`/`min_charge`/`max_charge`/`server_url` had all moved into
their own `[recalibration.rt]`/`[recalibration.iim]` tables (the entries
just below), `dimensions` was carrying the exact same fact a second time —
and the name `[recalibration.rt_iim]` read as "both RT and IIM" while
actually meaning "the RT/IIM subsystem in general," which was a real
footgun reading a config cold. `dimensions` is now computed in the pipeline
factory (`tuple(d for d in ("rt", "iim") if d in cfg.recalibration)`), not
read from config; `[recalibration.rt_iim]` no longer exists as a concept.
`NoPrediction`/the `run_sage_with_predicted`/`update_sage_config_rt_iim`
Python-callback-command pattern is unaffected by this — that part remains
necessary because necroflow's `@command` templates are static and RT/IIM
are independently optional (4 real combinations), even though necroflow
*does* support skipping a rule/node entirely via plain `if/else` branching
in the factory (its `docs/rules.md`'s "Conditional pipelines") — properly
removing the sentinel would mean splitting `run_sage_with_predicted`/
`update_sage_config_rt_iim` into up to 3 static-template rule variants
instead; not done, deliberately deferred as a larger, separate change.

`final_mz_pmsms`/`final_precursors` are three-way selections feeding
`convert_search_pmsms_to_mzml`/`convert_search_pmsms_to_mgf` — whichever
mode ran, exported MGF/mzML headers match exactly what SAGE1 actually
searched against. Before B.6, `final_precursors` didn't exist at all
(exports always used raw `search_precursors`, even in mode 2) — a
pre-existing gap this closes for mz too, not just RT/IIM.

MGF export is optional: `convert_search_pmsms_to_mgf` is wired only when
`[mgf].config_path` is present. The path is a scalar rule input rather than a typed
artifact, so different paths produce different nodes; changing one config file in
place does not invalidate an existing node. `fragpipe_synthetic_pipeline` likewise
requires `[mgf].config_path` when `output_format = "mgf"`.

**`update_sage_config_rt_iim` writes `rt_sigma_sec`/`iim_sigma`, not just
`rt_tol_sec`/`mobility_tol` (2026-08-26).** SAGE's `Input::build()` requires
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
configs used so far all resolved to nodes upstream of `run_sage_with_predicted`.

**Config derivation, not restatement** (see `git/featureprediction`'s
`AI.md` for the numbers): `[recalibration.iim]`'s `min_charge`/
`max_charge` default to `cfg.sage.get("precursor_charge", (2, 4))` (SAGE's
own compiled-in default) rather than crashing when `sage.precursor_charge`
is unset, which real job configs often do.

**`tolerance_percentiles`/`tolerance_method`: one table per dimension, not
shared (2026-08-25)** — `[recalibration.mz]`, `[recalibration.rt]`, and
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

**`[recalibration.rt]`/`[recalibration.iim].server_url` (2026-08-25)** —
optional, either a plain string or a TOML array (`server_url = ["ip0:8500",
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
templates — see `_run_sage_with_predicted_command`'s docstring on that
limitation), computed by Python first so the flag is never omitted, just
resolves to the same default `koina_client` itself would use when the job
config doesn't set it.

`ssl`/`timeout`/`retries`/`backoff_factor` (see `git/featureprediction`'s
`AI.md`) are **not** yet wired through `pipelines.py` — only `server_url`
was added at the pipeline-config level so far; those four are available at
the CLI/`feature_prediction.predict` config-dict level today, add them here
the same way if a job actually needs to set one.

## `run_sage`/`run_sage_with_predicted` merged into one rule (2026-08-26)

`run_sage_with_predicted` (a separate rule, `--predicted-rt`/`--predicted-iim`
always real DAG edges resolving to the `NoPrediction` sentinel when a
dimension was inactive, flag inclusion driven by a scalar `dimensions`
tuple) no longer exists — `run_sage` itself now takes
`predicted_rt: PredictedRt | None = None`/`predicted_iim: PredictedIim | None
= None` directly, using necroflow's "mixed Node/value inputs" support
(`docs/rules.md`) added after the original RT/IIM independent-dimensions
split. One rule now covers pass-1/mode-1/mode-2 (neither prediction,
`predicted_rt`/`predicted_iim` simply omitted) and mode 3's final pass
(either/both, passed as real Nodes) — the pipeline factory's mode-3 `else`
branches assign `P.predicted_rt = None`/`P.predicted_iim = None` directly
(a plain Python attribute, not a DAG edge — confirmed safe by reading
`Pipeline.__setattr__` itself: non-`Node` values always fall through to a
plain instance attribute, no validation, no DAG registration) instead of
`write_no_prediction_marker(...)`.

`_run_sage_command` (a Python command callback, not a static template,
since necroflow's `@command` string templates can't conditionally include
a flag) decides which `--predicted-*` flags to add purely from
`args.inputs.predicted_rt/predicted_iim is not None` — `args.inputs.*`
preserves a plain `None` verbatim for an unresolved mixed-input default
(only managed Nodes resolve to a real `Path`), so this needs neither the
old `dimensions` scalar nor the `NoPrediction` sentinel just to keep the
parameter a real DAG edge. `NoPrediction` itself is not removed — it's
still required for `update_sage_config_rt_iim`'s `rt_tolerance`/
`mobility_tolerance` inputs, which stayed plain required `NodeType`s (out
of scope for this merge; properly avoiding the sentinel there would mean
splitting that rule into up to 3 static-template variants instead of one
Python-callback command).

Verified against real F9477 data end-to-end (`ab_evict`-style RT-recalibration
job, both `run_sage` invocations it produces): the filtered/mode-2 pass's
realized Sage config has `"predicted_rt": null, "predicted_iim": null,
"rt_tol": null` (flag genuinely absent, 26,673 target PSMs at 1% FDR); the
final mode-3 pass's has `"predicted_rt": ".../predicted_rt.parquet"` with a
real path, Sage's own log confirming `"loaded 6300882 predicted sequence ->
rt entries"`, and populated `rt_sigma`/`rt_tol` (90,125 target PSMs at 1%
FDR). Both branches of the merged rule's Python callback exercised for
real, not just dry-run/import-checked.

## Conventions to follow when editing this file

- Design rationale, migration history, and "why this and not that" reasoning belongs
  here (or in this repo's own future notes), **not** in `pipelines.py` docstrings —
  keep code docstrings about what the code does now, not why it changed.
- FragPipe's `.workflow` file is a hand-maintained flat `key=value` properties file
  (`configs/search/fragpipe/workflows/*.workflow`), symlinked in as-is by
  `source_fragpipe_workflow`. `patch_fragpipe_workflow` derives a copy with
  `database.db-path=` rewritten (via `sed`) to point at `generate_fragpipe_decoy_fasta`'s
  output, so `run_fragpipe` always searches against a decoy-augmented FASTA — this is
  the one exception to "never patched", added because FragPipe/Philosopher has no
  generate-decoys-at-search-time option (unlike Sage) and there's no CLI override flag
  for it.

## Fragment-intensity cache wired into `run_sage`, `[fragment_intensity]` gate (2026-08-31)

`export_fragment_intensity_for_sage`/`predict_fragment_intensity` (see their own
docstrings above) had no downstream consumer until now — the Rust-side reader
(`software/sage/devel_fixed`'s `--predicted-fragment-intensity-index`/
`-cache`, MSBooster-parity `ms2_*` features) existed but nothing in this
pipeline ever passed those flags, so every job's `ms2_*` PIN/TSV columns were
real column shapes computing on empty input — all `0.0`, verified directly on
real F9477 output before this change.

`run_sage` gained two more optional mixed Node/`None` inputs,
`predicted_fragment_intensity_index: FragmentIntensityForSage | None`/
`predicted_fragment_intensity_cache: FragmentIntensityCache | None`, both-or-
neither (mirrors Rust's own `Input::build()` validation). `_run_sage_command`
passes `--predicted-fragment-intensity-cache
{cache_node}/arrays.mmappet` — the Rust reader only ever wants the shared
cache directory's `arrays.mmappet` subdirectory, never `index.sqlite3`/
`write.lock`.

**Gated by a new `[fragment_intensity]` presence-only config table**, computed
once (`_final_pass_fragment_intensity_index`/`_cache`, `None`/`None` when the
table is absent) and threaded into all three *final*-pass `run_sage` call
sites (mode-1 no-recalibration, mode-2 mz-only, mode-3 rt/iim) — never the
calibration-anchor pass. Without this gate, threading the cache in
unconditionally would make every sage job transitively depend on
`predict_fragment_intensity`, breaking that rule's own documented "never runs
just because `dumped_peptides` exists" invariant. Same table-presence-as-flag
convention as `[recalibration.rt]`/`[recalibration.iim]`.

**Positional, not keyword, gotcha**: these two params are Node-shaped (mixed
Node/`None`) inputs, so necroflow's `Rule._validate_input_presence` requires
them passed *positionally* like `predicted_rt`/`predicted_iim` — passing them
as `keyword=` args (as a first attempt at the mode-1/mode-2 call sites did)
raises `TypeError: run_sage: unexpected inputs: [...]` at pipeline-factory
time, since necroflow classifies each rule input as strictly one or the other
class based on its type annotation (`_node_input_contract` returning
non-`None`), never both. When `predicted_rt`/`predicted_iim` also need to be
`None` at a call site that isn't mode-3 (to reach the trailing positional
fragment-intensity args), pass literal `None, None` positionally — passing an
explicit `None` for a mixed Node/`None` **positional** input is fine and
already exercised elsewhere (mode-3's single-dimension-active branches); only
plain scalar (`_kw_inputs`) config values can't safely be `None` (see below).

## `[mokapot].plugin`: config-driven model choice, not hardcoded (2026-08-31)

The plain-SAGE `mokapot(...)` call site's `plugin` kwarg used to be a hardcoded
`"xgboost"` literal. Now `cfg.mokapot.get("plugin", "") if "mokapot" in cfg
else ""` — a job can select `[mokapot] plugin = "xgboost"` or omit the table
entirely for mokapot's own default (linear SVM) model. The sagepy_rescore
branch's own separate `mokapot(...)` call is untouched, still always
`plugin="xgboost"` (unrelated call, no config knob added there).

**`""`, never `None`, for "no plugin" — necroflow can't serialize a bare
`None` scalar (`_kw_inputs`) value.** `plugin`/`rt_source`/`iim_source` are
plain `str | None`-typed scalar rule params, not Node-shaped, so they're
recorded verbatim in `dependencies.toml` (tomlkit-serialized) for provenance
whenever a call site actually passes the keyword (even a value that happens
to equal the Python-level default) — a real, explicit `None` there raises
`tomlkit.items._ConvertError: Invalid type <class 'NoneType'>` and crashes
the whole run *after* the underlying work (e.g. a real mokapot fit) already
completed. Found via a real 2x2x3 ablation grid job (the first time this
call site's `plugin` was ever computed as `None` instead of always
`"xgboost"`). `""` is falsy in `_mokapot_command`'s own `if
args.config.plugin` check, so the emitted command is identical (`--plugin`
flag omitted) — only the provenance-recording path differs.

## Leakage-safe mokapot PIN filtering (2026-08-31)

See `plans/mokapot_leakage_safe_pin.md` for the full design. Summary: SAGE's
raw `results.sage.pin` contains `posterior_error` (SAGE's own in-run,
label-conditioned LDA score) — feeding it to mokapot as a feature is real
leakage. `scripts/mokapot_pin_adapter.py --mode sage` projects the PIN down to
a fixed, hardcoded-safe column registry via DuckDB before mokapot ever sees
it; `--mode passthrough` (unchanged original behavior, just drops `FileName`)
still serves the sagepy_rescore branch's already-filtered PIN. `mokapot()`
gained `rt_source`/`iim_source: str | None = None` (kw_inputs, selects which
external-prediction column pair — if any — the safe-PIN filter includes;
`"external"` exactly when that dimension's real predicted-property Node
exists for this job, `"none"` otherwise) threaded through
`_mokapot_command`'s existing Python-callback pattern.

## Real F9477 measurements: mode-3 + fragment-intensity + mokapot (2026-08-31)

Full 2x2x3 ablation grid (`jobs/f9477_ablation/`: RT+IIM on/off x
fragment-intensity on/off x mokapot {off, default, xgboost}), all against the
same F9477 raw data, sharing one node store (IM2Deep/Chronologer/the
fragment-intensity cache reused across every combo, no recompute). Real
1%-FDR "ions" (distinct `(peptide, charge)`, `searchops`'s own definition):

| RT-IIM | intensity | mokapot | ions |
|---|---|---|---|
| off | off | off (SAGE own) | 21,873 |
| off | on | off (SAGE own) | 21,873 |
| on | off | off (SAGE own) | 22,500 |
| on | on | off (SAGE own) | 22,500 |
| on | on | default | 31,710 |
| **on** | **on** | **xgboost** | **32,088** |

SAGE's own native ranking/eviction is completely blind to the
fragment-intensity features (identical SAGE-own numbers with intensity on vs
off, both RT-IIM settings) — the entire intensity gain flows through
mokapot's learned model, never SAGE's own score. Best combination found:
RT+IIM + fragment-intensity + mokapot(xgboost), 32,088 ions — beats every
previously-recorded number for this dataset (old `ranking_score`-focused best
was 25,409 ions, RT-only, no mokapot, no real fragment-intensity signal —
see `software/sage/devel_fixed`'s `CLAUDE.md`).

**Follow-up finding**: a later RT-only run (fragment-intensity + mokapot
xgboost, no IIM at all) beat the RT+IIM number above — 34,133 ions vs.
32,088. IIM was actively hurting, not just failing to help: SAGE's own
ions dropped 24,642→22,500 when IIM was added on top of RT, consistent
with IM2Deep's IIM eviction quality being weaker than SAGE's own internal
model (a "roughly a wash" comparison noted elsewhere) — the `combined_score`
ranking penalty (`0.5·(z_rt²+z_iim²)`) adds noise from a mediocre dimension
without adding signal. A follow-up `[tof_score_filter]` ablation on top of
this RT-only best (`score_margin=0.05`, this pipeline's existing TOF-neighbor
score-competition filter) came out roughly neutral (34,020 vs 34,133 ions,
-0.3%) — not worth keeping for this config.

## `predict_rt`/`predict_iim` finally get real caching (2026-09-01)

Both rules gained a new required `RtPredictionCache`/`IimPredictionCache`
input (mutable NodeTypes, mirroring `FragmentIntensityCache`), produced by
new `fill_rt_prediction_cache`/`fill_iim_prediction_cache` rules and passed
through to `git/featureprediction`'s new `--cache-path` flag. Before this,
`predict_rt`/`predict_iim` had no way to pass a cache path at all — every
job made a full, uncached Koina call for every sequence, even when overlap
with a previous job's `dumped_peptides` was total. Split into two rules
(fill + predict) because necroflow's `mutable=True` requires a
single-output rule; `predict_rt`/`predict_iim` each already have three
outputs. See `git/featureprediction`'s AI.md for the CLI-side change and
the real `mmappeteer` bulk-scale bug this surfaced (fixed separately).

**Two new `_pos_inputs` gotchas hit while wiring this** (both apply
generally to any future mixed Node/`None` or scalar rule input, not just
this change):

1. Node-shaped rule inputs (matched via `_node_input_contract`) must be
   passed *positionally*, never as `keyword=` — necroflow's
   `Rule._validate_input_presence` raises `TypeError: <rule>: unexpected
   inputs: [...]` if a Node-typed param is passed by keyword, since it
   classifies every input as strictly `_pos_inputs` (Node-shaped) or
   `_kw_inputs` (plain scalar) at rule-declaration time, never both. Bit
   `run_sage`'s `predicted_fragment_intensity_index`/`_cache` first (see
   the section above) and would have bitten `predict_rt`/`predict_iim`'s
   new cache params the same way had they been added as keywords.
2. Plain scalar (`_kw_inputs`) rule params can never be passed a literal
   `None` at a call site, even when the function's own Python default is
   `None` — necroflow records the selected value in `dependencies.toml`
   for provenance (tomlkit-serialized), and tomlkit has no representation
   for Python's `None`/TOML's absent-value. Use `""` (or another
   falsy-but-serializable sentinel matching the command template's own
   truthiness check) instead. Found via `[mokapot].plugin` (see above);
   the same trap applies to any future optional string/scalar rule param.

Real, measured payoff (F9477, RT-only jobs sharing `dumped_peptides`):
first job pays a real fill (255.3s, 6.3M sequences, plus SQLite bulk-insert
overhead — *slower* than the old always-uncached 92.1s single call, since
caching adds real DB write cost on top of the same network call). A
second, independent job (different `[recalibration]`/`[fragment_intensity]`
config, same `dumped_peptides`) then gets `fill_rt_prediction_cache` *and*
`predict_rt` itself fully reused (confirmed via `./nf -n`, neither node
appears in the "would-run" list) — the entire RT-prediction cost becomes
zero for that job, not just the Koina call. The break-even is the second
job, not the first.
