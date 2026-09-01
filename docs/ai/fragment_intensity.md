# Fragment-intensity prediction cache

## This repo's first `mutable=True` rule (2026-08-31)

`predict_fragment_intensity` wraps `git/featureprediction`'s
`feature-prediction-generate-fragments` (see that repo's `AI.md` for the
Prosit/`compact_trt` predictor itself and the cache's own internal shape).
Independently requestable (same pattern as `ms2_tsf_events`/`ms2_tfs_events`,
see `ms1_ms2_extraction.md`) off `dumped_peptides` — never runs just because a
job dumps peptides, since a full human-proteome fill is a real, expensive,
network-bound operation (~43 minutes against the live Koina server measured
on F9477).

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

See `plans/fragment_intensity_cache.md` for the cache's own design (predictor
choice, slot-numbering scheme, dtype choices).

## `export_fragment_intensity_for_sage`: ordinary rule consuming a mutable parent (2026-08-31)

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
consumer at the time this rule was added, never runs just because
`dumped_peptides` exists.

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

## Wired into `run_sage`, `[fragment_intensity]` gate (2026-08-31)

`export_fragment_intensity_for_sage`/`predict_fragment_intensity` had no
downstream consumer until this point — the Rust-side reader
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
convention as `[recalibration.rt]`/`[recalibration.iim]` (see
`recalibration_modes.md`).

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
plain scalar (`_kw_inputs`) config values can't safely be `None` (see
`mokapot_integration.md`'s `[mokapot].plugin` section, and
`rt_iim_caching.md`'s gotcha list, for the general rule).

See `mokapot_integration.md` for the real F9477 measurements of how much
fragment-intensity signal is worth once it reaches mokapot.
