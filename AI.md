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

## Recalibration: three independently selectable modes (B.6, 2026-08-20)

`ionmaiden_pipeline`'s `if "sage" in cfg:` branch has three modes, gated by
two nested config keys so each is independently selectable for comparison
testing — see `plans/better_sage_filtering.md`'s B.6 for the full design:

1. **No recalibration** — no `[recalibration]` section at all. One SAGE
   pass, `run_sage` directly on `search_mz_pmsms`/`search_precursors`.
2. **mz recalibration alone** — `[recalibration]` present, no
   `[recalibration.rt_iim]`. Existing two-pass mz correction
   (`recalibrate_pmsms_mz`/`recalibrate_precursors`/`update_sage_config`),
   unchanged since before B.6.
3. **mz + RT + IIM recalibration** — `[recalibration.rt_iim]` also present.
   Nests *inside* mode 2's branch (chains onto mode 2's already
   mz-corrected outputs, doesn't duplicate them): `predict_rt_iim` (the
   `--predicted-properties` cache, `git/featureprediction`) runs off the
   same `filtered_sage_results_tsv` anchors the mz fits use, then
   `correct_precursors_rt_iim` (`git/featureprediction`'s B.6 precursor
   correction) chains onto `recalibrated_precursors`, then
   `update_sage_config_rt_iim` chains onto `update_sage_config`'s output to
   add `rt_tol_sec`/`mobility_tol`, then `run_sage_with_predicted_properties`
   (a distinct rule from plain `run_sage`, adding `--predicted-properties`
   — necroflow `@command` templates are static, can't conditionally add a
   flag, and `predictions` would need to become an optional DAG input
   either way; see `_mokapot_command` for the Python-callback alternative,
   not used here).

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
