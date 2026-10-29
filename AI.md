# AI.md — ionmaidentools

## What this repo is

`src/ionmaidentools/pipelines.py` holds the necroflow workflow for the ionmaiden
DIA timsTOF pipeline — the module almost every job TOML's `.pipeline` key points at
(e.g. `git/ionmaidentools/src/ionmaidentools/pipelines.py:ionmaiden_pipeline`): one large
module, `Node` subclasses + `@command`-decorated rule functions + the `@workflow`.
`src/ionmaidentools/fragpipe_pipeline.py` builds on it with the FragPipe types, rules and
two FragPipe workflows; `pipelines.py` imports nothing from it
(`plans/fragpipe_pipeline_split.md` in necromerge2: an import the other way would load a
second copy of `pipelines.py`, since necroflow loads a `.pipeline` file as an
unregistered module, and moving the shared types out would rehash most rules).

Package metadata (`pyproject.toml`): depends on `tomlkit` + `dictodot`, no CLI entry
points registered — this package is imported by necroflow's job runner (`./nf`, i.e.
`.venv/bin/necroflow`), not invoked directly.

## Necroflow rule conventions (as used here)

- A rule is a plain function decorated `@command`. Its body declares the outputs
  (`a, b = outputs(A, B)`), then the shell template `command = "... {placeholders}"`, then
  returns the outputs. When flags depend on which inputs are given, `command` is a nested
  `def command():` returning the template (`run_sage`, `update_sage_config`,
  `update_precursors`, `mokapot`).
  Parameters annotated with a `Node` subclass are DAG-edge inputs (parent nodes);
  plain-typed parameters (`str`, `int`, `float`) are scalar config/CLI-flag inputs, not
  DAG edges; `threads: Resource[int] = CORES` reserves cores, outside node identity.
- Workflows (`ionmaiden_pipeline(P, config)`, `fragpipe_synthetic_pipeline(P, config)`)
  are decorated `@workflow`, so `P` is the active Pipeline while one runs and rules are
  called without it: `P.some_label = some_rule(parent_node, ..., scalar_kwarg=...)`.
  Helpers defined inside a workflow (e.g. `_finalize_confident_psms`) use the same
  active Pipeline. Necroflow still accepts an explicit leading `P`; this module no
  longer passes it.
- Both workflows serialize their complete Necroflow-supplied config through
  `write_pipeline_config` and publish it as requestable `P.pipeline_config`
  (`pipeline_config.toml`). Necroflow strips metadata keys such as `.pipeline` and
  `.requests`, resolves `.extends`, and expands grids before calling the workflow, so
  the file records the effective workflow input.
- Every rule shells out to a real installed/patched CLI directly — no bespoke Python
  wrapper CLIs. External tools live at fixed, pre-installed paths (`git/sage/target/release`,
  `software/fragpipe/fragpipe-24.0`), same convention on both sides.
- See root `CLAUDE.md`'s "Precursor table format" — all precursor tables in this pipeline
  are `.mmappet` directories (`MmappetDataset` subclasses here), never intermediate parquet.

## Three workflows

- `pipelines.py:ionmaiden_pipeline`: the real pipeline — Bruker `.d` → MS1 peak picking →
  quadrupole transmission → pseudo-MS/MS → SAGE search (+ optional recalibration pass,
  optional MGF export when the job config has an `[mgf]` section). `[mgf].config_path` is
  passed as a scalar input to `convert_search_pmsms_to_mgf`; its value therefore selects
  and fingerprints the MGF flavour. Jobs without that key do not add `search_mgf` to the
  pipeline. A `[fragpipe]` section is an error here. Its body is
  `build_ionmaiden_pipeline(P, config)`, which the FragPipe workflow reuses.
- `fragpipe_pipeline.py:ionmaiden_fragpipe_pipeline`: `ionmaiden_pipeline` plus a FragPipe
  side-by-side search of its search mzML (`[fragpipe]` required: `workflow_path`,
  optional `ram`).
- `fragpipe_pipeline.py:fragpipe_synthetic_pipeline`: FragPipe smoke test on
  Koina-simulated peptides from a FASTA (no Bruker `.d` input, no Sage).

Both FragPipe workflows run the same search, `fragpipe_search` (`source_fragpipe_workflow`/
`generate_fragpipe_decoy_fasta`/`patch_fragpipe_workflow`/`write_fragpipe_manifest`/
`run_fragpipe`/`extract_fragpipe_log`/`summarize_fragpipe`, labelled `fragpipe_*`).

## Documentation index

This file stays a short overview; design rationale, migration history, and real
measurements for each feature live in `docs/ai/`, one file per topic:

| File | Covers |
|------|--------|
| `docs/ai/ms1_ms2_extraction.md` | C++ MS1 extraction, compact/query-oriented MS2 benchmark targets |
| `docs/ai/fragment_intensity.md` | Fragment-intensity `PredictionCache`, its export and `run_sage` wiring -- removed 2026-10-26; SPECLIB (`[speclib]`, `check_speclib`) is the only source now |
| `docs/ai/recalibration_modes.md` | The three mz/RT/IIM recalibration modes, table-presence gating, tolerance percentiles/method, fragment m/z correction applied inside SAGE (`fragment_shift_ppm`) |
| `docs/ai/run_sage_merge.md` | `run_sage`/`run_sage_with_predicted` merged into one mixed-Node/value rule |
| `docs/ai/mokapot_integration.md` | Config-driven mokapot plugin, leakage-safe PIN filtering, real F9477 ablation-grid measurements |
| `docs/ai/fdr_levels.md` | Spectra/peptidoforms/ions at 1% FDR as the last step of `run_sage`/`mokapot`, and recording runs in `$ION_MAIDEN_REGRESSION_DB` |
| `docs/ai/rt_iim_caching.md` | `predict_rt`/`predict_iim`'s `PredictionCache` wiring, necroflow gotchas, measured payoff -- removed 2026-10-27; RT/IIM come from the SPECLIB (`[recalibration.rt]`/`.iim]` require its parts) |
| `docs/ai/mapped_inputs.md` | `dump_compact_peptides`/`CompactPeptidesTable` for the sliced search, and the prediction nodes as mmappet datasets SAGE maps (2026-10-17) |
| `docs/ai/bestrun.md` | `jobs/f9477_best.toml` — the standing best-known-config job, current numbers, and the ask-before-updating rule |

Check freshness against `git log` on the file(s) a `docs/ai/*.md` entry names if something
looks stale — there is no automated freshness check for this repo's own docs (unlike the
top-level monorepo's `summarise/fingerprints.json`).

## Conventions to follow when editing this file

- Design rationale, migration history, and "why this and not that" reasoning belongs
  in this file or `docs/ai/*.md` (never in `pipelines.py` docstrings — keep code
  docstrings about what the code does now, not why it changed).
- Before appending a substantial new section to this root file, check whether it
  belongs in an existing `docs/ai/*.md` topic instead, or needs a new one — keep this
  file to overview + index, not accumulating topics again.
- FragPipe's `.workflow` file is a hand-maintained flat `key=value` properties file
  (`configs/search/fragpipe/workflows/*.workflow`), symlinked in as-is by
  `source_fragpipe_workflow`. `patch_fragpipe_workflow` derives a copy with
  `database.db-path=` rewritten (via `sed`) to point at `generate_fragpipe_decoy_fasta`'s
  output, so `run_fragpipe` always searches against a decoy-augmented FASTA — this is
  the one exception to "never patched", added because FragPipe/Philosopher has no
  generate-decoys-at-search-time option (unlike Sage) and there's no CLI override flag
  for it.
