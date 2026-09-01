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

## Documentation index

This file stays a short overview; design rationale, migration history, and real
measurements for each feature live in `docs/ai/`, one file per topic:

| File | Covers |
|------|--------|
| `docs/ai/ms1_ms2_extraction.md` | C++ MS1 extraction, compact/query-oriented MS2 benchmark targets |
| `docs/ai/fragment_intensity.md` | Fragment-intensity `PredictionCache` (`mutable=True`), its export, and wiring into `run_sage` |
| `docs/ai/recalibration_modes.md` | The three mz/RT/IIM recalibration modes, table-presence gating, tolerance percentiles/method, `server_url` |
| `docs/ai/run_sage_merge.md` | `run_sage`/`run_sage_with_predicted` merged into one mixed-Node/value rule |
| `docs/ai/mokapot_integration.md` | Config-driven mokapot plugin, leakage-safe PIN filtering, real F9477 ablation-grid measurements |
| `docs/ai/rt_iim_caching.md` | `predict_rt`/`predict_iim`'s `PredictionCache` wiring, necroflow gotchas, measured payoff |
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
