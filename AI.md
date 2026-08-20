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
- Every rule shells out to a real installed/patched CLI directly — no bespoke Python
  wrapper CLIs. External tools live at fixed, pre-installed paths (`software/sage/devel_fixed`,
  `software/fragpipe/fragpipe-24.0`), same convention on both sides.
- See root `CLAUDE.md`'s "Precursor table format" — all precursor tables in this pipeline
  are `.mmappet` directories (`MmappetDataset` subclasses here), never intermediate parquet.

## Two pipeline factories

- `ionmaiden_pipeline`: the real pipeline — Bruker `.d` → MS1 peak picking → quadrupole
  transmission → pseudo-MS/MS → SAGE search (+ optional recalibration pass, optional
  FragPipe side-by-side comparison when the job config has a `[fragpipe]` section).
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

**Config derivation, not restatement** (see `git/featureprediction`'s
`AI.md` for the numbers): `[recalibration.rt_iim]`'s `min_charge`/
`max_charge` default to `cfg.sage.get("precursor_charge", (2, 4))` (SAGE's
own compiled-in default) rather than crashing when `sage.precursor_charge`
is unset, which real job configs often do. `tolerance_percentiles` is
required explicitly in `[recalibration.rt_iim]` — deliberately *not*
inherited from `cfg.recalibration`'s own (mz-scoped) value, since real
values there (e.g. `[0, 100]`, no trimming) would be a bad default for RT
specifically.

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
