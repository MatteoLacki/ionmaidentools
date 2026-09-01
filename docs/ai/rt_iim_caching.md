# `predict_rt`/`predict_iim` finally get real caching (2026-09-01)

Both rules gained a new required `RtPredictionCache`/`IimPredictionCache`
input (mutable NodeTypes, mirroring `FragmentIntensityCache` — see
`fragment_intensity.md`), produced by new `fill_rt_prediction_cache`/
`fill_iim_prediction_cache` rules and passed through to
`git/featureprediction`'s new `--cache-path` flag. Before this,
`predict_rt`/`predict_iim` had no way to pass a cache path at all — every
job made a full, uncached Koina call for every sequence, even when overlap
with a previous job's `dumped_peptides` was total. Split into two rules
(fill + predict) because necroflow's `mutable=True` requires a
single-output rule; `predict_rt`/`predict_iim` each already have three
outputs. See `git/featureprediction`'s AI.md for the CLI-side change and
the real `mmappeteer` bulk-scale bug this surfaced (fixed separately).

## Two `_pos_inputs` gotchas hit while wiring this

Both apply generally to any future mixed Node/`None` or scalar rule input,
not just this change:

1. **Fixed upstream 2026-09-01** (necroflow's "Accept Node inputs by
   keyword, not just positional") — at the time this was found, Node-shaped
   rule inputs (matched via `_node_input_contract`) had to be passed
   *positionally*, never as `keyword=`: `Rule._validate_input_presence`
   raised `TypeError: <rule>: unexpected inputs: [...]` if a Node-typed
   param was passed by keyword, since it classified every input as strictly
   `_pos_inputs` (Node-shaped) or `_kw_inputs` (plain scalar) at
   rule-declaration time, never both. Bit `run_sage`'s
   `predicted_fragment_intensity_index`/`_cache` first (see
   `fragment_intensity.md`) and would have bitten `predict_rt`/
   `predict_iim`'s new cache params the same way had they been added as
   keywords. After pulling that necroflow release, `run_sage`'s call sites
   were cleaned up to pass `predicted_rt`/`predicted_iim`/the
   fragment-intensity pair by keyword (dropping the `None, None`
   positional-placeholder pattern the mode-1/mode-2 call sites needed
   before) — verified call-syntax changes alone don't affect provenance
   (necroflow's own commit message guarantees this; confirmed empirically
   too, zero hash changes across mode-1/2/3 dry-runs before vs. after).
2. Plain scalar (`_kw_inputs`) rule params can never be passed a literal
   `None` at a call site, even when the function's own Python default is
   `None` — necroflow records the selected value in `dependencies.toml`
   for provenance (tomlkit-serialized), and tomlkit has no representation
   for Python's `None`/TOML's absent-value. Use `""` (or another
   falsy-but-serializable sentinel matching the command template's own
   truthiness check) instead. Found via `[mokapot].plugin` (see
   `mokapot_integration.md`); the same trap applies to any future optional
   string/scalar rule param.

## Real, measured payoff

F9477, RT-only jobs sharing `dumped_peptides`: first job pays a real fill
(255.3s, 6.3M sequences, plus SQLite bulk-insert overhead — *slower* than
the old always-uncached 92.1s single call, since caching adds real DB
write cost on top of the same network call). A second, independent job
(different `[recalibration]`/`[fragment_intensity]` config, same
`dumped_peptides`) then gets `fill_rt_prediction_cache` *and* `predict_rt`
itself fully reused (confirmed via `./nf -n`, neither node appears in the
"would-run" list) — the entire RT-prediction cost becomes zero for that
job, not just the Koina call. The break-even is the second job, not the
first.
