# Mokapot integration

## Forked mokapot (`git/mokapot`), real CLI restored, model persistence dropped (2026-09-17)

Superseded most of the section below same-day. `scripts/run_mokapot.py` (a
Python-API reimplementation of the `mokapot` CLI, built to bypass the CLI's
slow PIN parser) is deleted -- the parser is now fixed at the source instead:
`git/mokapot` is a real fork (`github.com/MatteoLacki/mokapot`, pinned to a
`v0.10.0`-based branch `mokapot-necroflow-integration` -- **not** `main`,
which has since moved to a completely different `OnDiskPsmDataset`/
`TabularDataReader` architecture incompatible with both this fork's patches
and `mokapot-xgboost-plugin`'s `Model`/`brew()`-based API). Two patches:

- `read_percolator` now uses DuckDB's CSV reader instead of a hand-rolled
  line-by-line `str.split` + `pd.DataFrame.from_records` loop -- roughly 2x
  faster on a real PIN file, verified to produce identical row/column
  counts, dtypes, and values against the original.
- `--seed` now actually seeds `brew()`'s own rng (`rng=config.seed` threaded
  into `main()`'s `brew()` call) -- previously a no-op, since
  `np.random.seed(config.seed)` only touches the legacy global RNG that
  `brew()`'s `np.random.default_rng()`-based generator never reads. Verified:
  two independent real runs with the same `--seed` now produce byte-identical
  output files.

`Makefile`'s `venvs/mokapot/bin/mokapot` target installs this fork editable
(`pip install -e git/mokapot`) instead of plain `mokapot` from PyPI.

The real CLI is back in `_mokapot_command` (`venvs/mokapot/bin/mokapot`, not
`scripts/run_mokapot.py`) -- once the parser and seed were fixed at the
source, the CLI became strictly better than maintaining a parallel
reimplementation: less surface to keep in sync (the `--keep_decoys` default
mismatch bug from the `run_mokapot.py` era couldn't have happened here), and
any future upstream mokapot fix applies automatically. `mokapot()` is back to
its original 3-output shape (`used_pin, peptides, psms`) -- the
`MokapotModelFold1/2/3` NodeTypes and per-fold `Model.save()` calls from the
`run_mokapot.py` era are removed; not needed in the pipeline for now.

**Seeding, job-config-driven, for every random-state site in this part of
the pipeline**: new `[mokapot].seed` / `[sagepy_rescore].seed` key (default
1, matching mokapot's own CLI default -- omitting it changes nothing for
existing jobs), reused as `--xgboost_seed`/`--lightgbm_seed` too. Three
independent rng sites, now all covered: `brew()`'s own rng (mokapot core,
fixed above), and each plugin's estimator `random_state`
(`mokapot-xgboost-plugin`'s `XGBClassifier`/`BaggingClassifier`,
`mokapot-lightgbm-plugin`'s `LGBMClassifier`) -- previously none of the three
estimator classes set `random_state` at all.

**`[mokapot].plugin = "lightgbm"`** is now a real, separate mokapot plugin
(`git/mokapot-lightgbm-plugin`, sibling package to `mokapot-xgboost-plugin`,
same `BasePlugin`/entry-point pattern) instead of a `--model` flag inside
`run_mokapot.py` -- matches mokapot's own native multi-plugin mechanism.
Same GPU self-detection pattern as documented below (one-time smoke-test
fit, `LightGBMError` -> CPU fallback), now living in the plugin's own
`_select_device()` instead of a driver script.

Real F9477 validation after the full rewrite, `jobs/f9477_best.toml -call`,
only `mokapot` (and its immediate upstream, re-triggered by unrelated
earlier activity in the node store) re-ran:

| | wall-clock | PSMs (q≤0.01) | peptides | ions |
|---|---:|---:|---:|---:|
| xgboost (via real CLI, patched fork) | 46.0s | 99,147 | 31,802 | 34,486 |
| lightgbm (direct CLI test, same PIN) | 48.2s | 100,931 | 31,999 | 34,682 |

Both land inside the day's established noise band -- consistent with the
`scripts/run_mokapot.py`-era numbers below, confirming the rewrite changed
nothing about model behavior, only where the code lives.

## Model persistence, `[mokapot].model` choice, and GPU (2026-09-17)

**Superseded by the section above, same day** -- kept here as the record of
what was tried and why it changed, not current behavior.

**Trained fold models are no longer discarded.** `scripts/run_mokapot.py`
always calls `Model.save()` on each of the 3 fold models after `brew()`
returns (same `mokapot.model_fold-N.pkl` naming/pickle format mokapot's
own CLI uses for `--save_models`, never gated behind a flag here). Wired
into necroflow as 3 new fixed `NodeType`s (`MokapotModelFold1/2/3` --
fixed count because `_MOKAPOT_FOLDS` is fixed at 3), so `mokapot()` now
returns a 6-tuple (`used_pin, peptides, psms, model_fold_1, model_fold_2,
model_fold_3`) instead of 3 -- both call sites updated. Real reasons this
matters: inspecting a fitted model (feature importances, iteration count)
no longer requires re-running the whole `brew()`; mokapot's own
`--load_models`/passing pretrained `Model` objects into `brew()` lets a
saved model rescore a different PIN without retraining; and since
mokapot's `--seed` doesn't actually seed `brew()`'s own rng (see below),
the saved model is the only record of which model produced a given run's
numbers.

**`[mokapot].model`** (`""`/`"xgboost"`/`"histgb"`/`"lightgbm"`, only
meaningful together with `plugin = "xgboost"`) selects which estimator
`run_mokapot.py` builds, threaded through the same "" (not `None`)
sentinel convention as `plugin` itself. `histgb` exists in the script for
ad-hoc benchmarking but isn't exposed as a first-class job-config choice
(measured worse, see below) -- `xgboost`/`lightgbm` are the two real
options.

**GPU.** XGBoost's plugin (`mokapot_xgboost/plugin.py`) now passes
`device="cuda"` unconditionally. Confirmed safe on a real GPU-less
machine (AMD Renoir iGPU, no NVIDIA hardware, no `/dev/nvidia*`, no
`libcuda.so` anywhere) at realistic scale (180k rows x 36 features): a
`UserWarning` ("No visible GPU is found, setting device to CPU") and a
normal CPU fit, 2.29s for 200 estimators -- no exception, no silent
breakage. The same pip wheel bundles CUDA support and self-detects; no
separate install needed on a GPU machine.

LightGBM has no equivalent self-detection -- its standard pip wheel has
no GPU code compiled in at all (`device_type="cuda"` raises
`LightGBMError: CUDA Tree Learner was not enabled in this build`), and
even a GPU-enabled build raises rather than falling back on a
device-unavailable error. Two pieces compensate:
- `Makefile`'s `venvs/mokapot/bin/mokapot` target now branches on
  `command -v nvidia-smi` / `/dev/nvidia0` at `make` time: NVIDIA present
  -> `pip install lightgbm --config-settings=cmake.define.USE_CUDA=ON`
  (needs the CUDA toolkit/`nvcc` on that machine at install time, not
  just a driver); absent -> the plain CPU wheel (this dev machine's
  actual path, unchanged).
- `run_mokapot.py::_select_lightgbm_device` replaces the self-detection
  LightGBM doesn't have: one cheap smoke-test fit on trivial synthetic
  data with `device_type="cuda"`, catching `LightGBMError` to decide
  `cuda` vs `cpu` for the whole run. Decided once upfront, not per-fit,
  because mokapot's own bootstrap loop calls `.fit()` 30 times internally
  (3 folds x 10 iterations) with no hook to catch/retry mid-loop, and
  LightGBM's device-unavailable failures happen at device-init time
  (near-instant), not partway through a real fit, so nothing costly is
  ever discarded by probing first.

**AMD iGPU (this dev machine's Renoir)**: neither library can use it.
XGBoost's GPU path is CUDA-only, no AMD/ROCm backend exists at all.
LightGBM's OpenCL `device_type="gpu"` path is vendor-neutral in
principle, but `clinfo` reports zero OpenCL platforms on this machine (no
AMD OpenCL/ROCm runtime configured) -- not pursued; iGPU OpenCL compute
for gradient boosting is generally unreliable/low-value even where
available. The real GPU target is a separate server with an NVIDIA RTX
A4000 (16GB VRAM) -- untested from here (no CUDA toolkit, no GPU on this
machine to verify actual GPU execution or the CUDA-enabled LightGBM
build); only the CPU-fallback path is verified end-to-end.

**Real F9477 model-family benchmark (2026-09-17)**, same `used.pin`
(`nodes/mokapot/ac305abb...`), same `train_fdr=0.05`/`test_fdr=0.01`,
CPU-only (this machine), `xgboost_n_jobs=5`/`max_workers=3` where
applicable:

| model | wall-clock | ions |
|---|---:|---:|
| plain xgboost | ~44-46s | 33,636-34,764 |
| lightgbm | 44.2s | 34,665 |
| histgb | 81.5s | 34,450 |
| bagged xgboost (5x30k, see below) | ~130-140s | 32,828-34,212 |

LightGBM ties plain XGBoost on both wall-clock and ions (both inside
xgboost's own noise band) -- the only alternative tried that isn't
strictly worse. HistGradientBoostingClassifier is ~1.8x slower (no
`n_jobs`-equivalent thread control). Neither beats plain XGBoost outright
on CPU; the value of the choice is future GPU headroom, not a CPU-side
win today.

## Bagged XGBoost ensemble: `xgboost_bagging_*` (2026-09-17)

Tried per a direct request to explore ensembling as an alternative to a
single XGBoost fit: `sklearn.ensemble.BaggingClassifier(estimator=
XGBClassifier(...), n_estimators=5, max_samples=30_000, bootstrap=False,
random_state=<seed>)` -- 5 copies, each fit on a random 30k-row subsample
without replacement ("subagging", not bootstrap -- duplicated rows would
just make each copy as expensive as the unbagged fit). `predict_proba`
averages the 5 copies' probabilities (soft voting) automatically; this
matters because mokapot's `_get_scores` needs a continuous score for its
q-value ranking, and a hard majority-class vote would collapse to
`{0..5}` distinct values and wreck that ranking's resolution -- soft
voting was not optional here.

Drops straight into `Model(estimator, train_fdr=...)` with no changes to
mokapot itself (`BaggingClassifier` has no `decision_function` any more
than bare `XGBClassifier` does, so the existing `predict_proba` fallback
in `Model._get_scores` covers it unchanged). New CLI flags on the plugin:
`--xgboost_bagging_n_estimators` (0 default = off, plain `XGBClassifier`
exactly as before), `--xgboost_bagging_max_samples` (default 30,000),
`--xgboost_bagging_seed` (`random_state` -- the one part of this whole
setup that's actually reproducible, unlike mokapot's own unseeded
`brew()` rng).

**Real F9477 measurements**: ~3x slower than plain XGBoost (~130-140s vs.
~44-46s) -- each of the 30 bootstrap-loop calls (3 folds x 10 iterations)
now does 5 separate XGBoost fits serially (deliberately not parallelized
here, to avoid a 3rd nested-parallelism axis on top of fold-level
`--max_workers` and per-fit `--xgboost_n_jobs`), and per-fit fixed
overhead (DMatrix construction, booster init) doesn't shrink
proportionally with fewer rows. Tried reducing `max_samples` 10x (30k ->
3k) expecting a speed win: **made it worse on both axes** -- 140.9s
(slower, not faster: smaller samples don't reduce the fixed per-fit
overhead that dominates) and 32,828 ions (a real quality drop, clearly
outside the noise band -- 3,000 rows is too thin a sample for a
36-feature model). 30k confirmed as the better setting between the two
tried. Ion counts at 30k (32,828-34,212 across two seeds) sit inside or
just below plain XGBoost's noise band -- no demonstrated quality
improvement from bagging on this dataset, only cost.

## Auto-derived `--max_workers`/`--xgboost_n_jobs` threading (2026-09-17)

`mokapot.brew()`'s 3-fold CV runs strictly serial by default (`--max_workers`
CLI default 1, never overridden here), and the xgboost plugin forced
`n_jobs=1` inside `XGBClassifier` -- so a mokapot run was 3 folds x 10
`Model.fit` iterations = 30 full single-threaded XGBoost fits, one at a
time, on a 16-core/62GB box. That's the actual bottleneck behind mokapot
being one of the slower pipeline stages: not XGBoost's per-fit cost, but
fold- and thread-level serialization neither mokapot's defaults nor this
pipeline's command ever turned on.

`_mokapot_command` (`pipelines.py`) now computes core split at command-build
time: `max_workers = min(_MOKAPOT_FOLDS=3, n_cores)`,
`xgboost_n_jobs = max(1, n_cores // max_workers)` (`os.sched_getaffinity(0)`,
not `os.cpu_count()`, so a cgroup/taskset-constrained container -- see
`docker_present.sh` -- doesn't get oversubscribed). On this repo's dev box
(16 cores): 3 folds in parallel x ~5 threads/fit each, 15/16 cores used, one
core headroom. Applies to both mokapot call sites (plain-SAGE-PIN and the
sagepy_rescore branch) uniformly, plugin or no plugin -- `--max_workers`
helps the default linear-SVM path too; `--xgboost_n_jobs` only emitted when
`plugin == "xgboost"` since it's a plugin-specific flag.

Deliberately **not** a `mokapot()` rule kwarg / job-TOML knob -- it's pure
execution tuning (thread count), not model behavior, so it doesn't belong in
declared rule config. Confirmed safe against necroflow's fingerprinting
before wiring this in: `RuleCall.compute_provenance_hash` hashes declared
`config`/parent lineage/execution context, never the realized command
string (`git/necroflow/src/necroflow/rule_call.py`), so baking a
core-count-derived flag straight into the rendered shell command cannot
make node identity depend on which machine executes it -- same class of
hazard `mscaches` was designed to avoid, checked rather than assumed.

`git/mokapot-xgboost-plugin`'s `XGBoostPlugin.add_arguments` gained
`--xgboost_n_jobs` (default 1, unchanged for anyone invoking the plugin
outside this pipeline). Its sibling hardcoded `n_jobs=1` in
`git/sagepy-rescore/src/sagepy_rescore/mokapot_runner.py::_build_model`
(the Python-API path, not used by necroflow's DAG -- see "config-driven
model choice" below) was intentionally left untouched: out of scope for a
"current pipeline" optimization since nothing in the executed DAG calls it.

**Real F9477 measurement (2026-09-17)**, `jobs/f9477_best.toml`, `-call`,
same `SpecId`→`used.pin` charge-join counting method as elsewhere in this
file. Baseline node (`nodes/mokapot/0662eacc...`) had already been built
under the pre-optimization code the same day; only `mokapot` re-ran for the
optimized node (`nodes/mokapot/ac305abb...`) -- all 43 upstream nodes hit
cache, confirming this change doesn't touch node identity for anything
except the `mokapot` rule itself:

| | wall-clock | PSMs (q≤0.01) | peptides | ions |
|---|---:|---:|---:|---:|
| before (serial, n_jobs=1) | 161.8s | 100,545 | 31,940 | 34,693 |
| after (max_workers=3, xgboost_n_jobs=5) | **54.8s** | 100,722 | 32,014 | **34,764** |
| Δ | **-107.0s (-66%, 2.95x)** | +177 | +74 | +71 (+0.2%) |

Ion delta is an order of magnitude inside the ~500-ion/1.5% mokapot
run-to-run noise floor -- a pure speedup, no measurable effect on search
quality.

### Post-optimization profile: what's left after threading (2026-09-17)

Wall-clock instrumentation (not cProfile -- its per-call overhead would
distort the very Python-vs-native ratio being measured) against the same
real `used.pin`, production settings (`xgboost_n_jobs=5`), `max_workers`
forced to 1 only so the 3 folds run serially for a clean single timeline
(67.8s here vs. 54.8s parallel -- same total work, not overlapped):

| bucket | time | % of wall | calls |
|---|---:|---:|---:|
| `XGBClassifier.fit` | 27.2s | 40.1% | 30 (3 folds x 10 iters) |
| `PsmDataset._update_labels` | 16.2s | 23.9% | 335 |
| `PsmDataset._find_best_feature` | 14.9s | 22.0% | 4 |
| `read_percolator` (mokapot's own PIN parser) | 8.9s | 13.2% | 1 |
| `_get_scores` (rescoring after each fit) | 4.9s | 7.2% | 33 |
| `assign_confidence` (final q-value/PEP) | 4.7s | 6.9% | 1 |
| `StandardScaler` fit/transform + `_split` | 1.4s | 2.0% | -- |

`_find_best_feature` calls `_update_labels` internally (once per feature x
direction, swept over all ~36 features, 4 total invocations -- once per
fold's initial-direction pick, once at brew's end for the
best-feature-vs-learned-model check) -- those two rows' time overlaps, not
additive.

**XGBoost fit is the single largest bucket (~40% of wall, ~48% of
brew-only time) but does not dominate outright.** Mokapot's own
pure-Python/pandas bookkeeping -- relabeling targets/decoys against
`train_fdr` every iteration, plus the 4x initial best-feature sweep over
the full PSM x feature matrix -- is comparable in total size. mokapot's
hand-rolled PIN parser (`mokapot/parsers/pin.py::_parse_in_chunks`: manual
`str.split` per line + `pd.DataFrame.from_records`, not pandas' C
`read_csv`) costs a flat ~9s/13% independent of any thread tuning. None of
these three are reachable by `--max_workers`/`--xgboost_n_jobs` -- they're
inside mokapot's own library internals. No GPU on the dev box
(`nvidia-smi` not found), so CUDA XGBoost isn't an available lever either.
Going further here means patching mokapot's own parser/label-update code,
not this pipeline's rule -- not attempted.

### `scripts/run_mokapot.py`: bypass the CLI's PIN parser (2026-09-17)

Acted on the `pin_parse` finding above. `mokapot.parsers.pin.read_pin`
accepts a `pandas.DataFrame` directly (`pin.py:98-103`) and skips its own
`read_percolator` -- a hand-rolled per-line `str.split` + `pd.DataFrame
.from_records`, not pandas' C `read_csv` -- entirely when given one.
Measured on the real F9477 `used.pin` (362,288 rows): `pd.read_csv` 1.53s
vs. mokapot's own parser 8.94s. The DataFrame branch's
`pin_files.copy(deep=copy_data)` (default `copy_data=False`, i.e.
shallow) costs 0.08ms -- free; it exists only so `read_pin`'s later
in-place Label-column dtype/recode mutation can't leak back into the
caller's DataFrame (verified: reassigning a column on a shallow copy
doesn't touch the original's dtype).

New `scripts/run_mokapot.py` (top-level, next to `mokapot_pin_adapter.py`
-- same precedent: a bare script `venvs/mokapot/bin/python` runs directly,
not a package) replaces the `mokapot` CLI binary as the second half of
`_mokapot_command`'s shell chain: `pd.read_csv(used_pin)` →
`mokapot.parsers.pin.read_pin(df)` → `mokapot.brew.brew(...)` →
`results.to_txt(...)`. For tab-delimited input mokapot's own CLI
(`get_parser`) calls bare `read_pin(config.psm_files)` with every optional
column-detection kwarg left at its function default (`config.py` never
wires `group_column`/`filename_column`/etc. from any CLI flag) -- so
`read_pin(pin_df)` here is exact parity, not an approximation. The
xgboost plugin (`mokapot_xgboost.plugin.XGBoostPlugin`) is imported and
called directly instead of through mokapot's `entry_points(group=
"mokapot.plugins")` discovery, since the caller already knows whether
xgboost is wanted.

**Considered and rejected: putting this in `searchops`.** `searchops`
already depends on `xgboost` (`models.py::XGBoostDerivativePenalizedModel`,
an unrelated fragment-m/z-recalibration spline model) but is only
installed in `venvs/common` -- not `venvs/mokapot`, where the real
`mokapot`/`xgboost`/`mokapot_xgboost` packages actually live. Installing
those into `venvs/common` too just to host this script there would
duplicate an environment for no reason; `sagepy-rescore`'s own
`mokapot_runner.py::rescore()` (a second, independent Python-API mokapot
runner, same read-df/brew/to_txt shape) has the identical problem --
installed in `venvs/sagepy_rescore`, wrong environment for this call site.
A bare top-level script matches this rule's actual execution environment
with zero new packaging.

**Bug caught during validation, fixed before wiring in:** first draft
hardcoded `results.to_txt(dest_dir=..., decoys=True)`, copied from
`sagepy_rescore/mokapot_runner.py`'s pattern without checking against the
real CLI. The CLI actually does `to_txt(decoys=config.keep_decoys)`, and
`--keep_decoys` defaults `False` (never passed by this pipeline) --
`decoys=True` silently added two extra files
(`mokapot.decoy.{peptides,psms}.txt`, ~26MB) that neither prior CLI-based
node output nor the declared `MokapotPeptides`/`MokapotPsms` node types
expect. Confirmed via direct diff against a real prior CLI-based node's
file listing before fixing to `decoys=False`. Row counts of the two real
output files were identical either way -- `decoys` only gates the extra
files, not main-file content.

**Real F9477 measurement (2026-09-17)**, `jobs/f9477_best.toml -call`,
same job/counting method, only `mokapot` re-ran each time (rule identity
changed each edit) -- mokapot-rule-alone wall-clock as necroflow itself
reports it (`mokapot in Xs`), not full-pipeline or standalone-script time:

| stage | wall-clock (mokapot alone) | Δ vs. previous | ions |
|---|---:|---:|---:|
| original (serial folds, n_jobs=1, CLI parser) | 161.8s | -- | 34,693 |
| + fold/thread parallelism (`--max_workers`/`--xgboost_n_jobs`) | 54.8s | -66% (2.95x) | 34,764 |
| + `run_mokapot.py` (bypasses CLI's PIN parser) | **46.0s** | -16% (1.19x) | 33,636 |
| **total** | | **-71.6% (3.52x)** | |

(Two additional standalone `run_mokapot.py` runs outside necroflow, used
to catch the `decoys=True` bug below before wiring it in, timed
45.6s/43.8s -- consistent with the 46.0s necroflow figure, not added to
the table above since they didn't go through the actual rule.)

Four real xgboost-plugin runs recorded today total (the three table rows
above plus one pre-fix standalone run) span **33,636-34,764 ions (~3.2%
spread)** -- wider
than the single ~1.5%/500-ion pair this file previously characterized the
noise floor by. Not a regression from this change: `mokapot`'s CLI
`--seed`/`np.random.seed(config.seed)` only seeds the legacy global RNG;
`brew()`'s own `rng = np.random.default_rng(rng)` is never passed a seed
by the CLI (`brew(datasets, model=model, test_fdr=..., folds=...,
max_workers=...)`, no `rng=` kwarg) or by `run_mokapot.py` either (same
call shape) -- train/test-split randomness was never actually
reproducible despite the CLI having a `--seed` flag, on either code path.
Treat **~3%, not ~1.5%**, as this dataset's real mokapot(xgboost)
run-to-run floor going forward.

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
flag omitted) — only the provenance-recording path differs. Same trap
applies generally to any future optional string/scalar rule param; see
`rt_iim_caching.md`'s gotcha list for the Node-shaped-input counterpart.

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

## sagepy_rescore branch removed (2026-09-18)

The whole second `mokapot(...)` call site (`SagepyRescoreConfig`/
`SagepyRescorePredictions`/`SagepyRescorePin` NodeTypes,
`write_sagepy_rescore_config`/`run_sagepy_rescore_predict`/
`write_sagepy_rescore_pin` rules, the `if "sagepy_rescore" in cfg:` block)
is gone — `git/sagepy`/`git/sagepy-rescore`/`git/sagepy_ionmaiden_adapter`
were dropped from the project entirely (nothing in the live pipeline used
sagepy, and a prior evaluation had already found sagepy-rescore's
xgboost/rbf-svm rescoring "not worth it — GPU needed, poor results"; see
`plans/sagepy-rescore.md`). This corrects two now-stale claims above: the
"sagepy_rescore branch's own separate call is untouched" note under
`[mokapot].plugin` no longer applies (there is no second call), and
`--mode passthrough` has no live caller today — the plain-SAGE-PIN call
site always passes `rt_source`/`iim_source`, i.e. always resolves to `--mode
sage`. `passthrough` is left in `scripts/mokapot_pin_adapter.py` as a
generic default, not removed.

## Real F9477 measurements: mode-3 + fragment-intensity + mokapot (2026-08-31)

Full 2x2x3 ablation grid (`jobs/f9477_ablation/`: RT+IIM on/off x
fragment-intensity on/off x mokapot {off, default, xgboost}), all against the
same F9477 raw data, sharing one node store (IM2Deep/Chronologer/the
fragment-intensity cache reused across every combo, no recompute — see
`fragment_intensity.md`). Real 1%-FDR "ions" (distinct `(peptide, charge)`,
`searchops`'s own definition):

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
mokapot's learned model, never SAGE's own score. Best combination found at
this point: RT+IIM + fragment-intensity + mokapot(xgboost), 32,088 ions —
beats every previously-recorded number for this dataset (old
`ranking_score`-focused best was 25,409 ions, RT-only, no mokapot, no real
fragment-intensity signal — see `software/sage/devel_fixed`'s `CLAUDE.md`).

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

**Considered, built, and deleted**: an `IsolationForest`-based mokapot
plugin (was `git/mokapot-isolationforest-plugin`) — mokapot's `Model.fit`
retrains its estimator on the current confident-target/decoy subset each
iteration via `estimator.fit(samples, iter_targ)`, but `IsolationForest.fit(
X, y=None, ...)` is a plain unsupervised sklearn estimator that silently
discards `y`, so it never actually learns a target-vs-decoy boundary. Built
and tested in isolation (unit tests passed, entry point registered
correctly) but never wired into a real job or committed to a git repo of
its own — the user judged the underlying idea unsound before a real run
was attempted, and the directory was deleted outright (2026-09-01). Its
uninstall was incomplete, though: `venvs/mokapot`'s editable-install
registration (`__editable__.mokapot_isolationforest_plugin-0.1.0.pth` +
matching `.dist-info`) survived the directory deletion, pointing at a path
that no longer existed. mokapot enumerates every registered
`mokapot.plugins` entry point at startup regardless of which `--plugin` is
actually requested, so any real (non-cached) mokapot run —
`--plugin xgboost` included — crashed with
`ModuleNotFoundError: No module named 'mokapot_isolationforest'` until
those two stale files were removed (2026-09-01, found while verifying the
RT-heteroscedastic-spline change below, since it was the first real mokapot
recompute since the plugin's deletion).

## RT-heteroscedastic tolerance spline: real F9477 measurement (2026-09-01)

`plans/rt_heteroscedastic_tolerance_spline.md`: `correct_precursors_rt`'s
`rt_tol_sec`/`rt_sigma_sec` fit went from one flat (RT-independent) window
to a 10-knot RT-dependent spline (SAGE's `LinearSpline`/`ValueTolSpline`
already supported arbitrary node counts; `rt_sigma_sec` — and
`Feature::delta_rt_z2_external`'s scale — needed a small Rust change to
follow suit, see that plan for the exact diff). Motivated by a real
residual-vs-RT diagnostic showing genuine heteroscedasticity even after
`spectrum_q`-filtering (robust sigma ~3.7s at low RT vs ~6.4s at high RT,
quintile-binned) — quantifying, not contradicting, the small "genuine
widening at the tail" `plans/better_sage_filtering.md`'s B.5 (2026-08-20)
had already flagged but judged not worth modeling at the time.

Same `jobs/f9477_best.toml` config, same `SpecId`→`used.pin` charge-join
counting method, only the code changed:

| | PSMs | peptides | ions |
|---|---:|---:|---:|
| flat window (previous) | 96,909 | 31,125 | 33,640 |
| 10-knot spline | 99,027 | 31,755 | **34,431** |
| Δ | +2,118 | +630 | **+791 (+2.4%)** |

For scale: a comparably-sized flat-window change tested the same session
(`tolerance_percentiles` `[1,99]` → `[0.5,99.5]`, no spline) moved ions by
only +21/33,640 — noise-level. +791 is well above that, and above the
~1.5% (~500-ion) run-to-run mokapot stochasticity this session separately
observed on an *unchanged* config (mokapot's train/test split isn't
seeded). Single-run comparison, not yet repeated to harden the effect
size — see the plan file's verification section for what a fuller
confirmation would need.

## Considered, built, not adopted: fragment m/z RT-bias term (2026-09-02)

Analogous idea to the RT-tolerance spline above, but for `recalibrate_pmsms_mz`
(fragment m/z correction) instead of `correct_precursors_rt`: fragment ppm
residual after the existing mz-only P-spline fit still carries a small
RT-dependent bias (spectrum-level median drift from about -0.4ppm at low RT
to -1.1ppm at high RT, checked via the same confident-hit fragment set used
for the diagnostic above — real, survives aggregating to one row per PSM to
rule out a fragments-per-spectrum-count artifact, but an order of magnitude
smaller effect than the RT-tolerance-spline's heteroscedasticity finding).

Implemented as a second, sequential fit term: `f_mz(fragment_mz) +
f_rt(precursor_rt)`, since every fragment in a spectrum shares its
precursor's RT — a per-spectrum additive shift, not a per-fragment one.
`recalibrate_pmsms_mz` gained a new `precursors` input (`PreSageFilteredPrecursors`,
wired as `P.search_precursors`, sibling of the existing `mz_pmsms` input —
no new pipeline node, just a new edge into the existing rule) to supply both
the fit-time RT signal and, via `cut_and_index_precursors`'s existing
`fragment_spectrum_start`/`fragment_event_cnt` CSR columns, an apply-time
per-fragment RT broadcast (`timstofu.timstofmisc.broadcast_precursor_values_to_fragments`,
new). Applied via a new `apply_mz_recalibration_mz_rt` (sums both correctors)
alongside the existing single-dim `apply_mz_recalibration`.

**Found and fixed one real bug during verification**: `precursors_ds["rt"]`
is raw Bruker frame time in *seconds* (`timstofu.candidate_postprocessing
.annotate`'s `frame2rt` lookup, range ~0.6-508s on F9477), while
`sage_results_tsv`'s own `rt` column — what the RT term is fit on — is
*minutes* (~0.02-8.2 on F9477). A first end-to-end run silently applied a
near-constant ~-1ppm shift to almost every fragment (the fit's flat-clamped
edge value, since the mismatched-unit apply-time grid mostly fell outside
the fit's real domain) and cost 678 ions. Fixed by converting
`precursors_ds["rt"]` to minutes (`/ 60.0`) before use.

**Real F9477 measurement, after the fix**, same `jobs/f9477_best.toml`,
same counting method, against the RT-heteroscedastic-spline baseline above:

| | PSMs | peptides | ions |
|---|---:|---:|---:|
| baseline (no fragment RT term) | 99,027 | 31,755 | 34,431 |
| + fragment mz+rt correction | 97,826 | 31,326 | 33,905 |
| Δ | -1,201 | -429 | **-526 (-1.5%)** |

Right at the ~500-ion / ~1.5% mokapot run-to-run noise floor — not a clean
regression, but clearly not an improvement either, consistent with the
bias magnitude being an order of magnitude smaller than what made the
RT-tolerance-spline change worthwhile. **Decision: keep the code (it's
correct, tested, and the `precursors` wiring is harmless when unused) but
`jobs/f9477_best.toml` is not updated to use it** — same
ask-before-updating rule as `bestrun.md`, and this result doesn't clear
the bar.

### Follow-up: swapped the sequential fit for real backfitting (2026-09-07)

The "sequential" fit above (`f_mz.fit(mz, ppm)`, then `f_rt.fit(rt,
residual)` once) wasn't actually the GAM/backfitting this module already
had a real implementation of (`fit_additive_correction`'s
`pspline_additive` branch, present since before this feature but never
wired to any pipeline rule or test — genuinely dead code until now).
Sequential fitting only approximates the joint fit when the dims are
uncorrelated; real backfitting iterates both dims against each other's
current residual (`config.get("backfit_iters", 15)` rounds, re-centering
each component to zero mean every round so the split stays identifiable).

Two prerequisite fixes, both committed separately before the swap:
- Removed `fit_additive_correction`'s unused, untested `xgboost_additive`
  branch (dead code, never called, no config ever referenced it).
- `pspline_additive` hardcoded one shared `bin_width_da`/`lam1`/`lam2`/
  `degree` across every dim — wrong once dims stop sharing a scale
  (fragment m/z is Da-scale, ~200-1700; RT is minute-scale, ~0-8). Each
  now accepts either a scalar (unchanged behavior) or a `{dim: value}`
  dict. Verified on synthetic mz+rt data: a shared `bin_width_da=10`
  (sane for mz) gave ~16x worse RT-component fit than per-dim
  `{"mz": 10.0, "rt": 0.5}`.

`recalibrate_pmsms_mz` now calls `fit_additive_correction(dims=["mz",
"rt"], config={"model": "pspline_additive", ...})` directly instead of
two separate `PSplineModel.fit()` calls; `config["fragment_model"]["class"]`
must be `PSplineModel` now (backfitting needs the same smoother family for
both dims each round, unlike the old mz-only fit's pluggable `build_model`).

**Real F9477 measurement**, same job/counting method, against both
earlier numbers:

| | PSMs | peptides | ions |
|---|---:|---:|---:|
| baseline (no fragment RT term) | 99,027 | 31,755 | 34,431 |
| + fragment RT term, sequential (2026-09-02) | 97,826 | 31,326 | 33,905 |
| + fragment RT term, backfitting (2026-09-07) | 99,112 | 31,798 | **34,498** |
| Δ backfitting vs sequential | +1,286 | +472 | **+593 (+1.7%)** |
| Δ backfitting vs baseline | +85 | +43 | **+67 (+0.2%)** |

Backfitting recovers the sequential fit's regression entirely (the
one-pass approximation really was measurably worse here, not just
theoretically), but the net effect of the fragment RT term itself is
still a wash: +67/34,431 sits well inside the ~500-ion/1.5% noise floor.
**Same decision as before: keep the code, don't adopt into
`jobs/f9477_best.toml`** — now on firmer footing, since backfitting rules
out "the fit method itself is broken" as an explanation for the earlier
non-result.
