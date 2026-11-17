# Mokapot integration

## Forked mokapot (`git/mokapot`), real CLI

`mokapot()` runs the real `venvs/mokapot/bin/mokapot` CLI from our fork
(`github.com/MatteoLacki/mokapot`, branch `mokapot-necroflow-integration`, based on
`v0.10.0`; upstream `main` has a different `OnDiskPsmDataset` architecture incompatible
with the patches and with `mokapot-xgboost-plugin`'s `Model`/`brew()` API). The `Makefile`
installs it editable. Two patches:

- `read_percolator` reads the PIN with DuckDB (~2x upstream's per-line `str.split`), same
  rows, dtypes and values.
- `--seed` seeds `brew()`'s rng (upstream's `np.random.seed` never reaches
  `np.random.default_rng()`). Same seed, same PIN: byte-identical output.

**Seeding**: `[mokapot].seed` (default 1, mokapot's own) also feeds `--xgboost_seed`/
`--lightgbm_seed`, the plugins' estimator `random_state`. A seeded run is
byte-deterministic on an identical PIN; PSM row order or tiny feature changes still swing
~200 peptides/ions.

**Plugins**: `[mokapot].plugin` (`"xgboost"`, `"lightgbm"`, or omitted for mokapot's linear
SVM; `""`, never `None`, for none: see `AI.md`). `git/mokapot-lightgbm-plugin` is a
sibling of `mokapot-xgboost-plugin`. On F9477 CPU, lightgbm ties xgboost in time and ions;
`HistGradientBoosting` is 1.8x slower; bagged XGBoost (`xgboost_bagging_*`, soft-voting
`BaggingClassifier` of 5 subsamples, 30k rows each) is ~3x slower with no gain (3k rows
is worse on both).

**GPU**: xgboost's plugin passes `device="cuda"` and falls back to CPU with a warning.
LightGBM's pip wheel has no GPU code; the `Makefile` builds it with `USE_CUDA=ON` when
`nvidia-smi`/`/dev/nvidia0` exists, and the plugin's `_select_device()` decides cuda vs
cpu once by a smoke-test fit (mokapot calls `.fit()` 30 times with no retry hook). Only the
CPU path is verified; the AMD iGPU here is usable by neither.

**Threads**: `_mokapot_command` sets `--max_workers = min(folds, threads)` and
`--{xgboost,lightgbm}_n_jobs = threads // max_workers` from the threads the call reserves,
not job config: execution tuning only, and necroflow's provenance hash covers declared
config, not the rendered command. F9477: 161.8 -> 54.8 s, ions unchanged. What remains is
mokapot's own Python bookkeeping (`_update_labels`, `_find_best_feature`) about as much
as the XGBoost fits. `folds` is `[mokapot].folds` (default 3).

**`--keep_decoys`** writes `mokapot.decoy.{psms,peptides}.txt` (~26 MB on F9477), undeclared,
for ion-level FDR (`fdr_levels.md`); the target files are byte-identical without it.

## Leakage-safe PIN (`scripts/mokapot_pin_adapter.py`)

SAGE's `results.sage.pin` holds `posterior_error`, SAGE's in-run label-conditioned LDA
score: feeding it to mokapot is leakage. `--mode sage` projects the PIN onto a fixed safe
column registry with DuckDB; `rt_source`/`iim_source` (`"external"` when that dimension's
prediction exists for the job) select the external-prediction columns. `--mode passthrough`
(drop `FileName` only) has no caller. Design: necromerge2 `plans/mokapot_leakage_safe_pin.md`.

## A stale plugin registration breaks every mokapot run

mokapot enumerates every registered `mokapot.plugins` entry point at startup, whichever
`--plugin` is asked for. Deleting a plugin's directory without uninstalling it leaves its
editable `.pth` + `.dist-info` in `venvs/mokapot`, and every run then fails with
`ModuleNotFoundError`; remove both files.
