# AI_CHANGELOG — ionmaidentools

## 2026-08-20 — recalibration modes (B.6)
- What: three modes; `current_*` exports made MGF/mzML headers match what SAGE searched
  (before, exports always used raw `search_precursors`, even after m/z recalibration).

## 2026-08-25 — table presence instead of `[recalibration.rt_iim].dimensions`; per-table tolerances
- What: `[recalibration.rt_iim]` (with `dimensions`) removed; `tolerance_percentiles`/
  `tolerance_method` moved into `[recalibration.mz]`/`.rt`/`.iim` (were shared, mz's at the
  root). Jobs rewritten without a shim: `f9468_fragpipe`, `f9477_gam_test`, `b6699_gam_test`,
  `f9477_fragpipe`, `short_test_recal`. Jobs on the older flat `model = "..."` schema
  (`b6699_test_recal`, `quick_test_10k_recal`, `full_recal_grid`, `inspect_recalibrate`,
  `quick_test_10k_recal_grid`, `f9468_test_recal`, `short_test_recal_*grid`) were already
  broken and left alone.
- Why: `dimensions` restated the tables; the name read as "both RT and IIM".

## 2026-08-26 — `run_sage_with_predicted` merged into `run_sage`
- What: predictions became mixed Node/`None` inputs; the scalar `dimensions` tuple and
  the `NoPrediction` sentinel for `run_sage` (`write_no_prediction_marker`) went.
- Measured: F9477 RT job, both passes: pass 1 config `"predicted_rt": null` (26,673 target
  PSMs at 1%), final pass loaded 6,300,882 RT entries (90,125).

## 2026-08-26 — `update_sage_config_rt_iim` writes `rt_sigma_sec`/`iim_sigma`
- Why: with only `rt_tol_sec`/`mobility_tol`, SAGE rejected every mode-3 config; found on a
  real F9477 run, since smoke `--requests` never reached the final `run_sage`.
- Also: `[recalibration.iim]` charges defaulted to `cfg.sage.get("precursor_charge", (2, 4))`.

## 2026-08-31 — fragment-intensity cache in the pipeline
- What: `predict_fragment_intensity` (`FragmentIntensityCache`, this repo's first `mutable=True`
  rule) and `export_fragment_intensity_for_sage` (ordinary rule on a mutable parent), wired
  into the final-pass `run_sage` behind a presence-only `[fragment_intensity]` table, never
  the calibration pass.
- Why: SAGE's `ms2_*` features had computed on empty input (all `0.0`) on every job.
- Gotcha then: Node-shaped inputs had to be passed positionally (`TypeError: run_sage:
  unexpected inputs`), hence `None, None` placeholders; fixed in necroflow 2026-09-01
  ("Accept Node inputs by keyword"), call sites moved to keywords with zero hash changes.
- Measured: full F9477 fill ~43 min against Koina (~2.8 GB); export 116.6 s / 1.9 GiB, then
  18.0 s / 479 MiB as pointers; 2,028,474 of 18,902,646 (sequence, charge) pairs absent.

## 2026-08-31 — `[mokapot].plugin`, leakage-safe PIN, F9477 ablation grid
- What: `plugin` from config (was hardcoded `"xgboost"`); `""` for none after an explicit
  `None` crashed a 2x2x3 grid in tomlkit after the mokapot fit finished.
  `mokapot_pin_adapter.py --mode sage` drops `posterior_error` (`plans/mokapot_leakage_safe_pin.md`).
- Measured (ions at 1%, `jobs/f9477_ablation/`): SAGE alone 21,873 (no RT/IIM) / 22,500
  (RT+IIM), unchanged by fragment intensity; mokapot default 31,710, xgboost 32,088. RT-only
  later reached 34,133: IIM hurt (SAGE alone 24,642 -> 22,500 with IIM). `[tof_score_filter]`
  on top: 34,020 (neutral).

## 2026-09-01 — RT/IIM prediction caching
- What: `fill_rt_prediction_cache`/`fill_iim_prediction_cache` (`mutable=True`, single
  output) feeding `predict_rt`/`predict_iim` through featureprediction's `--cache-path`.
- Why: every job made a full uncached Koina call even when `dumped_peptides` overlapped fully.
- Measured: the first fill costs 255.3 s (vs 92.1 s uncached); a second job sharing
  `dumped_peptides` reuses both nodes, so RT prediction costs nothing from the second job on.

## 2026-09-01 — IsolationForest plugin deleted; 10-knot RT tolerance spline
- What: an `IsolationForest` mokapot plugin was built and deleted unused: mokapot's
  `estimator.fit(samples, iter_targ)` cannot teach an unsupervised estimator that drops `y`.
  Its leftover editable registration crashed every mokapot run until removed.
- What: `rt_tol_sec`/`rt_sigma_sec` became a 10-knot RT-dependent spline
  (`plans/rt_heteroscedastic_tolerance_spline.md`); robust sigma ~3.7 s at low RT vs ~6.4 s high.
- Measured: `f9477_best` 33,640 -> 34,431 ions (+2.4%); `[1,99]` -> `[0.5,99.5]` moved +21.

## 2026-09-01 — best: RT-only + fragment intensity + mokapot xgboost, no IIM
- Measured: 99,027 PSMs / 31,755 peptides / 34,431 ions (`SpecId`→`used.pin` charge-join
  counting) with the 10-knot RT spline; flat window 98,514 / 28,022 / 34,133.

## 2026-09-01 — mode-2/mode-3 branches flattened (708cc5b); `final_*` later renamed `current_*`
- What: one final `run_sage` for zero, one or both of RT/IIM; `update_sage_config_rt_iim`
  called unconditionally (`dimensions=()` was a plain copy).
- Measured: a mode-2 job's `sage_config.json` byte-identical, counts identical (75,888 PSMs /
  21,873 ions); every mode-2 final `run_sage` reran once (new parent).

## 2026-09-17 — mokapot fork, seeding, plugins, threading
- What, in order the same day: auto-derived `--max_workers`/`--xgboost_n_jobs`
  (161.8 -> 54.8 s, 34,693 -> 34,764 ions); `scripts/run_mokapot.py` bypassed the CLI's PIN
  parser (`pd.read_csv` 1.53 s vs 8.94 s; 46.0 s) and saved fold models
  (`MokapotModelFold1/2/3`), with a `[mokapot].model` choice; a first draft wrote decoy files
  by mistake (`decoys=True`). Then the parser and `--seed` were fixed in the fork, the real
  CLI came back, `run_mokapot.py` and the fold-model nodes were deleted, and lightgbm became
  its own plugin.
- Measured: seeded xgboost 46.0 s, 99,147 / 31,802 / 34,486; lightgbm 48.2 s, 100,931 /
  31,999 / 34,682. Unseeded xgboost runs had spread 33,636-34,764 ions (~3%).
- Profile (folds serial, 67.8 s): `XGBClassifier.fit` 27.2 s, `_update_labels` 16.2 s,
  `_find_best_feature` 14.9 s (overlaps), `read_percolator` 8.9 s, `_get_scores` 4.9 s.
- Bagging: 5x30k ~130-140 s, 32,828-34,212 ions; 5x3k 140.9 s, 32,828.

## 2026-09-18 — sagepy_rescore branch removed
- What: the second `mokapot(...)` call site and its NodeTypes/rules went with
  `git/sagepy`/`sagepy-rescore`/`sagepy_ionmaiden_adapter` (`plans/sagepy-rescore.md`:
  "not worth it — GPU needed, poor results"); `--mode passthrough` lost its caller.

## 2026-09-29 — `[recalibration_top_cell]`
- What: replaced the mkpmsms route (`[recalibration_pseudomsms]`, 81ec5cb: a second mkpmsms
  run, cut-and-index, materialized m/z). SAGE on raw top cells reproduced it on F9477 (all
  47,917 PSMs and 235,629 matched fragments). Comparison: pipeline_analysis
  `docs/ai/recalibration_paths.md`.

## 2026-10-05 — exports through the tof2mz table; `mz` materialization removed
- What: `materialize_pmsms_mz`/`materialize_recalibrated_pmsms_mz`, `MzPmsms`/
  `RecalibratedPmsms` removed; SAGE's `mz`-column input removed.
- Measured (F9477): MGF and `--indexed` mzML byte-identical to materialize-then-write;
  directory-mode mzML identical per precursor (923,634 spectra); mzML 47.4 s vs 8.8 + 47.8 s,
  MGF 88.5 s vs 8.8 + 84.9 s, no 6.5 GB `mz` column. `chattr +m`: MGF 90.9 -> 74.6 s.

## 2026-10-06 — `sage_map_to_pmsms` removed
- What: SAGE exports `fragment_pmsms_row`; `score_comparison` reads it. `confident_psms →
  score_comparison` runs once after the final `run_sage` (was a per-branch closure).

## 2026-10-07 — best: mkpmsms `ms2-peak-picker = "local_max_intensity"`
- What: each fragment peak at its most intense TOF bin (was `local_max_score`); same in
  `jobs/f9468_benchmark.toml` (necromerge2 `plans/tof_centroid_experiment.md`).
- Measured: 100,166 PSMs / 28,583 peptides / 34,995 ions (seed 1); seeds 1–3 mean 28,435 /
  34,802 vs 27,915 / 34,286; SAGE alone 21,475 / 25,142 (was 21,122 / 24,652); F9468 30,689 /
  38,454 (was 29,740 / 36,625).

## 2026-10-13 — FDR levels and run recording inside `run_sage`/`mokapot`
- What: `sage-fdr-levels`/`mokapot-fdr-levels` as the commands' last step; mokapot gets
  `--keep_decoys`.
- Measured: every SAGE/mokapot node reran once; `jobs/f9477_marginal_cosine.toml` 9 nodes in
  4 min, counts matching its earlier run (SAGE 21,461, mokapot 28,644 peptidoforms; run DB
  rows 1-3).

## 2026-10-17 — compact peptide table; prediction nodes as mmappet
- What: `dump_compact_peptides` → `CompactPeptidesTable` for the sliced search; `PredictedRt`,
  `PredictedIim`, `FragmentIntensityForSage` became mmappet datasets (renamed outputs reran
  their rules once, warm caches, no predictions).

## 2026-10-20 — `f9477_best`: `merge_fragment_windows = true`
- Why: correct prescore counting (git/sage `docs/ai/merged_fragment_windows.md`), not more ions.
- Measured (2026-10-15 inputs, one binary): SAGE 94,053 / 20,705 / 22,717 -> 95,216 / 21,046 /
  23,031 (spectra / peptidoforms / ions); mokapot 121,372 / 29,120 / 31,778 -> 122,144 /
  28,884 / 31,680.

## 2026-10-21 — `--compact-peptides` follows SAGE's sliced default
- What: passed also when `database.sliced` is unset and none of `prefilter`, `wide_window`,
  `chimera` is set.

## 2026-10-26 — fragment intensities only from a SPECLIB
- What: `predict_fragment_intensity`, `export_fragment_intensity_for_sage`,
  `FragmentIntensityForSage` and `[fragment_intensity]` removed (a job with that table is an
  error); `[speclib]`, `check_speclib`, `run_sage(speclib=...)` (necromerge2 `plans/speclib_only.md`).

## 2026-10-26 — `[mokapot].folds`
- What: `folds` is a rule parameter passed as `--folds` (was a module constant never passed).

## 2026-10-27 — RT/IIM only from the SPECLIB
- What: `_CACHE_ROOT`, the Koina model/server constants, `rt_cache_path`/`iim_cache_path`,
  `[recalibration.*] server_url` and `PredictedIim` removed; `predict_rt`/`predict_iim`/
  `correct_precursors_*` take the `Speclib` node and find PSMs by `peptide_id`
  (necromerge2 `plans/speclib_rt_iim.md`).

## 2026-10-28 — corrections as one-column datasets on the raw precursors
- What: m/z, RT and IIM corrections each evaluated on raw `search_precursors` and set by
  `update_precursors` (were chained, RT applied at 1/60 of its size, IIM evaluated at RT in s;
  necromerge2 `plans/precursor_corrections.md`). Runs before this are not comparable.

## 2026-11-17 — `fit_mz`: one m/z fit for precursors and fragments
- What: `fit_mz` replaces `fit_precursor_mz` and `fit_fragment_mz`; `FragmentShiftPpm`,
  `[recalibration.mz].fragment_model` and the fragment RT-bias term removed; jobs edited.
- Measured: F9477 `f9477_best` mokapot 29,920 / 32,463; separate fragment fit 29,800; the
  fragment RT-bias term +0.2% ions (noise); RT/1/K0 in the precursor fit −1.5% mokapot.
