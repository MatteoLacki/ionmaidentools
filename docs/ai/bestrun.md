# `jobs/f9477_best.toml`: the current best-known ion-wise config

`jobs/f9477_best.toml` (repo root, not under `git/ionmaidentools/`) is a
standing, deliberately-maintained job config — not a one-off experiment like
`jobs/f9477_ablation/*.toml`. It tracks whichever real, verified F9477
configuration has produced the highest **unique-ion** (`searchops`'s
definition: distinct `(peptide, charge)` at 1% FDR) count found so far,
via mokapot on top of SAGE. It exists so there's always one obvious,
current-best config to point people (or agents) at, instead of having to
dig through `docs/ai/mokapot_integration.md`'s ablation table history to
figure out which combination currently wins.

## Config change 2026-10-20: `merge_fragment_windows = true` (numbers pending)

`[sage] merge_fragment_windows = true` (git/sage `docs/ai/merged_fragment_windows.md`): the
prescore counts each fragment once however many peaks lie within tolerance of it. Enabled at the
user's request for correct counting, not for a higher ion count. The job's final SAGE and mokapot,
rerun by hand on the 2026-10-15 nodes' inputs with one binary (necromerge2
`temp/peptide_views/merge/best/`): SAGE 1% 94,053 / 20,705 / 22,717 -> 95,216 / 21,046 / 23,031
(spectra / peptidoforms / ions), mokapot 121,372 / 29,120 / 31,778 -> 122,144 / 28,884 / 31,680,
within mokapot's noise.

## Current best (2026-10-29)

`[score_competition]` added (each precursor's ToF scores against its neighbours'; the settings of
the former `f9477_best_plus_tof_filter`, `score_margin = 0.05`), on the SPECLIB store, after the
precursor RT unit fix (necromerge2 `plans/precursor_corrections.md`) and with the RT correction's
robust fit (featureprediction f4da93e, the default). Each level at 1% FDR from its own
competition (`*.fdr_levels.tsv`; run database rows 95/96):

| | spectra | peptidoforms | ions |
|---|---|---|---|
| SAGE | 78,912 | 19,711 | 21,312 |
| mokapot | 100,810 | **29,939** | **32,521** |

With DIA-NN 2.7.0 IIM instead of none: mokapot 29,101 / 31,519 (rows 90/91), plus
`[score_competition]` 28,923 / 31,174 (rows 92/93). A full node-store rebuild moves these by up
to ~2% (mkpmsms row order), so compare variants on shared upstream nodes. Rows 19–69 ran on MS1
nodes from before timstofu 42f4b0a and aren't comparable. These levels count differently from
the entries below; don't compare across.

## Previous best (2026-10-07)

Same config as below plus mkpmsms `"ms2-peak-picker" = "local_max_intensity"` (was
`local_max_score`): each fragment peak is written at its most intense TOF bin instead of its
best-scoring one. Chosen by the peak-position experiment on F9477 and F9468 (necromerge2
`plans/tof_centroid_experiment.md`; findings in git/ionmaidenmetal `docs/ai/tof_centroid.md` on
branch `experiment/tof-centroid`). `jobs/f9468_benchmark.toml`, the F9468 port, changed the same way.

Unique target peptides (`mokapot.peptides.txt`, q ≤ 1%) and ions (distinct (peptide, charge) of
target PSMs, q ≤ 1%):

| run | PSMs | peptides | ions |
|---|---|---|---|
| main, production binaries, mokapot seed 1 (2026-10-07) | 100,166 | 28,583 | **34,995** |
| experiment run, same config, mokapot seeds 1–3 mean | | 28,435 | 34,802 |
| previous config (`local_max_score`), seeds 1–3 mean | | 27,915 | 34,286 |

SAGE alone (rank 1, peptide q ≤ 1%): 21,475 peptides / 25,142 ions (was 21,122 / 24,652).
On F9468: mokapot seeds 1–3 mean 30,689 / 38,454 (was 29,740 / 36,625). One mkpmsms run per
config; mkpmsms' row order alone moves mokapot by ~1%. The peptide column is not comparable
with the 2026-09-01 entry below, which counted differently.

## Previous best (2026-09-01)

RT-only recalibration + fragment-intensity + mokapot(xgboost), no IIM:

| PSMs | peptides | ions |
|---|---|---|
| 99,027 | 31,755 | **34,431** |

`jobs/f9477_best.toml` itself did not change for this update — the RT
tolerance/sigma fit it drives went from one flat window to a 10-knot
RT-dependent spline (`plans/rt_heteroscedastic_tolerance_spline.md`), so the
same job config now produces a different, better number purely from the
code change. Previous recorded number: 98,514 / 28,022 / 34,133 (flat
window, `SpecId`→`used.pin` charge-join counting method, unchanged between
the two measurements). Single-run comparison — mokapot's own train/test
split isn't seeded, and this session separately measured ~1.5% (~500 ion)
run-to-run swing on an unchanged config, so treat this as a real but not
yet statistically hardened improvement.

See `mokapot_integration.md`'s full ablation table for how the underlying
config (RT-only + fragment-intensity + mokapot xgboost, no IIM) was
determined (IIM was found to actively hurt, not just fail to help) and the
`[tof_score_filter]` ablation that came out neutral on top of this exact
config.

## Maintenance rule

**When a new run beats `jobs/f9477_best.toml`'s current ion count, ask the
user before updating this file to match it.** Do not silently overwrite
`jobs/f9477_best.toml` (or this doc's recorded numbers) just because a new
experiment scored higher — confirm with the user first, since "best" here
is a deliberate, standing reference point other work may already be
built against, not just a running maximum to chase automatically.

When the user does confirm an update:

1. Overwrite `jobs/f9477_best.toml`'s content to match the new
   best-performing job exactly (copy it, don't hand-edit toward it).
2. Update the numbers in this file's "Current best" section and the date.
3. Leave the historical ablation entries in `mokapot_integration.md` alone
   — that table is a record of what was tried, not something to prune when
   the best changes.
