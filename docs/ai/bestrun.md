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

## Current best (as of 2026-09-01)

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
