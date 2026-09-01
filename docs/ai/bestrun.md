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
| 98,514 | 28,022 | **34,133** |

See `mokapot_integration.md`'s full ablation table for how this was
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
