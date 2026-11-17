# `jobs/f9477_best.toml`: the current best-known ion-wise config

`jobs/f9477_best.toml` (repo root, not under `git/ionmaidentools/`) is a
standing, deliberately-maintained job config — not a one-off experiment like
`jobs/f9477_ablation/*.toml`. It tracks whichever real, verified F9477
configuration has produced the highest **unique-ion** (`searchops`'s
definition: distinct `(peptide, charge)` at 1% FDR) count found so far,
via mokapot on top of SAGE. It exists so there's always one obvious,
current-best config to point people (or agents) at, instead of having to
dig through the ablation history to figure out which combination currently wins. Earlier bests: `AI_CHANGELOG.md`.

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
3. Append the replaced best (date, config change, numbers) to `AI_CHANGELOG.md`.
