# Mokapot integration

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
