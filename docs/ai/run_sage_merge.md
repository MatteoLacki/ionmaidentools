# `run_sage`/`run_sage_with_predicted` merged into one rule (2026-08-26)

`run_sage_with_predicted` (a separate rule, `--predicted-rt`/`--predicted-iim`
always real DAG edges resolving to the `NoPrediction` sentinel when a
dimension was inactive, flag inclusion driven by a scalar `dimensions`
tuple) no longer exists — `run_sage` itself now takes
`predicted_rt: PredictedRt | None = None`/`predicted_iim: PredictedIim | None
= None` directly, using necroflow's "mixed Node/value inputs" support
(`docs/rules.md`) added after the original RT/IIM independent-dimensions
split (see `recalibration_modes.md`). One rule now covers pass-1/mode-1/mode-2
(neither prediction, `predicted_rt`/`predicted_iim` simply omitted) and mode
3's final pass (either/both, passed as real Nodes) — the pipeline factory's
mode-3 `else` branches assign `P.predicted_rt = None`/`P.predicted_iim = None`
directly (a plain Python attribute, not a DAG edge — confirmed safe by reading
`Pipeline.__setattr__` itself: non-`Node` values always fall through to a
plain instance attribute, no validation, no DAG registration) instead of
`write_no_prediction_marker(...)`.

`_run_sage_command` (a Python command callback, not a static template,
since necroflow's `@command` string templates can't conditionally include
a flag) decides which `--predicted-*` flags to add purely from
`args.inputs.predicted_rt/predicted_iim is not None` — `args.inputs.*`
preserves a plain `None` verbatim for an unresolved mixed-input default
(only managed Nodes resolve to a real `Path`), so this needs neither the
old `dimensions` scalar nor the `NoPrediction` sentinel just to keep the
parameter a real DAG edge. `NoPrediction` itself is not removed — it's
still required for `update_sage_config_rt_iim`'s `rt_tolerance`/
`mobility_tolerance` inputs, which stayed plain required `NodeType`s (out
of scope for this merge; properly avoiding the sentinel there would mean
splitting that rule into up to 3 static-template variants instead of one
Python-callback command).

Verified against real F9477 data end-to-end (`ab_evict`-style RT-recalibration
job, both `run_sage` invocations it produces): the filtered/mode-2 pass's
realized Sage config has `"predicted_rt": null, "predicted_iim": null,
"rt_tol": null` (flag genuinely absent, 26,673 target PSMs at 1% FDR); the
final mode-3 pass's has `"predicted_rt": ".../predicted_rt.parquet"` with a
real path, Sage's own log confirming `"loaded 6300882 predicted sequence ->
rt entries"`, and populated `rt_sigma`/`rt_tol` (90,125 target PSMs at 1%
FDR). Both branches of the merged rule's Python callback exercised for
real, not just dry-run/import-checked.

`run_sage` was later extended again (`fragment_intensity.md`'s "Wired into
`run_sage`" section) to also take the fragment-intensity cache as two more
optional Node inputs, using the same mixed Node/`None` mechanism this merge
established.
