# C++ MS1/MS2 extraction targets

## C++ MS1 extraction

`ionmaiden_pipeline` uses `git/ionmaidenmetal/build/tdf2ms ms1` for `ms1_events`,
passing the Necroflow-allocated thread count, `--paced-writeback-mib 1024`, and
`--overwrite`. The 1024-MiB budget makes the production MS1 output clean before
exit, avoiding the measured 36-37 second dirty-page tail without a new pipeline
configuration knob. The rule validates both
mmappet split-index datasets and the event-dataset schema before completion. The in-memory
slice contract is unchanged from `d2ms1`; both indices now use the common mmappet-array
format, so all peak-picking and precursor consumers retain the same `(tof,urt,scan)` slices.

The converter decompresses each MS1 frame once, stable-count-sorts it to a `(tof,scan)` run, then
uses event-balanced TOF shards to merge runs into disjoint sequential mmap regions while building
the split index. Exact comparisons with `d2ms1` passed on F9477, F9468, and B6699; observed wall
speedups were 3.78x, 2.76x, and 3.43x. Python `d2ms1` remains installed as a reference/debug tool,
not the pipeline default.

## Compact TOF/scan/frame MS2 extraction benchmark target

`ionmaiden_pipeline` exposes `ms2_tsf_events` as an independently requestable target.
Its `tdf2ms2_tsf` rule runs `git/ionmaidenmetal/build/tdf2ms ms2-tsf` with the
Necroflow-allocated thread count and produces `events_ms2_tsf.mmappet`. The command
validates the payload and flat TOF/scan row-start mmappet schemas plus `stats.json`
before Necroflow marks the node complete.

This target deliberately has no downstream consumer yet. Existing pseudo-MS/MS stages
continue to use `ms2_events` and the established frame/scan layout; request
`ms2_tsf_events` explicitly to benchmark or inspect the compact `(tof,scan,frame)`
layout without changing production results. Necroflow records wall time and output size
automatically in the node `.rip/run.toml`.

## Query-oriented TOF/frame/scan MS2 extraction benchmark target

`ionmaiden_pipeline` also exposes `ms2_tfs_events` independently. Its
`tdf2ms2_tfs` rule runs `git/ionmaidenmetal/build/tdf2ms ms2-tfs` with the
Necroflow thread allocation and publishes `events_ms2_tfs.mmappet`. Completion
checks cover intensity data, packed-scan schema/shape, TOF index schema/shape,
frame index schema, and `stats.json`.

This node has no downstream consumer and does not replace `ms2_events` or
`ms2_tsf_events`. Request it to benchmark frame-selective box queries against the
`(tof,frame,scan)` layout. Scans use four little-endian uint10 values per five bytes;
TOF and frame headers locate the narrow ranges without expanding scan IDs. F9477
completed in 5.07 s and 1.554 GB, with all 160,689,740 events matching TSF exactly.
