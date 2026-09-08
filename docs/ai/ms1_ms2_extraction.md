# C++ MS1/MS2 extraction targets

## C++ MS1 extraction

`ionmaiden_pipeline` uses `git/ionmaidenmetal/build/tdf2ms ms1` for `ms1_events`,
passing the allocated thread count, `--paced-writeback-mib 1024`,
`--implicit-tof-urt`, and `--overwrite`. Its payload contains only uint32 scan and
intensity columns. Both TOF/URT split indexes remain unchanged. Completion guards
check both indexes, scan/intensity types, and absence of physical TOF/URT columns.

Timstofu derives raw-event centers from index spans; existing window processors
still receive complete local scratch tensors. Scale estimation retains its full
local tensor outputs for Python analysis. Selected precursor tables retain TOF/URT
coordinates. Legacy materialized datasets remain readable, and the standalone
converter still defaults to materialized output. No pipeline storage-mode knob.

The 1024-MiB budget makes MS1 conversion output clean before exit, avoiding the
previously measured dirty-page tail. See
`git/timstofu/docs/ai/compact_ms1.md` (relative to the pipeline root) for consumer
coverage and the 2026-09-08 full-chain equivalence/runtime measurements.

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
