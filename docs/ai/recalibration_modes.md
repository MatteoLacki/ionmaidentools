# Recalibration: m/z, RT and IIM, each gated by its own table

Without `[recalibration]` the workflow runs one SAGE pass on `search_pmsms`/
`search_precursors`. With it, a first (recalibration) pass searches a precursor sample;
then each of `[recalibration.mz]`, `[recalibration.rt]`, `[recalibration.iim]` that is
present fits its correction on that pass's confident PSMs, and one final `run_sage` searches
with all of them. Table presence is the switch; there is no `dimensions` list
(`plans/better_sage_filtering.md` B.6, `plans/rt_iim_independent_dimensions.md`).

- **m/z**: `fit_mz` (below).
- **RT**: `check_speclib_rt`, then `fit_precursor_rt` (featureprediction
  `correct_precursors_rt`): `predicted_rt`, the corrected `rt` column, `rt_tolerance`.
- **IIM**: `check_speclib_iim`, then `fit_precursor_iim`: the corrected `inv_ion_mobility`
  column and `mobility_tolerance`; SAGE reads the SPECLIB's IIM itself.
- Each correction is one column evaluated on the raw `search_precursors`;
  `update_precursors` sets them all into `recalibrated_precursors`
  (`plans/precursor_corrections.md`). `update_sage_config` writes each tolerance, and for
  RT/IIM also `rt_sigma_sec`/`iim_sigma`: SAGE's `Input::build()` requires
  `predicted_rt` + `rt_tol_sec` + `rt_sigma_sec` all or none (same for IIM).

**Tolerances are per table.** Each of the three tables carries its own required
`tolerance_percentiles` and optional `tolerance_method` (`"theoretic"`, default,
`median ± z·robust_sigma`; or `"empiric"`), never inherited: F9477 mass-error residuals are
right-skewed, RT residuals much less. `theoretic` needs `hi_pct` strictly in `(50, 100)`;
`[0, 100]` fails loudly. Implementations: featureprediction `tolerance.select_tolerance`,
searchops `recalibration._select_tolerance`.

**Exports match the final search.** `current_precursors`/`current_tof2mz_table` are
reassigned as stages run (hence `current_`, not `final_`) and feed
`convert_search_pmsms_to_mzml`/`convert_search_pmsms_to_mgf`, so MGF/mzML headers and m/z are
what the final SAGE pass searched. MGF export is wired only when `[mgf].config_path` is
present; the path is a scalar input, so editing that file in place does not invalidate a node.

## Recalibration spectra from raw MS2 top cells: `[recalibration_top_cell]`

Optional, gated by presence. The recalibration pass then searches each sampled precursor's
most probable (frame, scan) raw MS2 spectrum in place instead of the main pmsms:
`top_cell_precursors` points `fragment_spectrum_start`/`fragment_event_cnt` at that slice of
`events.ms2/data.mmappet`, and `run_sage` takes the `Ms2Events` store (recognised by its
`events.ms2` filename) with `--tof2mz`. These are the spectra mkpmsms'
`top_probable_frame_scan` would copy into a pmsms. The fits still apply to the main pmsms
and `search_precursors`.

## One m/z fit for precursors and fragments: `fit_mz`

Every `run_sage` reads fragment m/z as `tof2mz[tof]` from the pmsms' `tof` column
(`plans/fragment_mz_corrected_tof2mz_table.md`).

- `fit_mz` (searchops `recalibrate_mz`) fits `[recalibration.mz].precursor_model`, ppm
  error against m/z, on confident precursors and applies it to both the precursor `mz`
  column and `recalibrated_tof2mz_table` (float64, `table[t] / (1 + f(table[t])·1e-6)`).
  Each tolerance comes from its own residuals under that fit; the fragment one is the only
  reader of the first pass's matched fragments. A separate fragment fit gave the same
  F9477 mokapot counts; a fixed ±15 ppm fragment tolerance cost 1.5% of peptides.
- The table is float64: a float32 one rounds every m/z twice and moves 25% of F9477 peaks
  by one float32 ulp.
- Exports read m/z the same way (`--tof2mz`, `plans/exports_from_tof2mz_table.md`); nothing
  in the pipeline materializes `mz` (timstofu's `materialize_pmsms_mz` is a standalone
  tool). The pipeline writes no `fragment_shift_ppm`; SAGE and the writers apply one if
  the precursors carry it.
- `convert_search_pmsms_to_mgf` runs `chattr +m {workdir}` first (no-op off btrfs): zstd
  compression of ~24 GiB of text during writeback throttles the writer (F9477 90.9 ->
  74.6 s).
