# Peptide table and predictions mapped by SAGE (2026-10-17)

Plan: necromerge2 `plans/mapped_peptide_table.md`. SAGE side: `git/sage`
`docs/ai/mapped_inputs.md`.

## `dump_compact_peptides` → `CompactPeptidesTable`

`dump_compact_peptides(fasta, dump_peptides_config, dump_peptides_binary)` runs
`dump_peptides --compact {compact_peptides}`: the same binary and the same
digest-only `dump_peptides_config` as `dump_peptides`, so one table per FASTA and
digest setting, shared by every job and both SAGE passes. Output
`compact_peptides.mmappet` (rows as mmappet columns plus `mass_chunks.bin`,
`wide_masses.bin`, `meta.json`).

`run_sage` takes it as an optional input (`compact_peptides`, `--compact-peptides`),
passed only when `cfg.sage.database.sliced` is true: the monolithic search builds
its own index and does not read the table. The workflow declares the node for
every job; necroflow runs it only for jobs whose `run_sage` asks for it.

It is a separate node from `dump_peptides` on purpose: `peptides.parquet` still
comes from sage's full digest (its `proteins` column lists every protein holding a
peptide, which the table does not store), and keeping the two apart leaves
`dump_peptides` and everything downstream of it untouched.

## Prediction nodes became mmappet datasets

`PredictedRt` (`predicted_rt.mmappet`), `PredictedIim` (`predicted_iim.mmappet`)
and `FragmentIntensityForSage` (`fragment_intensity_for_sage.mmappet`) are now
`MmappetDataset`s: the same values and row order as the parquets they replace, with
the `dumped_peptides_sha256` fingerprint in a `meta.json` sidecar. SAGE maps them
read-only instead of decoding them into memory. The renamed outputs make
`predict_rt`, `predict_iim`, `export_fragment_intensity_for_sage` and their
consumers run once more; the Koina caches are warm, so that costs no predictions.
