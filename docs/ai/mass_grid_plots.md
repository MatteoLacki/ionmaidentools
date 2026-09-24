# Mass-grid heatmaps and the peptide length histogram (2026-09-16)

Request-only diagnostics over a job's digest and its two fragment indices.
Nothing in the search consumes them; like `fragment_intensity_for_sage`, they
exist in the DAG as sinks and only run when a job lists them in `.requests`
(every job TOML in `jobs/` sets `.requests`, and necroflow defaults to *all*
sinks otherwise). `jobs/f9477_mass_grid_plots.toml` is `f9477_best.toml`
with only `.requests` changed.

## Labels

| label | rule | output |
|---|---|---|
| `peptide_length_histogram` | `plot_peptide_length_histogram` | percent of target peptides per residue length |
| `fragment_index_grid_<v>` | `bin_fragment_index` | SAGE fragment index: peptide monoisotopic mass x fragment neutral mass |
| `predicted_fragment_grid_<v>` | `bin_predicted_fragments` | Koina predictions: precursor m/z (charges pooled) x fragment m/z |
| `fragment_index_heatmap_<v>` / `predicted_fragment_heatmap_<v>` | `plot_mass_grid` | PNG of the matching grid |

`<v>` is a key of `MASS_GRID_VARIANTS`, fixed in the workflow rather than job
config (no job has needed another): `5ppm` (geometric 5 ppm bins folded to
<= 8192 cells per axis, log colour) and `10da` (linear 10 Da bins, linear
colour capped at the 99th percentile of occupied cells). Both builders take
the same `--ppm`/`--bin-da` flag, so `MassGridVariant.binning` is that flag's
stem and one command template covers both kinds.

The scripts live in necromerge2's `scripts/` (`mass_grid.py`,
`mass_grid_heatmap.py`, `predicted_fragment_grid.py`,
`peptide_length_histogram.py`); the SAGE grid comes from `git/sage`'s
`dump_fragment_index` binary, sourced by `source_dump_fragment_index_binary`
like `dump_peptides`. As with the other `scripts/*.py` rules, node identity
covers the command string, not script contents: after editing a script,
`--invalidate` the affected labels.

## Why the Koina builder computes fragment m/z itself

The shared cache under `mscaches/fragment_intensity/` holds only
`predicted_intensity` and `annotation_id`. `git/featureprediction`'s
`sort_fragment_cache_by_mz` can add a sorted `mz` column, but
`PredictionCache.validate()` rejects a 3-column cache ("Unexpected mmappet
columns"), so the live cache the fill appends to never has one. The builder
therefore runs in `venvs/featureprediction` and derives m/z from each
peptide's residue masses and the slot's (kind, index, charge) with
`feature_prediction.fragment_mz`'s SAGE-ported tables, in a numba kernel
matching `sort_blocks_by_mz`'s formula.

Checked against the old 3-column cache's stored `mz` on 8,468,339
fragments from 200,000 random slots: max difference 1.2e-4 Da / 0.06 ppm,
i.e. float32 rounding of the stored values. Grids built this way vs the
earlier ones built from stored `mz`: same 357,211,338 pairs; 9,012 (5 ppm)
and 395 (10 Da) land in an adjacent cell, from that rounding at bin edges.
The SAGE-index grids from the pipeline are cell-for-cell identical to the
hand-run ones.

Requesting a Koina grid pulls in `export_fragment_intensity_for_sage`, fill
included -- cheap against a warm cache (160.6 s on F9477, mostly the
lookup), a real Koina fill against a cold one.

## Measured (F9477, human FASTA, f9477_best digest settings)

15 nodes from scratch in 283 s with `-call` (dump_peptides 26.9 s, export
160.6 s, each Koina grid ~73 s of which ~40 s is sequence parsing, each
SAGE grid ~10 s plus index build, heatmaps 2-24 s).

## Reading the pictures

What they showed on F9477:
- Per cell, density at a fixed fragment mass falls with precursor mass in
  step with peptide count (x0.49 vs x0.485 between 1100 and 3600 Da on the
  5 Da grid); per whole column it rises (+5.8 % from 2400-2600 to
  2900-3100 Da) because heavier peptides carry more fragments. Both agree
  with the length distribution.
- At sub-Da linear bins, peptide masses form a ~1 Da comb that drifts
  through the bin grid with the mass defect (~1000 Da beat), which shows up
  as broad light/dark patches when the image is downsampled. Use >= 5 Da
  bins to judge broad density.
- The Koina grid shows three wedges, one per precursor charge, ending at
  1001 / 1334 / 2000 m/z.
