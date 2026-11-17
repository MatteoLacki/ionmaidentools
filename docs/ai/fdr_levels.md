# FDR levels and run recording inside `run_sage` / `mokapot`

The last step of both commands counts the node's confident spectra, peptidoforms and
ions at 1% FDR, each level from its own target-decoy competition, into a declared output,
and records the run in `$ION_MAIDEN_REGRESSION_DB` when that is set.

| rule | step | output (label) |
|---|---|---|
| `run_sage` | `venvs/common/bin/sage-fdr-levels` | `results.sage.fdr_levels.tsv` (`sage_fdr_levels`; the recalibration pass: `filtered_sage_fdr_levels`) |
| `mokapot` | `venvs/common/bin/mokapot-fdr-levels` | `mokapot.fdr_levels.tsv` (`mokapot_fdr_levels`) |

Definitions: git/searchops `docs/ai/ion_level_fdr.md`. Database: git/searchops
`docs/ai/run_db.md`.

`mokapot` takes `sage_results_tsv` as a second input, because mokapot's tables carry no
charge and an ion is `(peptide, charge)`. It reads the decoy PSMs that `--keep_decoys`
writes (`mokapot.decoy.psms.txt`, undeclared, in the node directory).

`_fdr_levels_record_flags` passes the run database what it cannot find itself: the
node directory, the tool command (SAGE's call, or the pin adapter + mokapot call), and
every input path (`--parent`). The recorder walks those inputs' `.rip/dependencies.toml`
up to `symlink_BrukerD` for the configs and the dataset; the node's own metadata is
written only after its command finishes.

Why inside the tool commands rather than a separate rule: the counts and the record then
exist for every SAGE and mokapot node that ever runs, in its own directory, and a record
is made exactly when the tool ran. A job still has to request the labels to get them
copied to its results directory, so every job requesting SAGE or mokapot outputs
requests them (`sage_fdr_levels`, `mokapot_fdr_levels`, and
`filtered_sage_fdr_levels` where the job requests `filtered_sage_results_tsv`). The cost:
editing either command callback changes every SAGE or mokapot node's identity.

Counting is DuckDB SQL over the result tables and takes about 1 s per node on F9477
(git/searchops `fdr_levels.py`, `ion_level_fdr.py`).

The `sage_summarize` rule (`sage-summarize-raw`: PSMs, peptides, PSM-collapsed ions and
proteins at peptide-level FDR) counts differently and is requested by some jobs.
