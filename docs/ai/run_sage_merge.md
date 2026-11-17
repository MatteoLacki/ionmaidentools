# `run_sage`: one rule, optional predictions as mixed Node/`None` inputs

`run_sage` takes `predicted_rt`, `predicted_iim` and `speclib_fragments` as optional
`Node | None` inputs (necroflow's mixed Node/value inputs, `docs/rules.md`), so one rule
covers the first pass, modes 1/2 and the final pass. A workflow branch that has no
prediction assigns `P.predicted_rt = None` directly: `Pipeline.__setattr__` keeps a
non-`Node` value as a plain attribute, with no DAG edge.

`_run_sage_command` is a Python command callback, not a string template, because a
template cannot include a flag conditionally. It adds each `--predicted-*`/SPECLIB flag
when `args.inputs.<name> is not None`: `args.inputs` keeps an unresolved `None` as `None`
and resolves only managed Nodes to paths.

`sage_config` must already carry `rt_tol_sec`+`rt_sigma_sec` and/or
`mobility_tol`+`iim_sigma` for whichever prediction is passed (`update_sage_config`); SAGE's
`Input::build()` enforces the same.
