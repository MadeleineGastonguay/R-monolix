# To do later

## Settings to add to analysis.yaml
- **Importance sampling sample size.** Left out because the lixoftConnectors arguments
  couldn't be checked without Monolix; Monolix's default is used for now.
- (Done: the `scenario` section now requires `conditional_mode: true`, and the pipeline reads
  conditional-mode individual estimates.)

## Check on the computer where Monolix works
- `setup_project()` / `create_project()` with `distributions` and `estimation` from a config.
- `set_project_settings()` with SAEM and task settings from the config.
- `scenario`: what `getScenario()` shows right after `newProject()`, i.e. which tasks the SIV
  runs actually had. The package defaults to every task on except `plots`; check that reports
  still render with `plots: false`, and that `linearization: true` runs end to end.
- The Monolix version recorded in `run_log.yaml` (`getLixoftConnectorsState(quietly = TRUE)$version`).
- How Monolix names error-model parameters when there are several observation columns
  (`model.error_initial` keys aren't checked against the error model until then; the default
  starting values go by first letter: a/c = 1, b = 0.3).
- Whether a structural-model parameter can be named `a`, `b` or `c` (not on the reserved-keyword
  list, but it would clash with the error-model parameter names).
