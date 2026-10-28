# To do later

## Settings to add to analysis.yaml
- **Importance sampling sample size.** Left out because the lixoftConnectors arguments
  couldn't be checked without Monolix; Monolix's default is used for now.
- (Done: the `scenario` section now requires `conditional_mode: true`, and the pipeline reads
  conditional-mode individual estimates.)

## Check on the computer where Monolix works
Results are in `dev/monolix_check_results.md`; all checks 0–12 are done. Still open:
- `create_project()` stopping on a rejected model has only been tested with mocks; run it once
  through `setup_project()` with the `a`-named model from section 12.

## Found in the Monolix checks
- Monolix rejects a structural parameter named like an error parameter (`a`, `b`, `c`, or
  `a1`, `b1` … with several observations). Catch this in `validate_config()` /
  `draft_config()` before Monolix is called.
- `getPopulationParameterEstimationSettings()` returns lowercase names (`nbburningiterations`
  …), while the setter takes the camelCase ones. Use lowercase names in any code that reads
  settings back (e.g. the run log).
- `init_monolix()` drops the theophylline estimates from `check_installation()`, so nothing
  shows the demo ran. Print them, or return them invisibly.
- `model.error_initial` keys: with several observations Monolix names error parameters
  `a1`, `b1`, `c1`, `a2`, … (letter plus observation number). Check the keys against these
  names in `validate_config()`.
- Write `scm_truncated_model()` against `tests/testthat/fixtures/scm_results_synthetic.rds`:
  `tested` is the pair added at each step, `individualModels` the full covariate model, one
  entry has `bestModel = TRUE`.
- `initial_guess_sampling()`: an error in one optional record step (condition number, table)
  skips the record steps after it for that run. Consider a separate `tryCatch` per record step.
- In a new project, standard errors and log-likelihood are off by default; the package
  turns them on through `scenario`. Keep in mind when comparing with the SIV runs.
