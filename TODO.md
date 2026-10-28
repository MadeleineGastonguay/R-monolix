# To do later

## Settings to add to analysis.yaml
- **Importance sampling sample size.** Left out because the lixoftConnectors arguments
  couldn't be checked without Monolix; Monolix's default is used for now.
- (Done: the `scenario` section now requires `conditional_mode: true`, and the pipeline reads
  conditional-mode individual estimates.)

## Check on the computer where Monolix works
Results are in `dev/monolix_check_results.md`; all checks 0–15, 17 and 18 are done.

## Found in the Monolix checks
- (Done: `validate_config()` catches structural parameters named like error parameters, and
  checks `model.error_initial` keys against `a1`, `b1` … names with several observations;
  `init_monolix()` prints and returns the demo estimates; each optional record step in
  `initial_guess_sampling()` has its own error handler; `scm_truncated_model()` is written
  against the SCM fixture.)
- `getPopulationParameterEstimationSettings()` returns lowercase names (`nbburningiterations`
  …), while the setter takes the camelCase ones. Use lowercase names in any code that reads
  settings back (e.g. the run log).
- In a new project, standard errors and log-likelihood are off by default; the package
  turns them on through `scenario`. Keep in mind when comparing with the SIV runs.

## Check in the next Monolix round
- Rerun section 16 of `dev/monolix_checks.R` (Wald tests and `compare_models()`): it stopped at
  the first IGS run of the model with a treatment effect ("No sampling range for
  beta_ta_vac_1", now fixed). Confirms the column names of `getTests()$wald`.
- Rerun sections 15 and 18: `plot_individual_fits()` now draws every IGS run (grey) with the
  best run in colour, and censored observations as open circles (the IGS adds a `censored`
  column from the data); the model-building reports get the config's groups (`config =`).
- Whether Monolix rejects a structural parameter named `c`: the section 17 test model already
  defines a constant `c`, so it failed for that reason. `a` and `b` are rejected.

## Report (rendered with Monolix 2024R1 on the synthetic data; see dev/monolix_check_results.md)
- Simulx simulation of the model variables: not checked for models with doses or regressors
  (the synthetic model has neither); the imported simulation group's treatments/regressors
  are left as Simulx sets them.
- Several observation types: the prediction plots need one observation column in
  predictions.txt and are skipped otherwise. Check what Monolix writes with two observations.
- Report sections from the pipeline's own results (IGS ranges, SCM table, model comparison)
  aren't in the template yet; it reports one Monolix project.


## Covariate search options (method: none, SCM, COSSAC or both)
- (Done, check 18: SCM and COSSAC run through `covariate_model_building()` with lower-case
  `strategy`; COSSAC's results have SCM's shape with one best model; `useSambaBeforeCossac`
  gives no warning or message. It isn't a Monolix 2024R1 setting, so it could be removed.)
- (Done: the driver gives SCM and COSSAC their own projects, search/scm.mlxtran and
  search/cossac.mlxtran.)

## Step 8: run the pipeline for real (on the Monolix computer)
The driver (R/pipeline.R, inst/scripts/run_pipeline.R) is tested only with a fake Monolix.
- Run inst/examples/synthetic end to end through all three checkpoints (small settings already
  in its analysis.yaml), then once with `method: [SCM, COSSAC]` and once with `method: none`.
- Check that `saveProject()` to a new file moves the results folder with it, as
  `save_project_as()` assumes (the IGS reads predictions.txt from the project's folder).
- Step 3's `relationships` is passed to Monolix as written in the YAML; check the format
  `getModelBuildingSettings()$relationships` expects.
- Step 4 has no per-model cache: if it stops partway, the models already fitted are refitted.
  Consider resuming from the saved <model>_igs.rds files.
- The final model has no run log; after changing step 4, rerun it with rerun = "final".
- Generality: step 1 on Monolix's theophylline demo with only a new analysis.yaml.
- `devtools::check()`.
