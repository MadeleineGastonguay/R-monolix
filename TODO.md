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

## Check in the next Monolix round (step 5 code, tested only on fake results)
- `save_igs_diagnostics()` on a real IGS result: every output writes, and the shrinkage
  parameter names Monolix uses map onto the `omega_` rows of `identifiability_summary()`.
- `plot_individual_fits()` guesses the observation column of predictions.txt (the one column
  that isn't id, time or a prediction); check it on real output, including censored data.
- `compare_models()` reads Wald p-values from `getTests()$wald` by the `p.value` column (as the
  SIV script did); confirm the column name.
- Whether Monolix also rejects structural parameters named `b` or `c` with a constant error
  model. `validate_config()` rejects a, b and c to be safe.

## Report (rendered with Monolix 2024R1 on the synthetic data; see dev/monolix_check_results.md)
- Simulx simulation of the model variables: not checked for models with doses or regressors
  (the synthetic model has neither); the imported simulation group's treatments/regressors
  are left as Simulx sets them.
- Several observation types: the prediction plots need one observation column in
  predictions.txt and are skipped otherwise. Check what Monolix writes with two observations.
- Report sections from the pipeline's own results (IGS ranges, SCM table, model comparison)
  aren't in the template yet; it reports one Monolix project.

## Decide later
- `get_narrow_ranges()` now always keeps the best run. The SIV rule (cumulative weight
  <= 0.95) returned no runs when the best run alone carried more than 95% of the weight.

## Covariate search options (method: none, SCM, COSSAC or both)
- Check on Monolix: `covariate_model_building(method = "COSSAC")` runs, and
  `getModelBuildingResults()` from COSSAC has the same shape as SCM's (one entry with
  `bestModel = TRUE`), which `resolve_candidates()` assumes. `strategy` is now passed in lower
  case ("scm", "cossac") as `?getModelBuildingSettings` lists it; "SCM" also worked in check 9.
- `covariate_model_building()` sets `useSambaBeforeCossac`, which isn't among the setting names
  Monolix 2024R1 returned (check 8). Check it's ignored, or remove it.
- The driver (step 7) must give SCM and COSSAC their own project files so their ModelBuilding
  folders and final_model.mlxtran don't collide.
