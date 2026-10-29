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
- Real-Monolix test runs set up in dev/step8/ (copies of `identifiable`, run with dev/step8/run.R):
  - `cossac` (step 3 checked 2026-10-29): SCM and COSSAC with criterion LRT on vac only. tlr7
    was never tested; no scm_delta2 with LRT; SCM's LRT steps match the thresholds (betas-vac
    added at p = 0.005, m-vac rejected at p = 0.25); COSSAC's results have the same fields as
    SCM's (BICc, LL, tested, individualModels, bestModel) and reached the same model in 3 fits;
    cossac_best = scm_best was fitted once; scm_best was copied from the search. Neither search
    produced backward or rejected entries (Monolix reuses forward fits for backward tests), so
    `scm_results_table()` and `plot_model_building()` on COSSAC's alternating order are still
    untested: needs a search where COSSAC removes or rejects a covariate.
  - `no_omega`: `random_effects_excluded: [ta]`, from step 1 to the final model without stops
    (2 runs per IGS, all decisions filled in): checks only that every step runs with a parameter
    without a random effect. SCM and COSSAC search vac on ta (no random effect) and NPs.
- COSSAC with a parameter without a random effect has only that quick run: Monolix's COSSAC
  can't search such parameters and handles them with SCM at the end, so its results may mix the
  two orders. Check how `cossac_tested_models.csv` reads, then test it on a realistic search.
- Run inst/examples/synthetic end to end through all three checkpoints (small settings already
  in its analysis.yaml), then once with `method: [SCM, COSSAC]` and once with `method: none`.
- Step 3's `relationships` (`include`/`exclude`) goes to Monolix as a data frame of
  `parameters`, `covariates` and `locked` (logical), per the lixoftConnectors help. Check on
  real Monolix that included pairs are in every tested model and excluded pairs in none
  (`scm_tested_models.csv`), and that a logical `locked` column is accepted.
- Step 3 copies `scm_best`/`cossac_best` from the search's `final_model.mlxtran` instead of
  refitting them (`saveProject()` to a new file copies the results folder; confirmed). On real
  Monolix, check that step 4 loads and refits the copied `candidates/scm_best.mlxtran`.
- Step 4 has no per-model cache: if it stops partway, the models already fitted are refitted.
  Consider resuming from the saved <model>_igs.rds files.
- The final model has no run log; after changing step 4, rerun it with rerun = "final".
- Profile likelihoods (`run_profiles()`) on real Monolix: check that turning off standard errors
  through the scenario still gives the importance-sampling OFV, and that saving the narrow
  project as profile.mlxtran leaves step1_igs/narrow/ unchanged.
- Generality: step 1 on Monolix's theophylline demo with only a new analysis.yaml.
- `devtools::check()`.
