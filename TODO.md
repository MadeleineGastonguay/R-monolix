# To do

Checked on real Monolix (MonolixSuite2024R1) so far: every step from step 1 to the final model,
on a synthetic viral-load model (with SCM and COSSAC, BICc and LRT, locked and excluded
covariate relationships, a parameter without a random effect, profile likelihoods) and on
Monolix's phenobarbital demo (a PK model with doses, a random effect removed at checkpoint 1,
step 4 with one model). The tests use a fake Monolix, so anything new that talks to Monolix
needs a real run too.

## Not yet checked on real Monolix
- **`step3_model_building.method: none`** (only hand-written candidates) and
  **`[SCM, COSSAC]` together on a search where they disagree**.
- **COSSAC removing or rejecting a covariate.** No run so far produced backward or rejected
  entries, so `scm_results_table()` and `plot_model_building()` haven't seen COSSAC's
  alternating order. COSSAC with a parameter without a random effect (which Monolix handles with
  SCM at the end) has only had a quick run.
- **Logit-normal limits** (`model.distributions: {param: {distribution: logitnormal, limits: [min, max]}}`):
  `set_logit_limits()` calls `setIndividualLogitLimits()`, and moves the starting value inside
  the limits if Monolix doesn't take them at first. Also check that Monolix accepts
  `probitNormal`, which the connectors' help pages don't list.
- **Models with regressors** in the report's Simulx simulation (models with doses are checked).
- **Several observation types.** The report's prediction plots need one observation column in
  `predictions.txt` and are skipped otherwise.

## Improvements
- **Step 4 resumes from scratch.** If it stops partway, models already fitted are fitted
  again; it could load their saved `<model>_igs.rds`.
- **The final model has no run log**, so it isn't checked against changed settings. After
  changing step 4, rerun it (`rerun = "final"`).
- **Importance-sampling sample size** isn't a setting in `analysis.yaml` yet; Monolix's default
  is used.
- **The report covers one Monolix project.** The pipeline's own results (IGS ranges, the
  covariate search, the model comparison) aren't in it.
- `covariate_model_building()` sets `useSambaBeforeCossac`, which isn't a Monolix 2024R1
  setting (it's ignored without a warning) and could be removed.
