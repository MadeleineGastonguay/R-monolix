# Monolix check results

Run 2026-10-27 on Maddie's iMac: macOS 14.7.6, R 4.4.0 x86_64 (RStudio), MonolixSuite2024R1,
lixoftConnectors 2024.1, package commit `569dbd5`. Output was captured with `sink()` in two
rounds: all sections, then sections 0, 3, 4 and 12 again after the fixes below. `sink()` only
records printed output, so the first round has no errors, warnings or messages; the second
round prints them for section 12.

| # | Check | Result |
|---|---|---|
| 0 | Setup | Pass, with notes |
| 1 | Unit tests | Pass: 141 passed, 0 failed, 0 skipped (147 after this round's fixes) |
| 2 | Default scenario | Recorded |
| 3 | `setup_project()` | Pass (second round confirmed the SAEM settings) |
| 4 | 3-run IGS | **Failed, fixed, pass on rerun**: condition-number step errored in every run |
| 5 | Same seed twice | Pass: `all.equal` TRUE |
| 6 | Linearization | Pass |
| 7 | Run log | Pass: monolix `"2024R1"` |
| 8 | Model-building settings | Recorded |
| 9 | Small SCM, fixture | Pass, fixtures saved |
| 10 | Report, `plots: false` | Pass: rendered, no error chunks |
| 11 | Error-parameter names | Recorded: numbered by observation |
| 12 | Structural parameter `a` | Recorded: Monolix rejects it; package now stops clearly |

## 0. Setup
- `R.version$arch` is x86_64. `devtools::load_all()` failed at first because pkgload 1.3.4 is too
  old for rlang 1.1.6 (`env_unlock()` is defunct). Fixed by updating pkgload to 1.5.3.
- `MONOLIX_PATH` isn't set on this machine, so the path is passed to `init_monolix()` directly
  in the check script.
- `init_monolix()` ran without error, but the theophylline estimates aren't printed:
  `init_monolix()` returns `invisible(TRUE)` and drops the value from `check_installation()`.
  The "Expect" note was wrong, not the install.
- R processes started by Claude Code on this Mac are killed by macOS when lixoftConnectors
  loads Monolix ("SIGKILL (Code Signature Invalid)"), so the checks must be run in your own R
  session.

## 2. Default scenario of a new project
`getScenario()` right after `newProject()`:

| Task | Default |
|---|---|
| populationParameterEstimation | TRUE |
| conditionalDistributionSampling | TRUE |
| conditionalModeEstimation | TRUE |
| standardErrorEstimation | FALSE |
| logLikelihoodEstimation | FALSE |
| plots | TRUE |
| linearization | FALSE |

So the SIV runs had conditional mode, conditional distribution and plots on, unless their
scripts changed them. Standard errors and log-likelihood are off by default, so the SIV
scripts must have turned those on themselves.

## 3. Project setup from the config
- Scenario after `setup_project()`: every task TRUE except plots; linearization FALSE. As expected.
- SAEM settings: **pass (second round).** `setup_project()` sets them by the camelCase names
  in the help (`nbBurningIterations` …), and they reach Monolix, but
  `getPopulationParameterEstimationSettings()` returns **all-lowercase names**. That's why
  subsetting by camelCase gave `NULL` in the first round. Values read back:
  `nbburningiterations` 5, `nbexploratoryiterations` 50, `nbsmoothingiterations` 20,
  `smoothingalpha` 0.7, `exploratoryautostop` TRUE, `smoothingautostop` TRUE. The other
  defaults, which the package doesn't set: `simulatedannealing` TRUE, `exploratoryalpha` 0,
  `exploratoryinterval` 150, `smoothinginterval` 50, `smoothingratio` 0.1, `variability`
  "none", `tauomega` 0.95, `tauerrormodel` 0.95, `nboptimizationiterations` 20,
  `optimizationtolerance` 1e-04. Any package code that reads these settings back must use the
  lowercase names.
- Project settings: the seed field is `seed` and holds 123456. Other fields: `directory`,
  `exportresults`, `dataandmodelnexttoproject`, `nbsimulations` (500), `grid` (250), `project`.
- Distributions: all logNormal.
- Population parameters: `lambda_pop` … `m_pop`, `omega_*` and `a`, all MLE.

## 4. Initial guess sampling (3 runs)
- `error_starts()`: `a = 1`. Error parameter `a` starts at 1 with MLE in every run.
- `tabestimates`: columns `parameter, estimate, rse, run, ofv, aic, bic, bicc`; RSEs filled in
  (stochastic approximation); 3 runs.
- Shrinkage: one row per parameter per run, all three kinds filled in.
- `indivFit`: the predictions.txt columns, including `indivPred_mode`.
- **Fail:** every run printed `ERROR for i : infinite or missing values in 'x'`. With this short
  fit, `m_pop`'s standard error is `nan`, so its row and column of the correlation matrix are
  `nan` and `eigen()` errors. The error came after the estimates, shrinkage and predictions were
  stored, so those are complete. But `condition_number` and `table_predictions` came back NULL:
  the prediction table is read after the condition number and was never reached.
- **Fix:** new helper `condition_number_of()` (`R/scenario.R`) returns `NA` and prints which
  parameters have non-finite entries, instead of erroring. Test added in
  `tests/testthat/test-fitting.R`.
- **Rerun (second round):** no errors. Each run printed "Condition number not computed",
  `condition_number` is `NA` for runs 1–3, and every other result is the same as in the first
  round. `table_predictions` is now filled in: tables.txt has 28 columns, `id`, `time`, and
  `T/V/E/I/P` with `_mean`, `_mode`, `_popPred`, `_popPred_medianCOV` and `_SAEM`, plus `run`.
  The note listed every parameter, because one nan row and column touches every row; it now
  names only parameters whose own diagonal entry is nan (`m_pop` here).
- `pvals` isn't in the output: no covariates are in the model, so `getTests()$wald` is NULL
  and the element is dropped. Expected.

## 5. Same seed, same results
`all.equal(...)` is TRUE: the fixed Monolix seed makes runs reproducible.

## 6. Linearization
- `apply_scenario()` sets linearization to TRUE.
- 2 runs. RSE and BICc come from linearization (BICc 138.7 against 156.0 with importance
  sampling, same estimates).
- Condition numbers are finite: 3428.6 and 135976.8 (the linearization matrix has no `nan`).

## 7. Run log
`log$versions`: R 4.4.0, monolixPipelineR 0.0.0.9000 (commit 569dbd5), lixoftConnectors
2024.1, **monolix 2024R1**.

## 8. Model-building settings
- Setting names: `covariates, parameters, criterion, nexttoproject, relationships, strategy,
  threshold, useLin`.
- On this project the defaults were `criterion = "LRT"`, `strategy = "cossac"`, `useLin = TRUE`,
  and `threshold$lrt` / `threshold$correlation` of `c(0.01, 0.01)` / `c(0.30, 0.01)`.
- `?getModelBuildingSettings` lists `criterion` as `"BIC"` (documented default) or `"LRT"`, and
  `strategy` as `"cossac"` (default), `"samba"`, `"covsamba"` or `"scm"`. There's no separate
  `"BICc"` value, and section 9 confirms that `"BIC"` gives BICc.

## 9. Small SCM, saved as a test fixture
- Settings: strategy SCM, criterion BIC, `useLin = FALSE`, covariates `tlr7`, `vac`,
  parameters `ta`, `NPs`. All accepted.
- `getModelBuildingResults()` returns 8 entries, each with `BICc`, `LL`, `tested`,
  `individualModels` and `bestModel`:
  - Entry 1 is the base model (`tested = c("", "")`).
  - `tested` is the parameter–covariate pair added at that step.
  - `individualModels` is a logical data frame (rows = parameters, columns = covariates) holding
    the full covariate model at that step: entry 6 has both `ta~vac` and `ta~tlr7`.
  - Exactly one entry has `bestModel = TRUE`: entry 3 (`ta~vac`, BICc 166.7, against 178.7 for
    the base). The second forward step (entries 6–8) didn't improve BICc.
- Saved, both from synthetic data only: `tests/testthat/fixtures/scm_results_synthetic.rds`
  and `scm_results_synthetic_str.txt`.

## 10. Report with plots: false
`generate_report()` rendered `igs_check/report.html` (5.5 MB) with no error chunks. Sections:
Project Details, Project Setup, Predictions (population trends per group, individual
predictions), Results (population and individual estimates), Model Diagnostics, Performance of
Estimation Algorithms. I haven't checked the content by eye.

## 11. Error-parameter names with several observations
The script loaded the first match, `warfarin_PKPDseq_project` (one observation, `y2`), which
gave a constant error model with parameter `a2`. The two-observation demo
`warfarinPKPD_project.mlxtran` shows the pattern:
`y1 = combined1(a1, b1)` (plus `c1`), `y2 = constant(a2)`.
**Error parameters are named by letter plus observation number** (`a1`, `b1`, `c1`, `a2` …),
not `a_y1`. So `model.error_initial` keys for several observations should be `a1`, `b2` etc.,
and starting values by first letter (a/c = 1, b = 0.3) still work.

## 12. Structural parameter called `a`
**No.** Monolix refuses the model (second round):

    [ERROR] Name 'a' is already used as error parameter in the statistical model. Please use another one in the structural model.

`newProject()` returns `FALSE` and raises no R error, so the next call fails with
`[ERROR] No project loaded`. With one observation the error parameters are `a`, `b`, `c`; with
several they're `a1`, `b1` … (section 11), so those names would clash in turn.
**Fix:** `create_project()` now stops right away with "Monolix could not create the project from
<model file>" when `newProject()` returns `FALSE`. Test added in `tests/testthat/test-setup.R`.
An earlier check, in `validate_config()` or `draft_config()`, that structural parameter names
don't clash with error-parameter names is in TODO.md.

## Fixes made in this round
- `condition_number_of()`: returns NA with a note instead of erroring on `nan`
  (`R/scenario.R`, used in `initial_guess_sampling()`). Confirmed on a real run.
- `create_project()`: stops when `newProject()` returns `FALSE` (`R/setup.R`).
- Tests for both; 147 tests pass (mocked, run without Monolix).
- `dev/monolix_checks.R`: section 3 prints all SAEM settings; section 12 prints conditions.
  Your edits are kept: explicit Monolix path, `setwd()` and `sink()`.

The `create_project()` stop is tested with mocks only. Section 12 calls `newProject()`
directly, so it hasn't been run through Monolix yet.

# Round 3: report (sections 1, 10, 13, 14; log `dev/monolix_check_output4.log`, 2026-10-28)

## 1. Unit tests
282 pass, 0 fail, after updating GGally to 2.4.0 (the `test-plots.R` failure came from
GGally 2.2.1 with ggplot2 4.0.3).

## 10 and 13. Reports
Both rendered. Checked in the HTML:
- CV column filled for every omega (names `omega_<param>_CV` as expected).
- `NP_pop` rescaled by the config (1.5e+05); `beta_ta_vac_1` in the treatment-effects table.
- report_check: groups from `tlr7` x `vac` with labels; stratified VPC and prediction
  distribution drawn. igs_check (no config): one group, "No stratification" note.
- Condition-number notes name the non-finite parameters (`m_pop`; in report_check also
  `beta_ta_vac_1`, `omega_p`): expected with 3 short SAEM runs on synthetic data.
- **Bug found:** "Correlation of random effects" showed "Not available: arguments imply
  differing number of rows: 0, 1". The pairs plot kept only columns named after the parameters,
  but the random-effects table calls them `eta_<param>`. Fixed (`prefix = "eta_"`).

Fine-grained fits: `ChartsData/IndividualFits/vl_fits.txt` (and `vl_observations.txt`), 256
times per individual, columns `ID, time, pop, popPred, indivPredMean, indivPredMode`. `pop` uses
the median covariates (equals `popPred_medianCOV` in predictions.txt); `popPred` uses each
individual's covariates. **Bug found:** `read_indfits()` picked `pop`, so the population plot
ignored covariate effects. Fixed to use `popPred`.

## 14. Simulx model variables
Works. 5 variables x 12 individuals x 50 times; ids are the data's ids; lixoftConnectors is back
in Monolix mode afterwards. At each individual's first and last observation time the simulated
values match `<variable>_mode` in tables.txt to within 0.09% (ODE solver), so the EBEs are
paired with the right individuals. Simulx warns at import: "For mlx_EBEs and mlx_Cov, shared
IDs should be ticked ..." and "mlx_PopUncertainSA has not been created because there are NaNs
in the RSE". Neither matters here: the simulation uses the EBEs, which already include
covariate effects, not the covariate or uncertainty elements.

# Round 4: report fixes and step 5 code (sections 13, 15–18; log `dev/monolix_check_output5.log`, 2026-10-28)

## 13. Report fixes
Confirmed on the real render: the population trends plot has a separate curve for vac = 1
(popPred), the random-effects pairs plot is drawn, and no section shows "Not available".

## 15. Step 5 diagnostics on a real IGS result
`save_igs_diagnostics()`: all 15 outputs written (13 PDFs, 3 CSVs). `identifiability_summary()`
puts shrinkage on every `omega_` row (Monolix names them `lambda`, `betas`, ...).
`observation_column()` finds `vl`. ggplot warned that 4 points were removed as missing or out
of range in a plot; check `individual_fits.pdf` for censored observations.

## 16. Wald tests and compare_models()
**Failed:** the base model's IGS ran, but the first run of the ta ~ vac model stopped with
"No sampling range for beta_ta_vac_1". `sampling_ranges()` gives treatment effects no range, so
`draw_initial_values()` needed `set_treatment`. **Fixed:** a `beta_` without `set_treatment` or
a range now starts at 0, as in `set_base()`; the other draws are unchanged. Test added. Rerun
section 16 to check the Wald column names and `compare_models()`.

## 17. Rejected parameter names
- `setup_project()` with a structural parameter `a` now stops with "Monolix could not create
  the project from .../model_param_a.txt" (the create_project() fix, confirmed on Monolix).
- `b` is rejected too ("Name 'b' is already used as error parameter"), even with a constant
  error model.
- `c`: inconclusive. The test model already defines a constant `c = 23`, so Monolix reported
  "Conflicting variable definition for 'c'".

## 18. Covariate search through the package
SCM (8 models) and COSSAC (3 models) both ran through `covariate_model_building()` with
`strategy` in lower case; each has exactly one `bestModel = TRUE`, and COSSAC's entries have
the same fields as SCM's. No warning or message about `useSambaBeforeCossac`. Both chose
`ta ~ vac`; `resolve_candidates()` kept `scm_best` and dropped `scm_delta2` and `cossac_best`
as duplicates, with the message and a warning that only 1 distinct candidate is left. The
best-model reports are in `search_scm/final_model/SCM_final_model_report.html` and
`search_cossac/final_model/COSSAC_final_model_report.html` (the check's expectation said
`<search>/`; corrected).
