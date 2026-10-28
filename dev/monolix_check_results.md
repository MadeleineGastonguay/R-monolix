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
