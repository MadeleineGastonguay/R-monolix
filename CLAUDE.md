# monolixPipelineR

R package for fitting nonlinear mixed effects models in Monolix (through lixoftConnectors) with a
four-step workflow: (1) two-phase initial guess sampling (IGS) to check identifiability,
(2) IGS again with some population parameters fixed or random effects removed, (3) covariate
model building (SCM and/or COSSAC) and one fit per candidate model, (4) paired comparison of the
candidates; then the final model and its report. Works with any long-format dataset and Mlxtran
structural model. README.md is the user guide.

## Rules
- **The package is self-contained.** Code, tests, checks and examples use only files in this
  repo. Don't read from, write to or refer to paths in any folder outside the repo.
- **Keep everything general.** Help pages, messages, examples and the README use generic names
  (`param1`, `cov1`) or the user's own model and covariate names, never those of a particular
  study. The `# Ported from ...` comments at the top of the ported R files stay: they record
  where that code came from.
- **Never install, update or upgrade an R package without asking first**, including packages
  pulled in as dependencies. Say which package is missing and what it blocks, and offer the
  command.
- **lixoftConnectors only works with the x86_64 build of R**, not arm64. Check
  `R.version$arch` before Monolix work. Packages installed under one architecture aren't
  visible from the other.
- **No real study data in this repo.** Tests use `tests/testthat/fixtures/synthetic_vl_data.txt`.
  Ask before committing anything derived from real data (test runs on other data go in the
  gitignored `dev/step8/`).
- Commit only when asked. Work on `main`.

## Layout
- `R/`: package code. `fitting.R`, `setup.R`, `selection.R`, `plots_igs.R` and `plots_model_building.R` are ported from an
  earlier analysis's scripts (changes since porting noted at the top of each file). `config.R`, `draft_config.R`,
  `mlxtran.R`, `model_spec.R`, `run_log.R`, `history.R`, `scenario.R`, `schema.R`, `sampling.R`, `summaries.R`,
  `model_building.R`, `diagnostics.R`, `plots_fits.R`, `plots_model_comparison.R`, `plots_profile.R`, `plots_theme.R`,
  `report.R`, `simulx.R`, `pipeline.R`, `pages.R` are new.
- `R/pipeline.R`: one function per step (`run_step1()` … `run_final_model()`, `run_pipeline()`).
  `run_profiles()` runs the optional profile likelihoods after step 1, in a copy of the narrow
  project, with its own run log; it never refits or changes step 1. Step 2 runs when checkpoint 1 fixes parameters (`decisions.fixed`)
  or removes random effects (`decisions.random_effects_removed`, applied from step 2 on;
  `model.random_effects_excluded` is for ones left out from step 1); see `step2_runs()`. When
  step 2 runs, checkpoint 1b (`decisions.fixed_confirmed: true`) must be recorded before step 3. Checkpoint 2b
  (`step4_model_comparison.compare`: `all` or a list) comes after step 3's single candidate fits
  (`single_fits.csv`); `step4_models()` applies it in step 4, so step 3 always fits every candidate.
  `inst/scripts/run_pipeline.R` is the driver script; `inst/examples/synthetic/` is an example
  config on the synthetic data. Pipeline tests use a fake Monolix (`helper-pipeline.R`).
- `inst/rmarkdown/templates/results_page/` and `R/pages.R`: one HTML page per step, from saved
  results (blocks of figures, tables and notes): `render_igs_figures()` (steps 1 and 2, every figure
  of `igs_figures()`), `render_model_building_page()` (step 3: the search, then the candidate fits)
  and `render_comparison_page()` (step 4). The pipeline writes each at the end of its step; the
  PDFs and CSVs are the same figures and tables.
- `inst/rmarkdown/templates/monolix_report/`: HTML report template, rendered by `generate_report()`
  in an environment whose parent is the package namespace (it uses package helpers directly).
  Groups, labels and rescaling come from `report.*` and `model.report_transforms` in the config;
  table helpers are in `R/report.R`. Before rendering, model variables are simulated with Simulx
  (`R/simulx.R`), which switches lixoftConnectors to Simulx and back to Monolix.
- `tests/testthat/`: tests run without Monolix by mocking lixoftConnectors calls.
- `TODO.md`: open items, including what still needs a check on real Monolix.
- `dev/` (not in the built package): `monolix_checks.R` and `monolix_check_results.md`, the first
  round of real-Monolix checks; `step8/` (gitignored) holds real-Monolix test runs, each a folder
  with its own `analysis.yaml`, run with `dev/step8/run.R`.

## Conventions
- Settings live in `analysis.yaml` (see `config_defaults()`, `draft_config()`). `~` = not
  decided yet; `{}` = none. Run logs (`run_log.yaml`) are written by the pipeline, never edited.
  Rerunning a step deletes its old results; `<out_folder>/history.yaml` (`R/history.R`,
  `read_history()`) keeps one entry per finished step: what changed, the decisions, the main result.
- IGS plots about the consistency of the estimates (convergence, normalized boxplots, CV, paired
  estimates, correlations between estimates) show only the best-supported runs (`top_runs()`, the
  `p_threshold` BICc weight); the others show every run, coloured by BICc. Each caption says which
  runs it shows and what its lines and marks mean. The identifiability summary's parameter
  correlations also use the top runs.
- Every plot the package draws uses `theme_report()` (`R/plots_theme.R`; build plots with
  `ggplot_report()`), groups coloured with Dark2. Monolix's own plots in the report keep theirs.
- Plotting and summary functions never call lixoftConnectors (checked by
  `tests/testthat/test-no-monolix.R`), so they work on saved results without Monolix.
- Column names in every table the package makes: `parameter`, `run`, `estimate`, `rse`, `ofv`,
  `aic`, `bic`, `bicc`, `initial_value`, `model`, `profiled_parameter`, `profiled_value`.
  Tables passed to Monolix keep Monolix's names (`name`, `initialValue`).
- Monolix's model-building criterion value `"BIC"` means BICc.
- `step3_model_building.method` is SCM, COSSAC, [SCM, COSSAC] or none. Each search adds
  automatic candidates (`scm_best`, plus `scm_delta2` with criterion BIC only; `cossac_best`); `resolve_candidates()` merges
  them with `candidate_models` and fits identical covariate models only once.
- Error-model parameters are estimated but not sampled; they start at Monolix's defaults
  (a = 1, b = 0.3, c = 1) unless `model.error_initial` overrides them.
- The Monolix seed is fixed in code (`MONOLIX_SEED`); only the R sampling seed is configurable.

## Testing
- `devtools::test()` runs everything; Monolix-dependent tests use mocks, and the pipeline tests a
  fake Monolix (`helper-pipeline.R`). They check the package's logic, not Monolix's behaviour:
  anything new that talks to Monolix needs a real run too.
- `R CMD check` passes with no errors or warnings; keep it that way (every exported function
  documents all its arguments and its return value).
- `init_monolix()` reads the install path from the `MONOLIX_PATH` environment variable
  (e.g. set in `~/.Renviron`).
- On the Mac where this was developed, macOS kills an R process started by Claude Code once
  lixoftConnectors loads Monolix. Monolix code is run in the user's own R session (RStudio);
  `dev/step8/run.R` logs everything to `dev/step8/<run>_<time>.log` for Claude to read.
  Code that doesn't need Monolix (tests, plots from saved results) runs fine from Claude.
