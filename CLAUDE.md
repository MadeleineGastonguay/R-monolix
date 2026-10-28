# monolixPipelineR

R package for fitting nonlinear mixed effects models in Monolix (through lixoftConnectors) with a
four-step workflow: (1) two-phase initial guess sampling (IGS) to check identifiability,
(2) IGS again with some population parameters fixed, (3) covariate model building (SCM),
(4) paired comparison of candidate treatment-effect models. Built from the functions in the
SIV_HIV_rebound_kinetics project, generalized to any long-format dataset and Mlxtran model.

## Rules
- **The package is self-contained.** Code, tests, checks and examples use only files in this
  repo. Don't read from, write to or refer to paths in the SIV_HIV_rebound_kinetics project or
  any other folder outside the repo.
- **Never install, update or upgrade an R package without asking first**, including packages
  pulled in as dependencies. Say which package is missing and what it blocks, and offer the
  command.
- **lixoftConnectors only works with the x86_64 build of R**, not arm64. Check
  `R.version$arch` before Monolix work. Packages installed under one architecture aren't
  visible from the other.
- **No real study data in this repo.** Tests use `tests/testthat/fixtures/synthetic_vl_data.txt`.
  Ask before committing anything derived from real data.
- Commit only when asked. Work on the `package-skeleton` branch.

## Layout
- `R/`: package code. `fitting.R`, `setup.R`, `selection.R` and `plots_*.R` are ported from the
  SIV `functions.R` (bug fixes noted at the top of each file). `config.R`, `draft_config.R`,
  `mlxtran.R`, `model_spec.R`, `run_log.R`, `scenario.R`, `schema.R`, `sampling.R`, `summaries.R`,
  `model_building.R`, `diagnostics.R`, `plots_fits.R`, `plots_model_comparison.R`, `plots_theme.R`,
  `report.R`, `simulx.R`, `pipeline.R` are new.
- `R/pipeline.R`: one function per step (`run_step1()` … `run_final_model()`, `run_pipeline()`).
  `inst/scripts/run_pipeline.R` is the driver script; `inst/examples/synthetic/` is an example
  config on the synthetic data. Pipeline tests use a fake Monolix (`helper-pipeline.R`).
- `inst/rmarkdown/templates/monolix_report/`: HTML report template, rendered by `generate_report()`
  in an environment whose parent is the package namespace (it uses package helpers directly).
  Groups, labels and rescaling come from `report.*` and `model.report_transforms` in the config;
  table helpers are in `R/report.R`. Before rendering, model variables are simulated with Simulx
  (`R/simulx.R`), which switches lixoftConnectors to Simulx and back to Monolix.
- `tests/testthat/`: tests run without Monolix by mocking lixoftConnectors calls.
- `archive/preliminary/`: an older, superseded version of the pipeline, kept for reference.
- `TODO.md`: open items, including checks that need a working Monolix install.
- `HANDOFF.md` and `dev/monolix_checks.R`: the current round of Monolix checks.

## Conventions
- Settings live in `analysis.yaml` (see `config_defaults()`, `draft_config()`). `~` = not
  decided yet; `{}` = none. Run logs (`run_log.yaml`) are written by the pipeline, never edited.
- Every plot the package draws uses `theme_report()` (`R/plots_theme.R`; build plots with
  `ggplot_report()`), groups coloured with Dark2. Monolix's own plots in the report keep theirs.
- Plotting and summary functions never call lixoftConnectors (checked by
  `tests/testthat/test-no-monolix.R`), so they work on saved results without Monolix.
- Column names in every table the package makes: `parameter`, `run`, `estimate`, `rse`, `ofv`,
  `aic`, `bic`, `bicc`, `initial_value`, `model`, `profiled_parameter`, `profiled_value`.
  Tables passed to Monolix keep Monolix's names (`name`, `initialValue`).
- Monolix's model-building criterion value `"BIC"` means BICc.
- `step3_model_building.method` is SCM, COSSAC, [SCM, COSSAC] or none. Each search adds
  automatic candidates (`scm_best`, `scm_delta2`; `cossac_best`); `resolve_candidates()` merges
  them with `candidate_models` and fits identical covariate models only once.
- Error-model parameters are estimated but not sampled; they start at Monolix's defaults
  (a = 1, b = 0.3, c = 1) unless `model.error_initial` overrides them.
- The Monolix seed is fixed in code (`MONOLIX_SEED`); only the R sampling seed is configurable.

## Testing
- `devtools::test()` runs everything; Monolix-dependent tests use mocks.
- `init_monolix()` reads the install path from the `MONOLIX_PATH` environment variable
  (e.g. set in `~/.Renviron`).
