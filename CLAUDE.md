# monolixPipelineR

R package for fitting nonlinear mixed effects models in Monolix (through lixoftConnectors) with a
four-step workflow: (1) two-phase initial guess sampling (IGS) to check identifiability,
(2) IGS again with some population parameters fixed, (3) covariate model building (SCM),
(4) paired comparison of candidate treatment-effect models. Built from the functions in the
SIV_HIV_rebound_kinetics project, generalized to any long-format dataset and Mlxtran model.

## Rules
- **Never edit, move, delete or write into the SIV_HIV_rebound_kinetics project** (scripts,
  models, data, results). Reading it is fine. Only files in this repo may be changed.
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
  `mlxtran.R`, `model_spec.R`, `run_log.R`, `scenario.R`, `schema.R` are new.
- `inst/rmarkdown/templates/monolix_report/`: HTML report template (not generalized yet).
- `tests/testthat/`: tests run without Monolix by mocking lixoftConnectors calls.
- `archive/preliminary/`: an older, superseded version of the pipeline, kept for reference.
- `TODO.md`: open items, including checks that need a working Monolix install.
- `HANDOFF.md` and `dev/monolix_checks.R`: the current round of Monolix checks.

## Conventions
- Settings live in `analysis.yaml` (see `config_defaults()`, `draft_config()`). `~` = not
  decided yet; `{}` = none. Run logs (`run_log.yaml`) are written by the pipeline, never edited.
- Column names in every table the package makes: `parameter`, `run`, `estimate`, `rse`, `ofv`,
  `aic`, `bic`, `bicc`, `initial_value`, `model`, `profiled_parameter`, `profiled_value`.
  Tables passed to Monolix keep Monolix's names (`name`, `initialValue`).
- Monolix's model-building criterion value `"BIC"` means BICc.
- Error-model parameters are estimated but not sampled; they start at Monolix's defaults
  (a = 1, b = 0.3, c = 1) unless `model.error_initial` overrides them.
- The Monolix seed is fixed in code (`MONOLIX_SEED`); only the R sampling seed is configurable.

## Testing
- `devtools::test()` runs everything; Monolix-dependent tests use mocks.
- `init_monolix()` reads the install path from the `MONOLIX_PATH` environment variable
  (e.g. set in `~/.Renviron`).
