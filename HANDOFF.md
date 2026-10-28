# Handoff: Monolix checks done

The Monolix checks in `dev/monolix_checks.R` (sections 0–12) were run against real Monolix
(MonolixSuite2024R1, lixoftConnectors 2024.1, R 4.4.0 x86_64) on 2026-10-27. Everything passes
or is recorded. Two package bugs turned up and are fixed. Work can move on to steps 5–7
(summaries, SCM-based candidates, report, driver). Read `CLAUDE.md` first for the rules.

## Prompt to start the next session
> Read CLAUDE.md, HANDOFF.md and dev/monolix_check_results.md, then continue with the next
> step of the plan.

## Where things are
- `dev/monolix_check_results.md`: per-section results, including the default scenario,
  model-building criterion values, error-parameter names and the SCM results structure.
- `TODO.md`: what's still open and what the checks found.
- `tests/testthat/fixtures/scm_results_synthetic.rds` (and `_str.txt`): real
  `getModelBuildingResults()` output from the synthetic data, for writing `scm_truncated_model()`.

## Bugs fixed (commit 25dce6e)
1. **Condition number in `initial_guess_sampling()`.** When a parameter's standard error is
   `nan`, its row and column of the correlation matrix are `nan`, and `eigen()` errored. That
   also skipped reading `table_predictions` for the run. The new helper `condition_number_of()`
   (`R/scenario.R`) returns `NA` and prints which parameter caused it. Confirmed on a real run.
2. **`create_project()` ignored a failed `newProject()`.** Monolix reports a rejected model
   (e.g. a structural parameter named `a`) as an `[ERROR]` message and returns `FALSE`, with no
   R error, so the pipeline failed later with "No project loaded". It now stops right away.
   Tested with mocks only.

Both have tests; 147 tests pass.

## Findings to keep in mind for steps 5–7
- Error parameters are named `a`, `b`, `c` with one observation, and `a1`, `b1`, `c1`, `a2` …
  with several. Monolix rejects structural parameters with those names.
- `getPopulationParameterEstimationSettings()` returns lowercase names (`nbburningiterations`);
  the setter takes the camelCase ones.
- Model-building `criterion` is `"BIC"` (gives BICc) or `"LRT"`.
- In `getModelBuildingResults()`, `tested` is the pair added at each step, `individualModels`
  is the full covariate model, and one entry has `bestModel = TRUE`.

## Running Monolix code on the Mac (Maddie's iMac)
- macOS kills any R process that Claude Code starts once lixoftConnectors loads Monolix
  ("SIGKILL (Code Signature Invalid)"), so Claude can't run Monolix code there. Maddie runs it
  in her own R session (RStudio), from the repo root, and logs output with
  `sink("dev/monolix_check_output.log", split = TRUE)`. Claude reads the log. Logs aren't
  committed.
- `sink()` only captures printed output. Print errors, warnings and messages explicitly
  (see `show_conditions()` in section 12 of the check script).
- `MONOLIX_PATH` isn't set there; pass the path to `init_monolix()` or add it to `~/.Renviron`.
- pkgload was updated to 1.5.3 in the x86_64 R there, because 1.3.4 broke
  `devtools::load_all()` with rlang 1.1.6.
