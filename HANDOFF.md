# Handoff: Monolix checks

Steps 1-4 of the plan were built on a computer where lixoftConnectors doesn't work, so every
Monolix call has only been tested with mocks. This round checks them against real Monolix on a
computer where it works, before steps 5-7 (summaries, SCM-based candidates, report, driver) are
built on top. Read `CLAUDE.md` first for the rules.

## Prompt to start the session
> Read CLAUDE.md and HANDOFF.md, then help me run the Monolix checks in dev/monolix_checks.R
> section by section and record the results.

## Before running
1. `git pull` on the `package-skeleton` branch.
2. Open R and check `R.version$arch` is `x86_64`.
3. Check these packages are installed (don't install without asking): devtools, testthat,
   roxygen2, yaml, gt, rmarkdown, lixoftConnectors, plus the Imports in `DESCRIPTION`.
   `devtools::load_all()` lists any that are missing.
4. Set `MONOLIX_PATH` in `~/.Renviron` to the Monolix install folder, then restart R.
5. In `dev/monolix_checks.R`, set `siv_root` to where SIV_HIV_rebound_kinetics is on this computer.

## The checks (`dev/monolix_checks.R`)
Each section prints what to look for after `# Expect:` or `# Record:`.

| # | Check | Why |
|---|---|---|
| 1 | Unit tests with lixoftConnectors installed | mocked tests that were skipped or untested before |
| 2 | Default scenario of a new project | what the SIV runs used for tasks the scripts didn't set |
| 3 | `setup_project()` from a config | scenario, SAEM settings and seed reach Monolix |
| 4 | 3-run `initial_guess_sampling()` | results read correctly: RSE, BICc, shrinkage, condition number, predictions |
| 5 | Same seed twice | the fixed Monolix seed makes runs reproducible |
| 6 | `linearization: true` | results read from the linearization method |
| 7 | Run log | Monolix version recorded |
| 8 | Model-building settings | accepted `criterion` values ("BIC" = BICc) |
| 9 | Real `getModelBuildingResults()` | structure needed for `scm_truncated_model()` next; **ask before saving it to fixtures** |
| 10 | Report with `plots: false` | report renders; template not generalized yet, so note failures |
| 11 | Error-parameter names with 2 observations (optional) | so `model.error_initial` keys can be checked |
| 12 | Structural parameter called `a` (optional) | whether it clashes with the error parameter |

The script writes only to a temporary folder (and, in section 9, to `tests/testthat/fixtures/`
if Maddie agrees). Section 9 only loads an SIV project to read it; nothing there is run or saved.

## When a check fails
- If the cause is small and clearly in this package (e.g. a wrong lixoftConnectors field name),
  fix it, add or update a test, and list the fix in the results.
- Otherwise don't fix it; record it in the results and in `TODO.md`.

## Reporting back
1. Write `dev/monolix_check_results.md`: for each section, pass/fail and anything recorded
   (default scenario, criterion values, error-parameter names, failures with messages).
2. Update `TODO.md`: remove items that are now confirmed, add anything new.
3. Commit (with Maddie's OK) and push to `package-skeleton`, so the work can continue on the
   other computer with `git pull`.
