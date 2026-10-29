# monolixPipelineR

monolixPipelineR fits nonlinear mixed effects models in Monolix from R, through
lixoftConnectors. It's built to assess practical identifiability of model
parameters through initial guess sampling (ISG; also known as a multi-start
approach [1-3], and implement rigorous comparison of several candidate treatment
effect models. The steps in this workflow are:

1.  **Identifiability.** Fit the model without covariates from many starting
    values over wide ranges, then again over narrower ranges set by the best
    runs. See which parameters the data determine, which follow their starting
    values, and which are correlated with each other.
2.  **Fix what the data can't determine.** Fix those parameters or remove their
    random effects, and repeat the IGS on the reduced model.
3.  **Covariate model.** Search covariate effects with Monolix's SCM and/or
    COSSAC to define data-driven treatment effect models, and fit any alternative 
    models you define. 
4.  **Model comparison.** Fit the candidates from the same sampled starting
    values and compare them run by run (BICc, how often each fits best,
    condition numbers, treatment effect significance).

Then the chosen model is refitted from its best run, and a report is written.

The pipeline works with any long-format dataset and any Mlxtran structural
model. Every setting and decision is in one YAML file, `analysis.yaml`. The
pipeline stops at each point where a decision is yours, and every result it
saves records the settings it was made with.

## Requirements

- MonolixSuite, with its lixoftConnectors R package (tested with 2024R1).
  lixoftConnectors comes with Monolix: install it from the Monolix folder with
  `install.packages("<Monolix folder>/connectors/lixoftConnectors.tar.gz", repos = NULL, type = "source")`.

- R 4.1 or later. **On a Mac with Apple silicon, lixoftConnectors only works in
  the x86_64 (Intel) build of R**, not the arm64 build. Check with
  `R.version$arch`.

- The `MONOLIX_PATH` environment variable set to the Monolix install folder,
  e.g. in `~/.Renviron`:

  ```         
  MONOLIX_PATH=/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/
  ```

- pandoc for the HTML pages and reports. 

Install the package from a copy of this repository:

``` r
remotes::install_github("MadeleineGastonguay/monolixPipelineR")
```

## Quick start

1.  Put your structural model (an Mlxtran file) and your data (long format, one
    row per observation or dose) in a folder for the analysis.

2.  Copy the driver script next to it and run it one section at a time:

    ``` r
    file.copy(system.file("scripts", "run_pipeline.R", package = "monolixPipelineR"), ".")
    ```

    Each section runs one step, prints its main result and lists where its files
    are. When a decision is needed, the step stops and says what to record in
    `analysis.yaml`. Record it and run the next section.
    
    If `anlaysis.yaml` does not already exist, the first step writes `analysis.yaml`, 
    with every key filled in. The column types are guessed from the data and marked 
    `# guessed`, so check them. Values only you can give are marked `# REQUIRED`, 
    mainly a plausible range for each parameter. The examples in its comments use 
    your model's parameter and covariate names.

`validate_config(read_analysis_config("analysis.yaml"), 1)` checks the file
without running anything. It reports every problem at once: typos in key names,
unknown parameters or columns, impossible values.

## The steps and their checkpoints

| Step | What it does | Where to look | Then decide (in `analysis.yaml`) |
|--------------------|--------------------|--------------------|--------------------|
| 1\. IGS, wide then narrow | Fits the base model from `step1_igs.wide_nits` random starts, then from narrower ranges | `step1_igs/igs_figures.html`, `identifiability_summary.csv` | **Checkpoint 1:** `decisions.fixed`, `decisions.fixed_rationale`, `decisions.random_effects_removed` |
| Profile likelihoods (optional) | Fits the model with one parameter fixed at each value on a grid | `step1_igs/profile_likelihood/` | (informs checkpoint 1) |
| 2\. IGS with the decisions | Refits with those parameters fixed and random effects removed; skipped if there are none | `step2_igs/igs_figures.html` | **Checkpoint 1b:** `decisions.fixed_confirmed` |
| 3a. Covariate search | SCM and/or COSSAC from step 2's best run | `step3_model_building/search/` | **Checkpoint 2:** `step3_model_building.candidate_models` |
| 3b. Candidate fits | One fit per candidate | `step3_model_building/candidates/single_fits.csv` | **Checkpoint 2b:** `step4_model_comparison.compare` |
| 4\. Paired comparison | All compared models from the same starting values | `step4_model_comparison/` | **Checkpoint 3:** `step4_model_comparison.final_model`, `step4_model_comparison.final_rationale` |
| Final model | Refits the chosen model from its best run; writes the report | `final_model/final_model/report.html` |  |

In the YAML, `~` means "not decided yet" and `{}` means "none". A decision left
at `~` stops the pipeline at its checkpoint, with a message saying what to
record.

### What the steps report

- **Identifiability summary** (steps 1 and 2): one row per parameter, with its
  spread among the best runs, its correlation with its starting value and with
  other parameters, its RSE and its shrinkage. Flags point out parameters worth
  a closer look. The best runs are the lowest-BICc runs that together make up
  95% of the BICc weight (`p_threshold`).
- **Figures** (steps 1 and 2): one HTML page per step with every figure. The
  plots about the consistency of the estimates (convergence, normalized
  estimates, CV, paired estimates, correlations) show only the best runs. The
  others (BICc against each estimate, estimates against starting values, SAEM
  trajectories) show every run, coloured by BICc. Each caption says which runs
  it shows.
- **Single fits** (step 3b): BICc, the condition number of the correlation
  matrix of the estimates, and a note naming covariate effects that aren't
  significant (Wald test, p ≥ 0.05).
- **Comparison** (step 4): how often each model fits best, BICc differences per
  pair, the distribution of BICc and condition numbers, how often each treatment
  effect is significant, and estimates against BICc for every model.

The search adds its own candidates: `scm_best` (SCM's best model), `scm_delta2`
(SCM stopped once BICc improves by less than `scm_stop_rule_delta_bicc`, with
criterion BIC only) and `cossac_best`. Candidates with the same covariate model
are fitted once.

## Example analysis.yaml

A model with three parameters (`param1`, `param2`, `param3`), one observation
(`y`) and three covariates, with every decision made. Keys you leave out take
their defaults (`config_defaults()`), including Monolix's SAEM settings
(`estimation`) and the tasks run for each fit (`scenario`).

``` yaml
analysis:
  name: example
  out_folder: results/example   # relative to this file
  r_seed: 1                     # seed for sampling starting values

data:
  file: data.csv
  header_types: {id: id, time: time, y: observation, cov1: catcov, cov2: catcov, cov3: contcov}

model:
  structural_model: model.txt
  error_model: {y: combined1}   # constant, proportional, combined1 or combined2
  random_effects_excluded: []   # parameters with no random effect at all
  distributions: {}             # e.g. {param3: normal}; default lognormal

step1_igs:
  initial_ranges:               # plausible range of each parameter, sampled uniformly
    param1: {min: 0.1, max: 100, log10: true}
    param2: {min: 0.01, max: 10, log10: true}
    param3: {min: 0.1, max: 10, log10: true}
  omega_range: [1, 5]           # starting random-effect SDs in the wide IGS
  wide_nits: 100                # runs in the wide IGS
  narrow: {nits: 100, p_threshold: 0.95, min_fold_width: ~, omega_initial: 1}
  profile_likelihood:           # optional, after step 1
    param3_pop: {range: [0.1, 10], n: 10, log10: true}

decisions:                      # filled in at the checkpoints
  fixed: {param3_pop: 1}        # CHECKPOINT 1
  fixed_rationale: "param3_pop spread over 10-fold among the best runs and its profile is flat"
  random_effects_removed: [param3]
  fixed_confirmed: true         # CHECKPOINT 1b

step2_igs: {nits: 100}

step3_model_building:
  method: SCM                   # SCM, COSSAC, [SCM, COSSAC] or none
  criterion: BIC                # BICc in Monolix; or LRT
  covariates: [cov1, cov2, cov3]
  parameters: [param1, param2]
  relationships: {include: {param1: [cov1]}}   # in every model; `exclude` for none
  candidate_models:             # CHECKPOINT 2: models of your own, besides scm_best
    model1: {param1: [cov1], param2: [cov2]}

step4_model_comparison:
  nits: 50
  compare: all                  # CHECKPOINT 2b: all, or a list of candidates
  final_model: scm_best         # CHECKPOINT 3
  final_rationale: "lowest BICc in 46 of 50 paired runs"

report:
  stratify_by: [cov1]           # groups in the report's plots
  group_labels: {"0": control, "1": treated}
```

`draft_config()` writes every key, with a comment on what each one does.

## Results

Everything goes under `analysis.out_folder`:

```         
results/example/
  step1_igs/              wide.mlxtran, narrow.mlxtran, igs_figures.html,
                          identifiability_summary.csv, narrow_ranges.csv,
                          wide/ and narrow/ (figures as PDFs, estimates.csv),
                          profile_likelihood/
  step2_igs/              igs_figures.html, identifiability_summary.csv, diagnostics/
  step3_model_building/
    search/               <method>_tested_models.csv, <method>_model_building.pdf, <method>/
    candidates/           <candidate>.mlxtran, single_fits.csv
  step4_model_comparison/ summary.csv, pairwise.csv, model_comparison.pdf,
                          bicc_vs_estimates.pdf, condition_numbers.pdf
  final_model/            final_model.mlxtran, final_model/report.html, decision.yaml
  history.yaml            one entry per finished step
```

`pipeline_paths(config)` gives these folders in R.

## Rerunning steps

Each step saves a `run_log.yaml` with the settings it used, hashes of the data
and model files, the seeds and the software versions. Next time, a step whose
saved settings match the config is loaded instead of refitted. If you change a
setting it depends on, the step stops and lists the change, for example
`decisions.fixed.param3_pop: 1 -> 2`. To refit it, set that step's `RERUN` flag
in the driver script. This replaces the step's results, and later steps then
need rerunning too.

`history.yaml` keeps one entry for every version of every step: when it ran,
what changed from the version it replaced, the decisions it used and its main
result (best BICc, best model). `read_history(config)` reads it.

The R sampling seed is `analysis.r_seed`. The Monolix seed is fixed, so a rerun
with the same settings should give the same results.

## Working with saved results

The plotting and summary functions don't need Monolix, so saved results can be
explored on any machine:

``` r
config <- read_analysis_config("analysis.yaml")
step1 <- readRDS(file.path(pipeline_paths(config)$step1, "results.rds"))
identifiability_summary(step1$narrow)
render_igs_figures(config, 1)       # remake step 1's figures page
plot_estimates(step1$narrow$tabestimates)
read_history(config)
```

## References
[1] Duchesne et al., BMC Bioinform., 2021 

[2] Lavielle et al., JPKPD, 2016

[3] Cassidy et al., Bull. Math. Biol., 2026

## Attributions
Madeleine S. Gastonguay designed the workflow and wrote the original analysis code, which fit a nonlinear mixed effects model to viral rebound kinetics after immunotherapy. That code is the basis for the package’s core functions: multi-start fitting from sampled starting values, narrowing the sampling ranges, covariate model building, the IGS diagnostic plots and the HTML report for a fitted Monolix project. These functions already took the model, data and statistical model as arguments.

Claude Code (Claude Opus 5.5, Anthropic) was used to add a layer around these functions so a new analysis needs only a model, a dataset and a settings file. That layer is:

- the analysis.yaml settings and their checks;
- the step-by-step pipeline with its checkpoints;
- run logs and history;
- the figure html pages;
- the tests and documentation.

It also removed the parts that were specific to the original study and replaced them with generalized functions. Madeleine S. Gastonguay directed and reviewed this work and checked it with Monolix runs.
