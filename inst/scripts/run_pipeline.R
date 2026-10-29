# monolixPipelineR: an analysis, one step at a time
#
# Run this script one section at a time (in RStudio, put the cursor in a section and press
# Ctrl/Cmd + Alt + T, or run it line by line). Each section runs one step, prints its main
# result and lists where its files are. Between steps you record a decision in
# analysis.yaml (a CHECKPOINT); each section reads the file again, so edit it and carry on.
#
# A step whose saved results match the current settings is loaded, not refitted. If you
# changed a setting it depends on, the step stops and lists what changed: set its RERUN flag
# to TRUE to refit it (this replaces its results; history.yaml keeps a record of each version).
#
# Needs x86_64 R with lixoftConnectors, and MONOLIX_PATH set (e.g. in ~/.Renviron).

library(monolixPipelineR)

# Draft analysis.yaml (first time only) ####
# Writes analysis.yaml from the model and data, with column types guessed from the data and
# every setting filled in. Check the guessed values (marked `# guessed`) and fill in every
# value marked REQUIRED, mainly a plausible range for each parameter. An existing
# analysis.yaml is never overwritten.
CONFIG <- "analysis.yaml"
MODEL <- "model.txt"   # Mlxtran structural model
DATA <- "data.csv"     # long-format data
NAME <- "my_analysis"  # results go to results/<NAME>
if (!file.exists(CONFIG)) draft_config(MODEL, DATA, name = NAME, path = CONFIG)

# Setup ####
# Paths in analysis.yaml are relative to the file. validate_config() lists every problem at once.
init_monolix()
config <- read_analysis_config(CONFIG)
validate_config(config, 1)
paths <- pipeline_paths(config)   # where each step writes its results
paths$root


# Step 1: initial guess sampling (IGS), wide then narrow ####
# Fits the base model (no covariates) from many random starting values over the ranges in
# step1_igs.initial_ranges, then again over narrower ranges set by the best wide runs.
# Results in paths$step1:
#   igs_figures.html             every figure: wide IGS, the ranges, narrow IGS
#   identifiability_summary.csv  one row per parameter, with flags for a closer look
#   wide/, narrow/               the same figures as PDFs, and estimates.csv
config <- read_analysis_config(CONFIG)
RERUN_STEP1 <- FALSE
step1 <- run_step1(config, rerun = RERUN_STEP1)
step1$identifiability

# Profile likelihoods (optional) ####
# For parameters that look poorly determined: list them in step1_igs.profile_likelihood,
# e.g. {param1_pop: {range: [0.1, 10], n: 10, log10: true}}. Doesn't change step 1, and
# only profiles whose grid changed are refitted.
# Results in paths$step1/profile_likelihood/ (also added to step 1's igs_figures.html)
config <- read_analysis_config(CONFIG)
RERUN_PROFILES <- FALSE
profiles <- run_profiles(config, step1, rerun = RERUN_PROFILES)
attr(profiles, "summary")

# CHECKPOINT 1 ####
# In analysis.yaml, under decisions:
#   fixed: {param1_pop: 10}       parameters the data can't determine, at a fixed value; {} for none
#   fixed_rationale: "..."        why
#   random_effects_removed: [param2]   optional: random effects to drop from step 2 on


# Step 2: IGS with the checkpoint 1 decisions ####
# Refits the base model with those parameters fixed and random effects removed, over step 1's
# narrow ranges. Skipped when nothing is fixed or removed (step 1's narrow IGS is used).
# Results in paths$step2: igs_figures.html, identifiability_summary.csv, diagnostics/
config <- read_analysis_config(CONFIG)
RERUN_STEP2 <- FALSE
step2 <- run_step2(config, rerun = RERUN_STEP2)
if (!isTRUE(attr(step2, "skipped"))) identifiability_summary(step2, config$step1_igs$narrow$p_threshold)

# CHECKPOINT 1b (only when step 2 ran) ####
#   decisions.fixed_confirmed: true    if nothing else needs fixing or removing. Otherwise change
#   fixed or random_effects_removed, set RERUN_STEP2 <- TRUE and run step 2 again.


# Step 3a: covariate search ####
# Searches covariate effects with step3_model_building.method (SCM, COSSAC, both, or none),
# starting from step 2's best run. The first time, this stops after the search at
# checkpoint 2.
# Results in paths$step3_search:
#   <method>_tested_models.csv    every model tested, with its BICc
#   <method>_model_building.pdf   the search path
#   <method>/                     Monolix's model-building folders and best model
config <- read_analysis_config(CONFIG)
RERUN_SEARCH <- FALSE
step3 <- run_step3(config, rerun = RERUN_SEARCH)

# CHECKPOINT 2 ####
#   step3_model_building.candidate_models: models of your own to compare with the search's
#   best (scm_best, scm_delta2, cossac_best), e.g. {model1: {param1: [cov1]}}; {} for none

# Step 3b: fit each candidate once ####
# One fit per candidate model, from step 2's best run.
# Results in paths$step3_candidates:
#   single_fits.csv               BICc, condition number and Wald-test note per candidate
#   <candidate>.mlxtran           the fitted projects
config <- read_analysis_config(CONFIG)
RERUN_CANDIDATES <- FALSE   # e.g. after changing candidate_models
step3 <- run_step3(config, rerun_candidates = RERUN_CANDIDATES)
step3$single_fits

# CHECKPOINT 2b ####
#   step4_model_comparison.compare: all, or a list of candidates, e.g. [scm_best, model1].
#   Leave out candidates that clearly fit worse. With one model, step 4 is just its IGS.


# Step 4: paired comparison of the candidates ####
# Fits every compared model from the same sampled starting values (step4_model_comparison.nits
# runs), so the models are compared run by run.
# Results in paths$step4:
#   summary.csv, pairwise.csv     how often each model fits best, BICc differences
#   model_comparison.pdf, bicc_vs_estimates.pdf, condition_numbers.pdf
config <- read_analysis_config(CONFIG)
RERUN_STEP4 <- FALSE
step4 <- run_step4(config, rerun = RERUN_STEP4)
step4$comparison$summary

# CHECKPOINT 3 ####
#   step4_model_comparison.final_model: the chosen model, and final_rationale: why


# Final model and report ####
# Refits the chosen model from its best step 4 run and writes its report.
# Results in paths$final:
#   final_model/report.html       the report
#   decision.yaml                 the choice, the rationale and its best run
#   paired_SAEM_population_estimates.pdf, normalized_estimate_boxplots.pdf, rses_vs_bicc.pdf
config <- read_analysis_config(CONFIG)
RERUN_FINAL <- FALSE
final <- run_final_model(config, rerun = RERUN_FINAL)


# History ####
# One entry per finished step: what changed from the version it replaced, its decisions and
# its main result (also in paths$root/history.yaml).
read_history(config)
