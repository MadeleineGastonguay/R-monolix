# Run the monolixPipelineR workflow on one analysis.
#
# 1. Make analysis.yaml next to this script:
#      monolixPipelineR::draft_config("model.txt", "data.txt", name = "my_analysis")
#    and fill in every value marked REQUIRED.
# 2. Source this script. It stops at each checkpoint and says which decision to record in
#    analysis.yaml; record it and source the script again. Steps whose saved results match the
#    config are loaded, not refitted.
# 3. To refit a step (e.g. after changing Monolix or wanting fresh results), set its RERUN flag
#    to TRUE. This deletes that step's saved results.
#
# Needs x86_64 R with lixoftConnectors, and MONOLIX_PATH set (e.g. in ~/.Renviron).

library(monolixPipelineR)

# Setup ####
CONFIG <- "analysis.yaml"
init_monolix()
config <- read_analysis_config(CONFIG)

# Step 1: two-phase initial guess sampling ####
# Wide IGS, narrow ranges, narrow IGS, diagnostics, identifiability summary, profile likelihoods.
RERUN_STEP1 <- FALSE
step1 <- run_step1(config, rerun = RERUN_STEP1)
# CHECKPOINT 1: look at <out_folder>/step1_igs/identifiability_summary.csv, the narrow/
# diagnostics and any profile likelihoods, then record decisions.fixed ({} to fix nothing)
# and decisions.fixed_rationale.

# Step 2: IGS with the fixed parameters ####
RERUN_STEP2 <- FALSE
step2 <- run_step2(config, rerun = RERUN_STEP2)

# Step 3: covariate search and candidate models ####
RERUN_STEP3 <- FALSE
step3 <- run_step3(config, rerun = RERUN_STEP3)
# CHECKPOINT 2 (stops inside run_step3 after the search): look at
# <out_folder>/step3_model_building/search/, then record step3_model_building.candidate_models
# ({} to compare only the automatic candidates).

# Step 4: paired comparison of the candidates ####
RERUN_STEP4 <- FALSE
step4 <- run_step4(config, rerun = RERUN_STEP4)
step4$comparison$summary
# CHECKPOINT 3: look at <out_folder>/step4_model_comparison/, then record
# step4_model_comparison.final_model and final_rationale.

# Final model and report ####
RERUN_FINAL <- FALSE
final <- run_final_model(config, rerun = RERUN_FINAL)
