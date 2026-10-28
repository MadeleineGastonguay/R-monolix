# Checks that need a working Monolix install. Run section by section from the repo root,
# in x86_64 R, and compare with "Expect:" in each section. See HANDOFF.md.
# Self-contained: uses only files in this repo. Writes to a temporary folder, and in section 9
# to tests/testthat/fixtures/.

# Optional: setwd() to the repo root, and log printed output to a file Claude can read
# (not committed): sink("dev/monolix_check_output.log", split = TRUE)
# 0. Setup ####
stopifnot(R.version$arch == "x86_64")   # lixoftConnectors needs x86_64 R
Sys.getenv("MONOLIX_PATH")              # Expect: the Monolix install path, e.g.
# "/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/"
# If empty: add MONOLIX_PATH=... to ~/.Renviron and restart R.

devtools::load_all()
init_monolix()   # Expect: no error (runs the theophylline demo; estimates aren't printed)
# Without MONOLIX_PATH, pass the path, e.g.
# init_monolix("/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/")

# sink() only captures printed output, so errors, warnings and messages are printed explicitly
show_conditions <- function(expr) withCallingHandlers(
  tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n")),
  warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
  message = function(m) { cat("MESSAGE:", conditionMessage(m)); invokeRestart("muffleMessage") })

# A config for the test model + synthetic data, with short SAEM settings so fits are quick
work <- file.path(tempdir(), "monolix_checks")
dir.create(work, showWarnings = FALSE)
file.copy(c("tests/testthat/fixtures/2regime_model_delay.txt",
            "tests/testthat/fixtures/synthetic_vl_data.txt"), work, overwrite = TRUE)
cfg_path <- file.path(work, "analysis.yaml")
draft_config(file.path(work, "2regime_model_delay.txt"), file.path(work, "synthetic_vl_data.txt"),
             name = "checks", path = cfg_path, overwrite = TRUE)
txt <- readLines(cfg_path)
ranges <- c(lambda = "{min: 10, max: 1000, log10: true}", betas = "{min: 0.001, max: 100, log10: true}",
            ta = "{min: 0.01, max: 10, log10: true}", NPs = "{min: 0.001, max: 100, log10: true}",
            p = "{min: 0.01, max: 10, log10: true}", m = "{min: 0.01, max: 100, log10: true}")
for (k in names(ranges)) {
  txt <- sub(paste0("^    ", k, ": \\{min: ~, max: ~, log10: true\\}.*$"), paste0("    ", k, ": ", ranges[[k]]), txt)
}
txt <- sub("^  saem: .*$", "  saem: {burn_in: 5, exploratory: 50, smoothing: 20, smoothing_alpha: 0.7, exploratory_autostop: true, smoothing_autostop: true}", txt)
txt <- sub("^  fixed: ~.*$", "  fixed: {}", txt)
writeLines(txt, cfg_path)
cfg <- validate_config(read_analysis_config(cfg_path), 2)

# 1. Unit tests with Monolix installed ####
devtools::test()
# Expect: 0 failures, 0 skips (tests mock Monolix, but skip when lixoftConnectors is missing)

# 2. Monolix's own default scenario for a new project ####
spec <- model_spec(cfg, use_fixed = FALSE)
lixoftConnectors::newProject(modelFile = spec$model_file,
                             data = list(dataFile = spec$data_file, headerTypes = spec$header_types))
lixoftConnectors::getScenario()
# Record: which tasks are TRUE by default and the linearization value. This tells us what the
# SIV runs used for the tasks its scripts didn't set (conditional mode, conditional
# distribution, plots). The package defaults to all TRUE except plots.

# 3. Project setup from the config ####
setup_project(spec)
lixoftConnectors::getScenario()
# Expect: tasks as in cfg$scenario (all TRUE except plots), linearization FALSE
str(lixoftConnectors::getPopulationParameterEstimationSettings())
# Expect: burn-in 5, exploratory 50, smoothing 20, smoothing alpha 0.7, both autostops TRUE.
# The names come back differently from the camelCase ones setup_project() sets
# (first run gave NULL when subsetting by camelCase), so this prints them all; note the names.
lixoftConnectors::getProjectSettings()
# Expect: seed = 123456 (if the field is named differently, note its name)
lixoftConnectors::getIndividualParameterModel()$distribution
lixoftConnectors::getPopulationParameterInformation()
# Expect: lambda_pop ... m_pop, omega_*, and the error parameter "a"

# 4. A small initial guess sampling run ####
lixoftConnectors::saveProject(file.path(work, "igs_check.mlxtran"))
pp <- lixoftConnectors::getPopulationParameterInformation()$name
ranges_tbl <- sampling_ranges(spec, pp)
starts <- error_starts(spec, pp)
starts   # Expect: a = 1
res <- initial_guess_sampling(ranges_tbl, log10 = attr(ranges_tbl, "log10"), nits = 3, seed = 1,
                              set_error = starts, record_shrinkage = TRUE, record_pvals = TRUE,
                              record_indiv_params = TRUE, record_condition_number = TRUE, record_table = TRUE)
names(res)
head(res$tabestimates)
# Expect: columns parameter, estimate, rse, run, ofv, aic, bic, bicc; rse not all NA; 3 runs
subset(res$initialvalues, parameter == "a")
# Expect: initial_value 1 in every run, method MLE
res$shrinkage            # Expect: one row per parameter with random effects per run, numbers not all NA
res$condition_number     # Expect: one finite number per run
head(res$indivFit)       # Expect: predictions.txt columns incl. indivPred_mode
head(res$table_predictions)

# 5. Same seed, same results (Monolix seed is fixed) ####
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
res2 <- initial_guess_sampling(ranges_tbl, log10 = attr(ranges_tbl, "log10"), nits = 3, seed = 1, set_error = starts)
all.equal(res$tabestimates[order(res$tabestimates$run, res$tabestimates$parameter), c("parameter", "run", "estimate", "bicc")],
          res2$tabestimates[order(res2$tabestimates$run, res2$tabestimates$parameter), c("parameter", "run", "estimate", "bicc")],
          check.attributes = FALSE)
# Expect: TRUE. If not, the seed isn't being applied (or Monolix isn't deterministic here).

# 6. Linearization ####
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
sc <- cfg$scenario; sc$linearization <- TRUE
apply_scenario(sc)
lixoftConnectors::getScenario()$linearization   # Expect: TRUE
res_lin <- initial_guess_sampling(ranges_tbl, log10 = attr(ranges_tbl, "log10"), nits = 2, seed = 1,
                                  set_error = starts, record_condition_number = TRUE)
head(res_lin$tabestimates)   # Expect: rse and bicc filled in (from the linearization method)
res_lin$condition_number

# 7. Monolix version in the run log ####
log <- yaml::read_yaml(write_run_log(cfg, 1, file.path(work, "log_check")))
log$versions
# Expect: monolix is a version string (e.g. "2024R1"), not missing

# 8. Model-building settings ####
lixoftConnectors::getModelBuildingSettings()
# Record: the value of $criterion, and which values ?getModelBuildingSettings / ?runModelBuilding
# say are allowed. The package passes "BIC" for BICc.

# 9. Model-building results as a test fixture ####
# scm_truncated_model() (next step) reads getModelBuildingResults(). Run a small SCM on the
# synthetic data (2 covariates x 2 parameters, short SAEM) and save its results as a fixture.
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
apply_scenario(cfg$scenario)
lixoftConnectors::saveProject(file.path(work, "scm_check.mlxtran"))
set <- lixoftConnectors::getModelBuildingSettings()
names(set)   # Record: the setting names (the package sets strategy, criterion, covariates, parameters, useLin)
set$strategy <- "SCM"
set$criterion <- "BIC"
set$useLin <- FALSE
set$covariates <- c("tlr7", "vac")
if ("parameters" %in% names(set)) set$parameters <- c("ta", "NPs")
lixoftConnectors::runModelBuilding(settings = set)
mb <- lixoftConnectors::getModelBuildingResults()
length(mb)
str(mb[1:2], max.level = 3)
# Expect: one entry per tested model, with the tested parameter/covariate, BICc and which model is best.
saveRDS(mb, "tests/testthat/fixtures/scm_results_synthetic.rds")
capture.output(str(mb, max.level = 4), file = "tests/testthat/fixtures/scm_results_synthetic_str.txt")
# Both files come from the synthetic data only.

# 10. Report with plots: false ####
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
lixoftConnectors::runScenario()
lixoftConnectors::saveProject(file.path(work, "igs_check.mlxtran"))
generate_report(file.path(work, "igs_check.mlxtran"))
# Expect: an HTML report in <work>/igs_check/report.html, with no groups (no config given).
# A section Monolix can't produce for this project shows "Not available for this project: ..."
# instead of failing; note which.

# 11. Error-parameter names with several observations (optional) ####
demos <- list.files(lixoftConnectors::getDemoPath(), pattern = "PKPD.*\\.mlxtran$", recursive = TRUE, full.names = TRUE)
demos
lixoftConnectors::loadProject(demos[1])
names(lixoftConnectors::getContinuousObservationModel()$errorModel)
lixoftConnectors::getContinuousObservationModel()$errorModel
grep("^(a|b|c)", lixoftConnectors::getPopulationParameterInformation()$name, value = TRUE)
# Record: how error parameters are named per observation (e.g. a1/a2 or a_y1)

# 12. Can a structural parameter be called `a`? (optional) ####
m <- readLines(file.path(work, "2regime_model_delay.txt"))
m <- sub("input = {lambda,betas,ta,NPs,p,m}", "input = {lambda,betas,ta,NPs,p,a}", m, fixed = TRUE)
m <- gsub("\\bm\\b", "a", m)
writeLines(m, file.path(work, "model_param_a.txt"))
print(show_conditions(lixoftConnectors::newProject(modelFile = file.path(work, "model_param_a.txt"),
                                                   data = list(dataFile = spec$data_file, headerTypes = spec$header_types))))
print(show_conditions(lixoftConnectors::getPopulationParameterInformation()$name))
# Record: whether newProject fails, or what the error parameter and a_pop are called

# 13. Generalized report: groups, labels, rescaling, covariate effect ####
# A fit with one covariate effect (ta ~ vac), reported by tlr7 x vac groups with NPs_pop rescaled.
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
lixoftConnectors::setCovariateModel(ta = c(vac = TRUE))
# Save under the new name before running: a run writes results to the loaded project's folder,
# so running first would put this model's results in igs_check/ ("old inconsistent project").
lixoftConnectors::saveProject(file.path(work, "report_check.mlxtran"))
lixoftConnectors::runScenario()
lixoftConnectors::saveProject(file.path(work, "report_check.mlxtran"))
names(read_monolix_data(spec$data_file))   # Record: id column name (predictions.txt joins on it)
rep_path <- show_conditions(generate_report(file.path(work, "report_check.mlxtran"), config = cfg,
                            stratify_by = c("tlr7", "vac"),
                            group_labels = list("0.0" = "control", "1.1" = "tlr7 + vac"),
                            report_transforms = list(NPs_pop = list(label = "NP_pop", expr = "x * 10^5"))))
rep_path
# Expect: report in <work>/report_check/report.html. Open it and check:
# - "Report groups: tlr7, vac"; prediction plots with 4 groups (control, tlr7=1, vac=0, ...)
# - Population parameters: NP_pop = NPs_pop x 10^5, with its limits
# - Treatment Effects: beta_ta_vac_1 with a fold change
# - Standard deviations of random effects: one row per omega, CV column filled (if empty,
#   note what getEstimatedPopulationParameters(coefficientsOfVariation = TRUE) calls the CVs)
# - Prediction plots drawn as smooth curves (fine grid), not only at observation times, and no
#   "No individual fits chart data" note
# - Stratified VPC and prediction distribution split by tlr7 and vac
list.files(file.path(work, "report_check", "ChartsData", "IndividualFits"))
head(read.csv(list.files(file.path(work, "report_check", "ChartsData", "IndividualFits"),
                         pattern = "_fits\\.txt$", full.names = TRUE)[1]))
# Record: file names and columns (read_indfits() expects <obs>_fits.txt with ID, time, popPred
# and indivPredMode/indivPredMean; if the names differ, note them)
# - "Model variables for each individual": subtitle "Simulated with Simulx ..." and smooth curves
#   (if it says tables.txt, a "Simulx simulation ... failed" WARNING above gives the reason)

# 14. Simulx simulation of the model variables ####
# generate_report() runs this; here it's run on its own to see each step's output.
sim_file <- show_conditions(simulate_model_variables(file.path(work, "report_check.mlxtran"), n_points = 50))
sim_file   # Expect: <work>/report_check/ModelVariables/simulated_model_variables.txt
sims <- read_simulated_variables(file.path(work, "report_check"))
str(sims)
table(sims$variable)   # Expect: T, V, E, I, P with 50 x (number of individuals) rows each
setequal(unique(sims$id), unique(read_monolix_data(spec$data_file)$id))
# Expect: TRUE (ids are the data's ids; if FALSE, note what getSimulationResults() calls them)
show_conditions(lixoftConnectors::getLixoftConnectorsState())
# Expect: software monolix (switched back after the simulation)
# If the simulation fails, record the error. Then, to see what Simulx made of the project:
# lixoftConnectors::initializeLixoftConnectors("simulx", force = TRUE)
# lixoftConnectors::importProject(file.path(work, "report_check.mlxtran"))
# str(lixoftConnectors::getGroups()); names(lixoftConnectors::getIndividualElements())
# lixoftConnectors::initializeLixoftConnectors("monolix", force = TRUE)
# Record: any "Not available for this project" notes and their messages.

# 15. Step 5 diagnostics on a real IGS result (TODO: "Check in the next Monolix round") ####
# Uses `res` from section 4.
diag_dir <- file.path(work, "igs_diagnostics")
diag <- show_conditions(save_igs_diagnostics(res, diag_dir))
print(diag)   # Expect: ok TRUE for every output; record any message
list.files(diag_dir, recursive = TRUE)
ids <- identifiability_summary(res)
print(ids[, intersect(c("parameter", "shrinkage", "flags"), names(ids))])
res$shrinkage[res$shrinkage$run == 1, ]
# Expect: shrinkage filled in on every omega_ row (the shrinkage table's names, e.g. "lambda",
# map onto omega_lambda); record any omega_ row left NA
observation_column(res$indivFit)   # Expect: "vl" (the one column that isn't id, time or a prediction)
names(res$indivFit)
# individual_fits.pdf in diag_dir: open it and check fits and observations look right, and how
# censored observations (censor = 1 in the data) are drawn

# 16. Wald tests and model comparison (compare_models()) ####
# Two models fitted from the same initial values: the base model (igs_check) and ta ~ vac
# (report_check, section 13). Each is saved under its own name first, so the runs don't
# overwrite igs_check or report_check results.
lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
lixoftConnectors::saveProject(file.path(work, "cmp_base.mlxtran"))
cmp_base <- initial_guess_sampling(ranges_tbl, log10 = attr(ranges_tbl, "log10"), nits = 3, seed = 1,
                                   set_error = starts, record_pvals = TRUE, record_condition_number = TRUE)
lixoftConnectors::loadProject(file.path(work, "report_check.mlxtran"))
lixoftConnectors::saveProject(file.path(work, "cmp_vac.mlxtran"))
cmp_vac <- initial_guess_sampling(ranges_tbl, log10 = attr(ranges_tbl, "log10"), nits = 3, seed = 1,
                                  set_error = starts, record_pvals = TRUE, record_condition_number = TRUE)
str(cmp_vac$pvals)
# Record: the column names of getTests()$wald (compare_models() looks for p.value, pvalue or
# p_value) and whether beta_ta_vac_1 has a p-value
cmp <- show_conditions(compare_models(list(base = cmp_base, vac = cmp_vac)))
cmp$summary   # Expect: 2 models, n_runs = number of runs that succeeded in both
cmp$wald      # Expect: rows for beta_ta_vac_1 (NA p-values if its standard error is nan)

# 17. Rejected structural parameter names ####
# a) create_project() stops on a rejected model (until now tested only with mocks)
m <- readLines(file.path(work, "2regime_model_delay.txt"))
writeLines(gsub("\\bm\\b", "a", m), file.path(work, "model_param_a.txt"))
spec_a <- spec; spec_a$model_file <- file.path(work, "model_param_a.txt")
show_conditions(setup_project(spec_a))
# Expect: ERROR: Monolix could not create the project from .../model_param_a.txt
# b) Does Monolix also reject `b` or `c` with a constant error model (error parameter `a` only)?
for (nm in c("b", "c")) {
  writeLines(gsub("\\bm\\b", nm, m), file.path(work, paste0("model_param_", nm, ".txt")))
  ok <- show_conditions(lixoftConnectors::newProject(modelFile = file.path(work, paste0("model_param_", nm, ".txt")),
                                                     data = list(dataFile = spec$data_file, headerTypes = spec$header_types)))
  cat(nm, "accepted:", isTRUE(ok), "\n")
}
# Record: accepted TRUE/FALSE for b and c (validate_config() rejects a, b and c to be safe)

# 18. Covariate search through the package: SCM and COSSAC (TODO: "Covariate search options") ####
# Each search gets its own project file, so their ModelBuilding folders don't collide.
# covariate_model_building() also renders a report of the best model (with the Simulx step).
# The best-model reports are grouped by tlr7 x vac, as in section 13
cfg_groups <- cfg
cfg_groups$report$stratify_by <- list("tlr7", "vac")
cfg_groups$report$group_labels <- list("0.0" = "control", "1.1" = "tlr7 + vac")
# Rerunning: covariate_model_building() won't overwrite an earlier search, so remove this
# section's own earlier projects first
unlink(file.path(work, c("search_scm", "search_cossac", "search_scm.mlxtran", "search_cossac.mlxtran")),
       recursive = TRUE)
searches <- list()
for (method in c("SCM", "COSSAC")) {
  lixoftConnectors::loadProject(file.path(work, "igs_check.mlxtran"))
  apply_scenario(cfg$scenario)
  proj <- file.path(work, paste0("search_", tolower(method), ".mlxtran"))
  show_conditions(covariate_model_building(proj, method = method, covariates = c("tlr7", "vac"),
                                           parameters = c("ta", "NPs"), scenario = cfg$scenario,
                                           config = cfg_groups))
  # Record: any WARNING or MESSAGE about useSambaBeforeCossac, or about the strategy name
  lixoftConnectors::loadProject(proj)
  searches[[method]] <- show_conditions(lixoftConnectors::getModelBuildingResults())
  cat(method, ": ", length(searches[[method]]), " models, best: ",
      sum(vapply(searches[[method]], function(x) isTRUE(x$bestModel), logical(1))), "\n", sep = "")
}
# Expect: for each method, models > 0 and exactly 1 best (resolve_candidates() assumes this)
str(searches$COSSAC[[1]], max.level = 2)   # Expect: same fields as SCM (BICc, LL, tested, individualModels, bestModel)
list.files(file.path(work, c("search_scm", "search_cossac")), pattern = "\\.(html|mlxtran)$", recursive = TRUE)
# Expect: final_model.mlxtran and final_model/<SCM or COSSAC>_final_model_report.html in each
# Open both reports: "Report groups: tlr7, vac", prediction plots coloured by group, stratified VPC
cfg_both <- cfg; cfg_both$step3_model_building$method <- list("SCM", "COSSAC")
cands <- show_conditions(resolve_candidates(cfg_both, searches))
str(cands)   # Expect: scm_best, scm_delta2, cossac_best (duplicates dropped, with a message)

# sink()   # if logging
