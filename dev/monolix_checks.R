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
# Expect: an HTML report in <work>/igs_check/report.html. Sections tied to the SIV study columns
# may fail; the template isn't generalized yet. Note what fails.

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
# sink() only captures printed output, so errors, warnings and messages are printed explicitly
show_conditions <- function(expr) withCallingHandlers(
  tryCatch(expr, error = function(e) cat("ERROR:", conditionMessage(e), "\n")),
  warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") },
  message = function(m) { cat("MESSAGE:", conditionMessage(m)); invokeRestart("muffleMessage") })
print(show_conditions(lixoftConnectors::newProject(modelFile = file.path(work, "model_param_a.txt"),
                                                   data = list(dataFile = spec$data_file, headerTypes = spec$header_types))))
print(show_conditions(lixoftConnectors::getPopulationParameterInformation()$name))
# Record: whether newProject fails, or what the error parameter and a_pop are called
# sink()   # if logging
