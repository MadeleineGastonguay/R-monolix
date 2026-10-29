# A fake Monolix for running the pipeline end to end in tests. Projects are tracked by file;
# each fit returns made-up estimates for the parameters of the loaded project (base model plus
# betas for its covariates). `env$calls` records what the pipeline asked for.
BASE_PARAMS <- c("lambda_pop", "betas_pop", "ta_pop", "NPs_pop", "p_pop", "m_pop",
                 "omega_lambda", "omega_betas", "omega_ta", "omega_NPs", "omega_p", "omega_m", "a")

local_fake_monolix <- function(env, .local_envir = parent.frame()) {
  env$projects <- list()     # file -> treatment effects
  env$current <- NULL        # treatment effects of the loaded project
  env$calls <- list()
  log_call <- function(what, x = NULL) env$calls[[length(env$calls) + 1]] <- list(what = what, x = x)
  params_now <- function() {
    te <- env$current
    betas <- unlist(lapply(names(te), function(p) paste0("beta_", p, "_", te[[p]], "_1")))
    c(BASE_PARAMS, betas)
  }
  fake_igs <- function(initialValue_ranges, log10 = NULL, nits = 100, set_variance = NULL, set_treatment = NULL,
                       set_error = NULL, record_shrinkage = FALSE, fixed = NULL, predefined_initialvalues = NULL,
                       record_pvals = FALSE, record_indiv_params = FALSE, record_condition_number = FALSE,
                       record_table = FALSE, seed = NULL) {
    log_call("igs", list(ranges = initialValue_ranges, nits = nits, fixed = fixed, set_variance = set_variance,
                         predefined = predefined_initialvalues, seed = seed, params = params_now()))
    pp <- params_now()
    runs <- if (is.null(predefined_initialvalues)) seq_len(nits) else sort(unique(predefined_initialvalues$run))
    nb <- sum(grepl("^beta_", pp))
    te <- do.call(rbind, lapply(runs, function(r) data.frame(
      parameter = pp, estimate = ifelse(pp %in% names(fixed), unlist(fixed)[pp], 1 + r / 100), rse = 10, run = r,
      ofv = 90, aic = 95, bic = 99, bicc = 100 + (r %% 5) - nb)))
    te$estimate[is.na(te$estimate)] <- 1
    iv <- do.call(rbind, lapply(runs, function(r) data.frame(
      parameter = pp, initial_value = 1 + (r %% 3), run = r, method = ifelse(pp %in% names(fixed), "FIXED", "MLE"))))
    out <- list(tabestimates = te, initialvalues = iv, tabiters = NULL, phase = NULL,
                condition_number = data.frame(run = runs, condition_number = 10))
    if (nb) out$pvals <- data.frame(parameter = rep(pp[grepl("^beta_", pp)], length(runs)),
                                    p.value = 0.01, run = rep(runs, each = nb))
    out
  }
  testthat::local_mocked_bindings(
    setup_project = function(spec) { log_call("setup", spec); env$current <- spec$treatment_effects; invisible(NULL) },
    set_base = function(initial_values, omega_under_1 = FALSE, fixed = NULL) log_call("set_base", initial_values),
    initial_guess_sampling = fake_igs,
    save_igs_diagnostics = function(igs, out_dir, ...) { dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); invisible(NULL) },
    profile_likelihood = function(parameter, range, igs, nits, log10 = FALSE, values = NULL, ...) {
      v <- values %||% profile_grid(range, nits, log10)
      log_call("profile", list(parameter = parameter, values = v, igs = igs, scenario = env$scenario))
      # one row per estimated parameter and grid value; OFV lowest in the middle of the range
      g <- expand.grid(parameter = c(setdiff(BASE_PARAMS, parameter)), profiled_value = v, stringsAsFactors = FALSE)
      data.frame(g[, "parameter", drop = FALSE], estimate = 1, rse = NA, profiled_parameter = parameter,
                 profiled_value = g$profiled_value, ofv = 100 + 10 * (log(g$profiled_value) - mean(log(range)))^2,
                 aic = 1, bic = 1, bicc = 1)
    },
    covariate_model_building = function(projectFile, method, covariates = "all", parameters = NULL, relationships = NULL, ...) {
      log_call("search", list(method = method, file = projectFile, covariates = covariates, parameters = parameters,
                              relationships = relationships))
      mb <- readRDS(testthat::test_path("fixtures", "scm_results_synthetic.rds"))
      # like Monolix: the best model saved as <project>/final_model.mlxtran
      final <- file.path(sub("\\.mlxtran$", "", projectFile), "final_model.mlxtran")
      dir.create(dirname(final), recursive = TRUE, showWarnings = FALSE)
      writeLines("fake project", final)
      env$projects[[normalizePath(final)]] <- list(te = scm_best_model(mb)[])
      mb
    },
    read_likelihood = function(methods) c(ofv = 1, aic = 2, bic = 3, bicc = 100),
    result_methods = function(linearization = NULL) list(se = "stochasticApproximation", ll = "importanceSampling"),
    generate_report = function(project_file, ...) { log_call("report", project_file); invisible(NULL) },
    .env = .local_envir
  )
  testthat::local_mocked_bindings(
    saveProject = function(projectFile = NULL) {
      log_call("save", normalizePath(projectFile, mustWork = FALSE))
      dir.create(dirname(projectFile), recursive = TRUE, showWarnings = FALSE)
      writeLines("fake project", projectFile)
      env$projects[[normalizePath(projectFile)]] <- list(te = env$current)
      invisible(TRUE)
    },
    loadProject = function(projectFile) {
      log_call("load", normalizePath(projectFile))
      p <- env$projects[[normalizePath(projectFile)]]
      if (is.null(p)) stop("no such fake project: ", projectFile)
      env$current <- p$te
      invisible(TRUE)
    },
    getPopulationParameterInformation = function() data.frame(name = params_now(), initialValue = 1, method = "MLE"),
    runScenario = function() invisible(TRUE),
    getScenario = function() env$scenario %||% list(tasks = c(populationParameterEstimation = TRUE, standardErrorEstimation = TRUE),
                                                    linearization = FALSE),
    setScenario = function(scenario) { env$scenario <- scenario; invisible(TRUE) },
    computeChartsData = function(...) invisible(TRUE),
    # condition number 3 (eigenvalues 1.5 and 0.5); Wald p-values from env$wald_p, else 0.01
    getCorrelationOfEstimates = function(...) list(stochasticApproximation = matrix(c(1, 0.5, 0.5, 1), 2)),
    getTests = function(...) {
      b <- grep("^beta_", params_now(), value = TRUE)
      p <- vapply(b, function(x) env$wald_p[[x]] %||% 0.01, numeric(1))
      list(wald = data.frame(parameter = b, statistics = 3, `p-value(sa)` = unname(p), check.names = FALSE))
    },
    .package = "lixoftConnectors", .env = .local_envir
  )
  invisible(env)
}

calls_of <- function(env, what) Filter(function(x) x$what == what, env$calls)

# Edit keys of an analysis.yaml written by make_config()
set_yaml <- function(path, replace) {
  txt <- readLines(path)
  for (old in names(replace)) {
    hit <- grepl(old, txt, fixed = TRUE)
    if (!any(hit)) stop("no line matches: ", old)
    txt[hit] <- replace[[old]]
  }
  writeLines(txt, path)
}
