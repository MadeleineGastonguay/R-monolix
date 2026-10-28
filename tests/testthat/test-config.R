test_that("a drafted config reads back with the guesses and asks for the ranges", {
  dir <- withr_tempdir()
  file.copy(fixture(c("synthetic_vl_data.txt", "2regime_model_delay.txt")), dir)
  path <- file.path(dir, "analysis.yaml")
  expect_message(draft_config(file.path(dir, "2regime_model_delay.txt"), file.path(dir, "synthetic_vl_data.txt"),
                              name = "study2", path = path), "REQUIRED")
  expect_error(draft_config(file.path(dir, "2regime_model_delay.txt"), file.path(dir, "synthetic_vl_data.txt"),
                            name = "study2", path = path), "already exists")
  cfg <- read_analysis_config(path)
  expect_equal(cfg$data$file, normalizePath(file.path(dir, "synthetic_vl_data.txt")))
  expect_equal(unlist(cfg$data$header_types)[["censor"]], "cens")
  expect_equal(cfg$model$error_model, list(vl = "constant"))
  expect_equal(cfg$model$report_transforms$betas_pop$expr, "x * 10^(-7)")
  expect_equal(cfg$estimation, config_defaults()$estimation)
  err <- tryCatch(validate_config(cfg), error = function(e) conditionMessage(e))
  expect_match(err, "6 problem")
  for (p in c("lambda", "betas", "ta", "NPs", "p", "m")) expect_match(err, paste0("initial_ranges.", p, " min and max"))
})

test_that("checkpoints stop with a pipeline_checkpoint condition naming the key", {
  cfg <- read_analysis_config(make_config())
  expect_invisible(validate_config(cfg, 1))
  e <- tryCatch(validate_config(cfg, 2), error = function(e) e)
  expect_s3_class(e, "pipeline_checkpoint")
  expect_match(conditionMessage(e), "decisions.fixed")

  cfg <- read_analysis_config(make_config(list("  fixed: ~" = "  fixed: {lambda_pop: 50, p_pop: 0.7}")))
  expect_invisible(validate_config(cfg, 2))
  expect_invisible(validate_config(cfg, 3, decisions_needed = FALSE))
  e <- tryCatch(validate_config(cfg, 3), error = function(e) e)
  expect_s3_class(e, "pipeline_checkpoint")
  expect_match(conditionMessage(e), "candidate_models")
})

test_that("fixed: {} means fix nothing, not undecided", {
  cfg <- read_analysis_config(make_config(list("  fixed: ~" = "  fixed: {}")))
  expect_invisible(validate_config(cfg, 2))
  expect_null(model_spec(cfg)$fixed)
})

test_that("a full config validates through step 4 and catches bad candidates", {
  repl <- list("  fixed: ~" = "  fixed: {lambda_pop: 50, p_pop: 0.7}",
               "  candidate_models: ~" = "  candidate_models: {am1: {ta: [vac], NPs: [vac, intTLR7VAC]}, am2: {ta: [vac]}}",
               "  final_model: ~" = "  final_model: am1")
  expect_invisible(validate_config(read_analysis_config(make_config(repl)), 4))

  repl[["  candidate_models: ~"]] <- "  candidate_models: {scm_best: {ta: [vac]}, am3: {ta: [weight]}, am4: {q: [vac]}}"
  repl[["  final_model: ~"]] <- "  final_model: am9"
  err <- tryCatch(validate_config(read_analysis_config(make_config(repl)), 4), error = function(e) conditionMessage(e))
  expect_match(err, "scm_best are reserved for the automatic candidates")
  expect_match(err, "am3' uses non-covariate column\\(s\\): weight")
  expect_match(err, "am4' has unknown parameter\\(s\\): q")
  expect_match(err, "final_model 'am9'")
})

test_that("typos, bad values and unknown parameters are all reported together", {
  repl <- list("  wide_nits: 100" = "  wide_nit: 100",
               "    vl: constant" = "    vl: exponential",
               "  distributions: {}" = "  distributions: {tau: normal}",
               "  fixed: ~" = "  fixed: {lambda: 50}")
  err <- tryCatch(validate_config(read_analysis_config(make_config(repl)), 2), error = function(e) conditionMessage(e))
  expect_match(err, "Unknown keys \\(typo\\?\\): step1_igs.wide_nit")
  expect_match(err, "exponential")
  expect_match(err, "distributions has unknown parameter\\(s\\): tau")
  expect_match(err, "decisions.fixed has unknown parameter\\(s\\): lambda \\(use the _pop names\\)")
})

test_that("numbers written as 1e-7 are read as numbers, and y/on/no stay text", {
  repl <- list("    betas: {min: 0.001" = "    betas: {min: 1e-7, max: 1e2, log10: true}")
  cfg <- read_analysis_config(make_config(repl))
  expect_equal(cfg$step1_igs$initial_ranges$betas$min, 1e-7)
  expect_equal(cfg$step1_igs$initial_ranges$betas$max, 100)

  dir <- withr_tempdir()
  d <- data.frame(id = rep(1:2, each = 3), time = rep(0:2, 2), y = c(1.2, 2.5, 3.1, 0.7, 1.9, 2.2),
                  on = rep(c(0, 1), each = 3), no = rep(c(5, 6), each = 3))
  readr::write_tsv(d, file.path(dir, "d.txt"))
  file.copy(fixture("2regime_model_delay.txt"), dir)
  suppressMessages(draft_config(file.path(dir, "2regime_model_delay.txt"), file.path(dir, "d.txt"),
                                name = "t", path = file.path(dir, "a.yaml")))
  cfg <- read_analysis_config(file.path(dir, "a.yaml"))
  expect_equal(unlist(cfg$data$header_types),
               c(id = "id", time = "time", y = "observation", on = "catcov", no = "catcov"))
  expect_equal(names(cfg$model$error_model), "y")
})

test_that("model_spec resolves one project's settings", {
  cfg <- read_analysis_config(make_config(list("  fixed: ~" = "  fixed: {p_pop: 0.7, lambda_pop: 50}",
                                               "  distributions: {}" = "  distributions: {ta: normal}")))
  expect_error(model_spec(read_analysis_config(make_config())), "checkpoint 1")
  s <- model_spec(cfg, treatment_effects = list(ta = "vac", NPs = c("vac", "intTLR7VAC")))
  expect_equal(s$header_types, c("id", "time", "observation", "cens", rep("catcov", 5)))
  expect_equal(s$fixed, c(p_pop = 0.7, lambda_pop = 50))
  expect_equal(s$distributions, c(ta = "normal"))
  expect_equal(s$treatment_effects$NPs, c("vac", "intTLR7VAC"))
  expect_equal(s$error_model, list(vl = "constant"))
  expect_equal(s$initial_ranges$parameter, paste0(c("lambda", "betas", "ta", "NPs", "p", "m"), "_pop"))
  expect_null(model_spec(cfg, use_fixed = FALSE)$fixed)

  r <- sampling_ranges(s, c("lambda_pop", "betas_pop", "ta_pop", "NPs_pop", "p_pop", "m_pop",
                            "beta_ta_vac_1", "omega_lambda", "omega_ta", "a"))
  expect_false(any(c("beta_ta_vac_1", "a") %in% rownames(r)))
  expect_equal(unlist(r["omega_ta", ]), c(min = 1, max = 5))
  expect_equal(r["lambda_pop", "max"], 1000)
  expect_equal(attr(r, "log10"), paste0(c("lambda", "betas", "ta", "NPs", "p", "m"), "_pop"))
  expect_error(sampling_ranges(s, c("lambda_pop", "k_pop")), "No sampling range for: k_pop")
})

test_that("error parameters start at Monolix's defaults unless overridden", {
  cfg <- read_analysis_config(make_config(list("  fixed: ~" = "  fixed: {}")))
  expect_equal(error_starts(model_spec(cfg), c("lambda_pop", "omega_lambda", "a", "b", "c", "beta_ta_vac_1")),
               c(a = 1, b = 0.3, c = 1))
  cfg <- read_analysis_config(make_config(list("  fixed: ~" = "  fixed: {}", "  error_initial: {}" = "  error_initial: {a: 0.5}")))
  expect_equal(error_starts(model_spec(cfg), c("lambda_pop", "a")), c(a = 0.5))
})

test_that("error parameters can't be given sampling ranges, and their starts are checked", {
  repl <- list("    m: {min: 0.01" = "    m: {min: 0.01, max: 100, log10: true}\n    a: {min: 0.1, max: 2, log10: false}",
               "  error_initial: {}" = "  error_initial: {b: 0.3, a: -1}")
  path <- make_config(repl)
  writeLines(gsub("\\n", "\n", readLines(path), fixed = TRUE), path)
  err <- tryCatch(validate_config(read_analysis_config(path)), error = function(e) conditionMessage(e))
  expect_match(err, "error-model parameter\\(s\\) a aren't sampled")
  expect_match(err, "error_initial values must be positive")
  expect_match(err, "not in the error model \\(a\\): b")
})

test_that("the scenario has defaults and the pipeline's required tasks are enforced", {
  cfg <- read_analysis_config(make_config())
  expect_equal(cfg$scenario, config_defaults()$scenario)
  expect_false(cfg$scenario$linearization)

  repl <- list("  log_likelihood: true" = "  log_likelihood: false",
               "  conditional_mode: true" = "  conditional_mode: false",
               "  plots: false" = "  plots: maybe",
               "  linearization: false" = "  linearisation: true")
  err <- tryCatch(validate_config(read_analysis_config(make_config(repl))), error = function(e) conditionMessage(e))
  expect_match(err, "the pipeline needs conditional_mode, log_likelihood")
  expect_match(err, "scenario.plots must be true or false")
  expect_match(err, "Unknown keys \\(typo\\?\\): scenario.linearisation")

  ok <- read_analysis_config(make_config(list("  linearization: false" = "  linearization: true",
                                              "  standard_errors: true" = "  standard_errors: false")))
  expect_invisible(validate_config(ok))
  expect_true(model_spec(ok, use_fixed = FALSE)$scenario$linearization)
})

test_that("error parameters are numbered by observation when there are several", {
  expect_equal(error_param_names(list(vl = "constant")), "a")
  expect_equal(error_param_names(list(y1 = "combined1", y2 = "constant")), c("a1", "b1", "a2"))
  expect_equal(error_param_names(c(y1 = "proportional", y2 = "proportional")), c("b1", "b2"))
})

test_that("a structural parameter named like an error parameter is caught before Monolix", {
  # Monolix rejects it: "Name 'a' is already used as error parameter in the statistical model"
  path <- make_config(list("    m: {min: 0.01" = "    a: {min: 0.01, max: 100, log10: true}"))
  model <- file.path(dirname(path), "2regime_model_delay.txt")
  m <- readLines(model)
  m <- sub("input = {lambda,betas,ta,NPs,p,m}", "input = {lambda,betas,ta,NPs,p,a}", m, fixed = TRUE)
  writeLines(gsub("\\bm\\b", "a", m), model)
  err <- tryCatch(validate_config(read_analysis_config(path)), error = function(e) conditionMessage(e))
  expect_match(err, "Structural parameter\\(s\\) a clash with Monolix's error-model parameter names")
  expect_equal(reserved_error_names(2), c("a1", "a2", "b1", "b2", "c1", "c2"))
})
