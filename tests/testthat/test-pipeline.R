# The pipeline end to end with a fake Monolix (helper-pipeline.R): checkpoints, cached
# results, stale-settings checks and what each step hands the next.
small_config <- function(extra = list()) {
  make_config(c(list("  wide_nits: 100" = "  wide_nits: 6",
                     "  narrow: {nits: 100," = "  narrow: {nits: 5, p_threshold: 0.95, min_fold_width: ~, omega_initial: 1}",
                     "step2_igs: {nits: 100}" = "step2_igs: {nits: 4}",
                     "  nits: 50" = "  nits: 3"), extra))
}
quietly <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("the pipeline stops at each checkpoint and continues once the decision is recorded", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config()

  # checkpoint 1
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "decisions.fixed")
  root <- read_analysis_config(path)$analysis$out_folder
  s1 <- file.path(root, "step1_igs")
  expect_true(all(file.exists(file.path(s1, c("results.rds", "run_log.yaml", "narrow_ranges.csv",
                                              "identifiability_summary.csv", "wide.mlxtran", "narrow.mlxtran")))))
  igs <- calls_of(env, "igs")
  expect_length(igs, 2)
  expect_equal(igs[[1]]$x$nits, 6); expect_equal(igs[[1]]$x$seed, 1)
  expect_equal(igs[[2]]$x$nits, 5); expect_equal(igs[[2]]$x$seed, 2)
  expect_equal(igs[[2]]$x$set_variance, 1)                       # narrow: omegas start at 1
  expect_true(all(grepl("_pop$", rownames(igs[[2]]$x$ranges))))  # narrow ranges for _pop only

  # nothing fixed: step 2 is skipped, step 3 runs SCM and stops at checkpoint 2; step 1 is reused
  set_yaml(path, list("  fixed: ~" = "  fixed: {}"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "candidate_models")
  expect_length(calls_of(env, "igs"), 2)
  expect_length(calls_of(env, "search"), 1)
  expect_false(dir.exists(file.path(root, "step2_igs")))
  expect_true(file.exists(file.path(root, "step3_model_building", "search", "scm_tested_models.csv")))

  # candidates: scm_delta2 is the same model as scm_best on the fixture, so 3 models are fitted
  # once in step 3, then checkpoint 2b shows their single fits
  set_yaml(path, list("  candidate_models: ~" = "  candidate_models: {am1: {ta: [vac], NPs: [tlr7]}, am2: {m: [vac]}}"))
  env$wald_p <- list(beta_NPs_tlr7_1 = 0.3)  # am1's NPs-tlr7 effect isn't significant
  msgs <- character()
  out <- withCallingHandlers(suppressWarnings(run_pipeline(path)),
                             message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_match(attr(out, "checkpoint"), "compare \\(checkpoint 2b\\)")
  expect_true(any(grepl("Step 3 single fits", msgs)))
  expect_length(calls_of(env, "search"), 1)  # the search isn't rerun after checkpoint 2
  expect_equal(names(out$step3$candidates), c("scm_best", "am1", "am2"))
  sf <- read.csv(file.path(root, "step3_model_building", "candidates", "single_fits.csv"))
  expect_equal(names(sf), c("model", "covariate_model", "bicc", "delta_bicc", "condition_number", "wald_note", "source"))
  expect_equal(sf$condition_number, rep(3, 3))
  expect_equal(sf$wald_note[sf$model == "am1"], "not significant at 5%: beta_NPs_tlr7_1 (p = 0.3)")
  expect_true(all(sf$wald_note[sf$model != "am1"] %in% c("", NA)))
  expect_setequal(sf$model, c("scm_best", "am1", "am2"))
  # scm_best is copied from the search's final model, not fitted again; the others are fitted
  expect_equal(sf$source[sf$model == "scm_best"], "copied from the search")
  expect_equal(sf$source[sf$model != "scm_best"], rep("fitted from step 2's best values", 2))
  loads <- vapply(calls_of(env, "load"), function(x) x$x, "")
  expect_true(any(endsWith(loads, file.path("search", "scm", "final_model.mlxtran"))))
  expect_true(file.exists(file.path(root, "step3_model_building", "candidates", "scm_best.mlxtran")))
  expect_length(calls_of(env, "igs"), 2)  # step 4 hasn't started

  # compare leaves am2 out: only scm_best and am1 are fitted in step 4
  set_yaml(path, list("  compare: ~" = "  compare: [scm_best, am1]"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "final_model")
  expect_setequal(as.character(out$step4$comparison$summary$model), c("am1", "scm_best"))
  step4_igs <- tail(calls_of(env, "igs"), 2)
  # paired: the same initial values for every parameter the models share
  pv <- lapply(step4_igs, function(x) x$x$predefined)
  shared <- intersect(pv[[1]]$parameter, pv[[2]]$parameter)
  expect_equal(pv[[1]][pv[[1]]$parameter %in% shared, ], pv[[2]][pv[[2]]$parameter %in% shared, ], ignore_attr = TRUE)
  expect_true("beta_NPs_tlr7_1" %in% pv[[2]]$parameter)
  expect_equal(unique(pv[[2]]$initial_value[grepl("^beta_", pv[[2]]$parameter)]), 0)  # treatment_initial
  expect_equal(unique(pv[[1]]$initial_value[grepl("^omega", pv[[1]]$parameter)]), 1)  # omega_initial
  expect_true(file.exists(file.path(root, "step4_model_comparison", "summary.csv")))

  # final model
  set_yaml(path, list("  final_model: ~" = "  final_model: am1"))
  out <- quietly(run_pipeline(path))
  expect_null(attr(out, "checkpoint"))
  expect_true(file.exists(file.path(root, "final_model", "final_model.mlxtran")))
  expect_length(calls_of(env, "report"), 1)
  expect_length(calls_of(env, "igs"), 4)  # nothing refitted along the way
  fplots <- file.path(root, "final_model", c("paired_SAEM_population_estimates.pdf", "normalized_estimate_boxplots.pdf",
                                                    "rses_vs_bicc.pdf"))
  expect_true(all(file.exists(fplots)))

  # running again keeps the final model but redraws its plots, from step 4's saved results
  file.remove(fplots)
  out <- quietly(run_pipeline(path))
  expect_true(all(file.exists(fplots)))
  expect_length(calls_of(env, "report"), 1)
  expect_length(calls_of(env, "igs"), 4)
})

test_that("changed settings stop instead of reusing old results; rerun refits", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config()
  quietly(run_pipeline(path))
  set_yaml(path, list("  wide_nits: 6" = "  wide_nits: 7"))
  expect_error(quietly(run_pipeline(path)), "step1_igs.wide_nits: 6 -> 7")
  quietly(run_pipeline(path, rerun = 1))
  expect_equal(tail(calls_of(env, "igs"), 2)[[1]]$x$nits, 7)
})

test_that("fixed parameters go through step 2 with the fixed values", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config(list("  fixed: ~" = "  fixed: {lambda_pop: 50}"))
  # checkpoint 1b: step 2 runs, then stops before step 3 for the confirmation
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "decisions.fixed_confirmed")
  expect_length(calls_of(env, "search"), 0)
  s2dir <- file.path(read_analysis_config(path)$analysis$out_folder, "step2_igs")
  summ <- read.csv(file.path(s2dir, "identifiability_summary.csv"))
  expect_false("lambda_pop" %in% summ$parameter)   # fixed, so not summarised
  expect_true("p_pop" %in% summ$parameter)
  set_yaml(path, list("  fixed_confirmed: ~" = "  fixed_confirmed: true"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "candidate_models")
  expect_length(calls_of(env, "igs"), 3)   # step 2 loaded, not refitted
  s2 <- calls_of(env, "igs")[[3]]$x
  expect_equal(s2$fixed, c(lambda_pop = 50))
  expect_equal(s2$seed, 3)
  expect_false(isTRUE(attr(out$step2, "skipped")))
  expect_true(file.exists(file.path(read_analysis_config(path)$analysis$out_folder, "step2_igs", "best_estimates.csv")))
  # SCM starts from step 2's best run with lambda fixed
  sb <- calls_of(env, "set_base")[[1]]$x
  expect_equal(sb$estimate[sb$parameter == "lambda_pop"], 50)
})

test_that("method none skips the search and compares the hand-written models", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config(list("  method: SCM" = "  method: none", "  fixed: ~" = "  fixed: {}", "  compare: ~" = "  compare: all",
                            "  candidate_models: ~" = "  candidate_models: {am1: {ta: [vac]}, am2: {NPs: [tlr7]}}"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "final_model")
  expect_length(calls_of(env, "search"), 0)
  expect_equal(names(out$step3$candidates), c("am1", "am2"))
  expect_setequal(out$step4$comparison$summary$model, c("am1", "am2"))
})

test_that("SCM and COSSAC each get their own project", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config(list("  method: SCM" = "  method: [SCM, COSSAC]", "  fixed: ~" = "  fixed: {}"))
  quietly(run_pipeline(path))
  s <- calls_of(env, "search")
  expect_equal(vapply(s, function(x) x$x$method, ""), c("SCM", "COSSAC"))
  expect_equal(basename(vapply(s, function(x) x$x$file, "")), c("scm.mlxtran", "cossac.mlxtran"))
})

test_that("the covariate search gets the listed parameters and covariates", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  quietly(run_pipeline(small_config(list("  fixed: ~" = "  fixed: {}"))))
  expect_null(calls_of(env, "search")[[1]]$x$parameters)  # all: Monolix's default
  expect_null(calls_of(env, "search")[[1]]$x$relationships)
  expect_equal(calls_of(env, "search")[[1]]$x$covariates, c("tlr7", "vac", "study", "intTLR7VAC", "intTLR7STUD"))

  # a subset of the covariate columns
  env <- new.env(); local_fake_monolix(env)
  quietly(run_pipeline(small_config(list("  fixed: ~" = "  fixed: {}", "  covariates: [" = "  covariates: [vac, study]"))))
  expect_equal(calls_of(env, "search")[[1]]$x$covariates, c("vac", "study"))

  env <- new.env(); local_fake_monolix(env)
  path <- small_config(list("  fixed: ~" = "  fixed: {}", "  parameters: all" = "  parameters: [ta, NPs]",
                            "  relationships: ~" = "  relationships: {include: {ta: [vac]}, exclude: {NPs: [tlr7]}}"))
  quietly(run_pipeline(path))
  s <- calls_of(env, "search")[[1]]$x
  expect_equal(s$parameters, c("ta", "NPs"))
  expect_equal(s$relationships, data.frame(parameters = c("ta", "NPs"), covariates = c("vac", "tlr7"), locked = c(TRUE, FALSE)))

  path <- small_config(list("  fixed: ~" = "  fixed: {}", "  parameters: all" = "  parameters: [ta, ta_pop]"))
  expect_error(quietly(run_pipeline(path)), "unknown parameter\\(s\\): ta_pop")
})

test_that("profile likelihoods run after step 1 from the best narrow run, without touching step 1", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config()
  quietly(run_pipeline(path))   # step 1 only, no profiles yet
  root <- read_analysis_config(path)$analysis$out_folder
  s1 <- file.path(root, "step1_igs")
  narrow_saves <- function() sum(vapply(calls_of(env, "save"), function(x) basename(x$x) == "narrow.mlxtran", TRUE))
  before <- narrow_saves()
  step1_log <- readLines(file.path(s1, "run_log.yaml"))

  # adding profiles after step 1 ran: step 1 is loaded, not refitted or stopped
  set_yaml(path, list("  profile_likelihood: {}" =
    "  profile_likelihood: {lambda_pop: {range: [20, 1000], n: 5, log10: true}, p_pop: {range: [1, 9], n: 3}}"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "decisions.fixed")
  expect_length(calls_of(env, "igs"), 2)
  expect_equal(narrow_saves(), before)
  expect_identical(readLines(file.path(s1, "run_log.yaml")), step1_log)

  pr <- calls_of(env, "profile")
  expect_equal(vapply(pr, function(x) x$x$parameter, ""), c("lambda_pop", "p_pop"))
  expect_equal(pr[[1]]$x$values, c(20, 50, 100, 200, 500, 1000))   # round numbers, log10
  expect_equal(pr[[2]]$x$values, c(1, 5, 9))                       # round step, linear
  expect_false(pr[[1]]$x$scenario$tasks[["standardErrorEstimation"]])   # no standard errors
  iv <- pr[[1]]$x$igs
  best <- best_run_estimates(out$step1$narrow$tabestimates)
  expect_equal(iv$initialValue[iv$name == "ta_pop"], best$estimate[best$parameter == "ta_pop"])
  expect_equal(unique(iv$initialValue[grepl("^omega", iv$name)]), 1)
  expect_equal(iv$initialValue[iv$name == "a"], 1)
  expect_true(all(iv$method == "MLE"))

  pdir <- file.path(s1, "profile_likelihood")
  expect_true(all(file.exists(file.path(pdir, c("lambda_pop.csv", "p_pop.csv", "profile_summary.csv",
                                                "profile_likelihood.pdf", "profile.mlxtran", "run_log.yaml")))))
  expect_equal(names(out$profiles), c("lambda_pop", "p_pop"))
  expect_equal(attr(out$profiles, "summary")$n_failed, c(0, 0))

  # unchanged grids are loaded; a changed grid refits only that parameter
  quietly(run_pipeline(path))
  expect_length(calls_of(env, "profile"), 2)
  set_yaml(path, list("  profile_likelihood: {lambda_pop" =
    "  profile_likelihood: {lambda_pop: {range: [20, 1000], n: 5, log10: true}, p_pop: {range: [1, 9], n: 4}}"))
  quietly(run_pipeline(path))
  pr <- calls_of(env, "profile")
  expect_length(pr, 3)
  expect_equal(pr[[3]]$x$parameter, "p_pop")
  expect_equal(pr[[3]]$x$values, c(1, 2.5, 5, 7.5, 9))
  expect_length(calls_of(env, "igs"), 2)
  # the same grid listed as values is not refitted
  set_yaml(path, list("  profile_likelihood: {lambda_pop" =
    "  profile_likelihood: {lambda_pop: {range: [20, 1000], n: 5, log10: true}, p_pop: {values: [1, 2.5, 5, 7.5, 9]}}"))
  quietly(run_pipeline(path))
  expect_length(calls_of(env, "profile"), 3)

  # rerun = "profiles" refits them all, still leaving step 1 alone
  quietly(run_pipeline(path, rerun = "profiles"))
  expect_length(calls_of(env, "profile"), 5)
  expect_length(calls_of(env, "igs"), 2)
})

test_that("a step folder with files but no run log isn't overwritten", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config()
  dir <- file.path(read_analysis_config(path)$analysis$out_folder, "step1_igs")
  dir.create(dir, recursive = TRUE); writeLines("x", file.path(dir, "old.txt"))
  expect_error(quietly(run_pipeline(path)), "has files but no run log")
})

test_that("checkpoint 2 warns when scm_delta2 keeps a covariate Monolix removed backward", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  local_mocked_bindings(covariate_model_building = function(projectFile, method, ...) fake_mb_backward())
  path <- small_config(list("  fixed: ~" = "  fixed: {}"))
  expect_warning(out <- suppressMessages(run_pipeline(path)), "scm_delta2 has ta-vac.*scm_delta2_reduced: \\{m: \\[tlr7\\]\\}")
  expect_match(attr(out, "checkpoint"), "candidate_models")

  # once candidate_models is decided, no more warning
  set_yaml(path, list("  candidate_models: ~" = "  candidate_models: {}"))
  expect_no_warning(suppressMessages(run_pipeline(path)))
})
