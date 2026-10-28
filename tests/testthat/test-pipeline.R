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

  # candidates: scm_delta2 is the same model as scm_best on the fixture, so 2 models are compared
  set_yaml(path, list("  candidate_models: ~" = "  candidate_models: {am1: {ta: [vac], NPs: [tlr7]}}"))
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "final_model")
  expect_length(calls_of(env, "search"), 1)  # the search isn't rerun after checkpoint 2
  expect_equal(names(out$step3$candidates), c("scm_best", "am1"))
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
  out <- quietly(run_pipeline(path))
  expect_match(attr(out, "checkpoint"), "candidate_models")
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
  path <- small_config(list("  method: SCM" = "  method: none", "  fixed: ~" = "  fixed: {}",
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

test_that("profile likelihoods run for the configured parameters", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config(list("  profile_likelihood: {}" = "  profile_likelihood: {lambda_pop: {range: [20, 1000], n: 3}}"))
  quietly(run_pipeline(path))
  expect_equal(calls_of(env, "profile")[[1]]$x, "lambda_pop")
  root <- read_analysis_config(path)$analysis$out_folder
  expect_true(file.exists(file.path(root, "step1_igs", "profile_likelihood", "lambda_pop.csv")))
})

test_that("a step folder with files but no run log isn't overwritten", {
  skip_if_not(monolix_available())
  env <- new.env(); local_fake_monolix(env)
  path <- small_config()
  dir <- file.path(read_analysis_config(path)$analysis$out_folder, "step1_igs")
  dir.create(dir, recursive = TRUE); writeLines("x", file.path(dir, "old.txt"))
  expect_error(quietly(run_pipeline(path)), "has files but no run log")
})
