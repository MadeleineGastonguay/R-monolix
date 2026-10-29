# step3_model_building.method: SCM, COSSAC, both, or none; and the candidate models each adds
cfg_with <- function(method, candidates = "~", compare = "all", final = "~") {
  read_analysis_config(make_config(list(
    "  method: SCM" = paste0("  method: ", method),
    "  fixed: ~" = "  fixed: {}",
    "  candidate_models: ~" = paste0("  candidate_models: ", candidates),
    "  compare: ~" = paste0("  compare: ", compare),
    "  final_model: ~" = paste0("  final_model: ", final))))
}
problems <- function(cfg, step = 4) tryCatch({ validate_config(cfg, step); "" }, error = function(e) conditionMessage(e))
scm_fixture <- function() readRDS(fixture("scm_results_synthetic.rds"))

test_that("method can be SCM, COSSAC, both or none, in any case", {
  expect_equal(model_building_methods(cfg_with("SCM")), "SCM")
  expect_equal(model_building_methods(cfg_with("cossac")), "COSSAC")
  expect_equal(model_building_methods(cfg_with("[SCM, COSSAC]")), c("SCM", "COSSAC"))
  expect_equal(model_building_methods(cfg_with("none")), character())
  expect_error(model_building_methods(cfg_with("[SCM, none]")), "must be SCM, COSSAC")
  expect_error(model_building_methods(cfg_with("samba")), "must be SCM, COSSAC")
})

test_that("each method adds its own automatic candidates", {
  expect_equal(candidate_names(cfg_with("SCM", "{am1: {ta: [vac]}}")), c("scm_best", "scm_delta2", "am1"))
  expect_equal(candidate_names(cfg_with("[SCM, COSSAC]", "{}")), c("scm_best", "scm_delta2", "cossac_best"))
  expect_equal(candidate_names(cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}}")), c("am1", "am2"))
})

test_that("with method none, at least two hand-written models are needed", {
  ok <- cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}}", final = "am2")
  expect_equal(problems(ok), "")
  expect_match(problems(cfg_with("none", "{am1: {ta: [vac]}}", final = "am1")), "at least two models when method is none")
  expect_match(problems(cfg_with("none", "{}", final = "am1")), "at least two models when method is none")
  # until candidates are written, it's a checkpoint naming what's needed
  err <- tryCatch(validate_config(cfg_with("none"), 3), pipeline_checkpoint = function(e) conditionMessage(e))
  expect_match(err, "at least two named treatment-effect models")
  # automatic names don't exist without a search
  expect_match(problems(cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}}", "[am1, scm_best]", "am1")),
               "compare has unknown model\\(s\\): scm_best")
})

test_that("search settings are only checked when a search runs", {
  cfg <- cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}}", final = "am1")
  cfg$step3_model_building$criterion <- "AIC"
  cfg$step3_model_building$covariates <- "weight"
  expect_equal(problems(cfg), "")
  cfg$step3_model_building$method <- "SCM"
  expect_match(problems(cfg), "criterion must be BIC")
})

test_that("cossac_best can be compared and chosen; compare needs two models", {
  expect_equal(problems(cfg_with("[SCM, COSSAC]", "{}", "[scm_best, cossac_best]", "cossac_best")), "")
  expect_match(problems(cfg_with("COSSAC", "{}", "[scm_best, cossac_best]", "cossac_best")), "unknown model\\(s\\): scm_best")
  expect_match(problems(cfg_with("SCM", "{}", "[scm_best]", "scm_best")), "at least two models")
})

test_that("resolve_candidates combines searches and hand-written models", {
  cfg <- cfg_with("SCM", "{am1: {ta: [vac], NPs: [vac]}, am2: {NPs: [tlr7]}}")
  cm <- resolve_candidates(cfg, list(SCM = scm_fixture()))
  # scm_delta2 is the same model as scm_best on this fixture (ta ~ vac), so it's fitted once
  expect_message(cm <- resolve_candidates(cfg, list(SCM = scm_fixture())), "scm_delta2 = scm_best")
  expect_equal(names(cm), c("scm_best", "am1", "am2"))
  expect_equal(cm$scm_best, list(ta = "vac"))
  expect_equal(cm$am1, list(ta = "vac", NPs = "vac"))
  expect_equal(attr(cm, "duplicates"), data.frame(dropped = "scm_delta2", same_as = "scm_best"))
})

test_that("SCM and COSSAC choosing the same model are fitted once; order of covariates doesn't matter", {
  cossac <- scm_fixture()  # stands in for COSSAC output with the same shape
  cfg <- cfg_with("[SCM, COSSAC]", "{am1: {NPs: [vac, tlr7]}, am2: {NPs: [tlr7, vac]}}")
  expect_message(cm <- resolve_candidates(cfg, list(SCM = scm_fixture(), COSSAC = cossac)),
                 "scm_delta2 = scm_best; cossac_best = scm_best; am2 = am1")
  expect_equal(names(cm), c("scm_best", "am1"))
  expect_error(resolve_candidates(cfg, list(SCM = scm_fixture())), "No COSSAC results")
})

test_that("with method none, resolve_candidates uses only the hand-written models", {
  cfg <- cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}, am3: {ta: [vac]}}")
  expect_message(cm <- resolve_candidates(cfg), "am3 = am1")
  expect_equal(names(cm), c("am1", "am2"))
})

test_that("step 3 fits every candidate; compare (checkpoint 2b) picks the ones for step 4", {
  hand <- "{am1: {ta: [vac], NPs: [vac]}, am2: {NPs: [tlr7]}, am3: {m: [vac]}}"
  cm <- suppressMessages(resolve_candidates(cfg_with("SCM", hand, "[am1, am2]"), list(SCM = scm_fixture())))
  expect_equal(names(cm), c("scm_best", "am1", "am2", "am3"))  # compare doesn't filter step 3
  expect_equal(step4_models(cfg_with("SCM", hand, "all"), cm), c("scm_best", "am1", "am2", "am3"))
  expect_equal(step4_models(cfg_with("SCM", hand, "[am2, scm_best]"), cm), c("am2", "scm_best"))
  # scm_delta2 was dropped as a duplicate of scm_best: it stands for scm_best
  expect_message(m <- step4_models(cfg_with("SCM", hand, "[scm_delta2, am1]"), cm), "scm_delta2 = scm_best")
  expect_equal(m, c("scm_best", "am1"))
  expect_error(suppressMessages(step4_models(cfg_with("SCM", hand, "[scm_delta2, scm_best]"), cm)), "only one distinct model")
  expect_error(step4_models(cfg_with("SCM", hand, "~"), cm), "checkpoint 2b")

  # until compare is decided, step 4 is a checkpoint
  err <- tryCatch(validate_config(cfg_with("SCM", hand, "~"), 4, decisions_needed = FALSE),
                  pipeline_checkpoint = function(e) conditionMessage(e))
  expect_match(err, "compare \\(checkpoint 2b\\): look at step3_model_building/candidates/single_fits.csv")
  # final_model must be among the compared models
  expect_match(problems(cfg_with("SCM", hand, "[am1, am2]", "am3")), "final_model 'am3' is not one of the compared models")
  expect_equal(problems(cfg_with("SCM", hand, "all", "am3")), "")
})

test_that("scm_delta2 is a BICc stop rule, so it isn't added with criterion LRT", {
  lrt <- read_analysis_config(make_config(list(
    "  fixed: ~" = "  fixed: {}", "  criterion: BIC" = "  criterion: LRT",
    "  scm_stop_rule_delta_bicc: 2" = "  scm_stop_rule_delta_bicc: ~", "  candidate_models: ~" = "  candidate_models: {}")))
  expect_equal(candidate_names(lrt), "scm_best")
  cm <- resolve_candidates(lrt, list(SCM = scm_fixture()))
  expect_equal(names(cm), "scm_best")
  expect_equal(nrow(attr(cm, "duplicates")), 0)
  # the stop rule isn't checked when it isn't used
  expect_invisible(validate_config(lrt, 3, decisions_needed = FALSE))
  # scm_delta2 stays reserved, so switching back to BIC can't clash with a hand-written model
  bad <- read_analysis_config(make_config(list(
    "  fixed: ~" = "  fixed: {}", "  criterion: BIC" = "  criterion: LRT", "  candidate_models: ~" = "  candidate_models: {scm_delta2: {ta: [vac]}}")))
  expect_match(problems(bad, 3), "scm_delta2 are reserved")
})
