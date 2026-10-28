# step3_model_building.method: SCM, COSSAC, both, or none; and the candidate models each adds
cfg_with <- function(method, candidates = "~", compare = "~", final = "~") {
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
  cfg <- cfg_with("none", "{am1: {ta: [vac]}, am2: {NPs: [tlr7]}, am3: {ta: [vac]}}", "[am1, am2]")
  cm <- resolve_candidates(cfg)
  expect_equal(names(cm), c("am1", "am2"))
  expect_length(attr(cm, "duplicates")$dropped, 0)
})
