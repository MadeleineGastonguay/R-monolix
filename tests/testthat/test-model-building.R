# scm_results_synthetic.rds is real getModelBuildingResults() output: a small SCM on the
# synthetic data (covariates tlr7 and vac, parameters ta and NPs; see dev/monolix_checks.R §9)
scm_fixture <- function() readRDS(fixture("scm_results_synthetic.rds"))

test_that("scm_results_table has one row per tested model", {
  tab <- scm_results_table(scm_fixture())
  expect_equal(nrow(tab), 8)
  expect_equal(tab$direction, c("base", rep("forward", 7)))
  expect_equal(tab$n_covariates, c(0, 1, 1, 1, 1, 2, 2, 2))
  expect_equal(tab$entry[tab$best_model], 3)
  expect_equal(tab$entry[tab$best_at_step], c(1, 3, 7))
  expect_equal(tab$covariate_model[3], "ta: vac")
  expect_equal(tab$covariate_model[1], "(none)")
})

test_that("scm_best_model returns Monolix's final covariate model", {
  b <- scm_best_model(scm_fixture())
  expect_equal(unclass(b)[names(b)], list(ta = "vac"))
  expect_equal(attr(b, "entry"), 3)
})

test_that("scm_truncated_model stops when BICc improves by less than delta", {
  mb <- scm_fixture()
  # Real path: base 178.7 -> ta~vac 166.7 (12.0 better) -> + NPs~tlr7 167.4 (worse): stop at step 1
  t2 <- scm_truncated_model(mb, delta_bicc = 2)
  expect_equal(attr(t2, "entry"), 3)
  expect_equal(unclass(t2)[names(t2)], list(ta = "vac"))
  # A delta bigger than the first improvement keeps the base model
  t20 <- scm_truncated_model(mb, delta_bicc = 20)
  expect_equal(attr(t20, "entry"), 1)
  expect_length(t20, 0)

  # Make the second step 1.5 better than the first: SCM would keep it, delta 2 wouldn't
  mb[[7]]$BICc <- mb[[3]]$BICc - 1.5
  t2b <- scm_truncated_model(mb, delta_bicc = 2)
  expect_equal(attr(t2b, "entry"), 3)
  expect_equal(attr(t2b, "path")$kept, c(TRUE, TRUE, FALSE))
  # ... while delta 1 adds NPs~tlr7
  t1 <- scm_truncated_model(mb, delta_bicc = 1)
  expect_equal(attr(t1, "entry"), 7)
  expect_equal(unclass(t1)[sort(names(t1))], list(NPs = "tlr7", ta = "vac"))
})

test_that("backward-elimination entries are ignored", {
  mb <- scm_fixture()
  back <- mb[[3]]  # after the 2-covariate step, a model with 1 covariate again
  back$BICc <- 100
  mb[[9]] <- back
  tab <- scm_results_table(mb)
  expect_equal(tab$direction[9], "backward")
  expect_false(tab$best_at_step[9])
  expect_equal(attr(scm_truncated_model(mb), "entry"), 3)
})

test_that("relationships become Monolix's parameters/covariates/locked table", {
  r <- model_building_relationships(list(include = list(ta = "vac"), exclude = list(m = c("tlr7", "vac"))))
  expect_equal(r, data.frame(parameters = c("ta", "m", "m"), covariates = c("vac", "tlr7", "vac"),
                             locked = c(TRUE, FALSE, FALSE)))
  expect_equal(model_building_relationships(list(exclude = list(NPs = "tlr7")))$locked, FALSE)
  expect_null(model_building_relationships(NULL))
  expect_null(model_building_relationships(list()))
  expect_null(model_building_relationships(list(include = list())))
})

test_that("a pair locked into every model doesn't change how the search is read", {
  # the fixture with m-tlr7 locked in: every model has one more covariate
  mb <- lapply(scm_fixture(), function(e) { e$individualModels["m", "tlr7"] <- TRUE; e })
  tab <- scm_results_table(mb)
  expect_equal(tab$direction, c("base", rep("forward", 7)))
  expect_equal(tab$n_covariates, c(0, 1, 1, 1, 1, 2, 2, 2) + 1)
  expect_equal(tab$entry[tab$best_at_step], c(1, 3, 7))
  expect_equal(unclass(scm_best_model(mb))[c("ta", "m")], list(ta = "vac", m = "tlr7"))
  ref <- scm_truncated_model(scm_fixture())
  expect_equal(attr(scm_truncated_model(mb), "entry"), attr(ref, "entry"))
})

test_that("scm_delta2 covariates removed by Monolix's backward step are flagged", {
  mb <- fake_mb_backward()
  expect_equal(scm_results_table(mb)$direction, c("base", rep("forward", 6), "backward", "backward"))
  expect_equal(unclass(scm_truncated_model(mb))[c("ta", "m")], list(ta = "vac", m = "tlr7"))
  expect_equal(scm_delta2_removed_backward(mb), data.frame(parameter = "ta", covariate = "vac"))
  note <- scm_delta2_backward_note(mb)
  expect_match(note, "scm_delta2 has ta-vac, which Monolix's backward elimination removed")
  expect_match(note, "scm_delta2_reduced: \\{m: \\[tlr7\\]\\}$")

  # no backward removal of a scm_delta2 pair: no note
  expect_null(scm_delta2_backward_note(scm_fixture()))
  expect_equal(nrow(scm_delta2_removed_backward(scm_fixture())), 0)
})
