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
