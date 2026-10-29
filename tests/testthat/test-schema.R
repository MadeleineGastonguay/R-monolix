test_that("standardize_columns renames old SIV columns", {
  old <- data.frame(parameter = "k_pop", run = 1, initialValue = 2, rses = 10, am = "am1")
  new <- standardize_columns(old)
  expect_named(new, c("parameter", "run", "initial_value", "rse", "model"))
  expect_identical(standardize_columns(new), new)
})

test_that("standardize_columns converts old profile likelihood tables", {
  old <- data.frame(parameter = c("k_pop", "v_pop"), estimate = c(1, 2), lambda_pop = 50, bicc = 100)
  new <- standardize_columns(old, profiled_parameter = "lambda_pop")
  expect_equal(names(new)[1:2], c("profiled_parameter", "profiled_value"))
  expect_equal(unique(new$profiled_value), 50)
  expect_error(standardize_columns(old, "p_pop"), "No column")
})

test_that("standardize_columns refuses tables with both old and new names", {
  expect_error(standardize_columns(data.frame(rses = 1, rse = 1)), "both old and new")
})

test_that("load_igs_results converts every table in an old RData file", {
  path <- tempfile(fileext = ".RData")
  tabestimates <- data.frame(parameter = "k_pop", estimate = 1, rses = 5, run = 1, bicc = 10)
  initialvalues <- data.frame(parameter = "k_pop", initialValue = 2, method = "MLE", run = 1)
  phase <- matrix(1:2, ncol = 2)
  save(tabestimates, initialvalues, phase, file = path)
  res <- load_igs_results(path)
  expect_setequal(names(res), c("tabestimates", "initialvalues", "phase"))
  expect_true("rse" %in% names(res$tabestimates))
  expect_true("initial_value" %in% names(res$initialvalues))
  expect_identical(res$phase, phase)
})
