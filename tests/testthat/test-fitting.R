# Fake Monolix: 2 parameters; runScenario() fails whenever k_pop starts at the `fail_at` value.
mock_monolix <- function(env, fail_at = NULL, linearization = FALSE) {
  env$set_values <- list()
  env$current <- NULL
  se_method <- if (linearization) "linearization" else "stochasticApproximation"
  ll_method <- if (linearization) "linearization" else "importanceSampling"
  local_mocked_bindings(
    getScenario = function() list(tasks = c(populationParameterEstimation = TRUE), linearization = linearization),
    getPopulationParameterInformation = function() {
      data.frame(name = c("k_pop", "omega_k", "a"), initialValue = c(1, 1, 1), method = "MLE")
    },
    setPopulationParameterInformation = function(x) {
      env$current <- x
      env$set_values[[length(env$set_values) + 1]] <- x$initialValue[x$name == "k_pop"]
    },
    runScenario = function() {
      if (!is.null(fail_at) && isTRUE(all.equal(env$current$initialValue[1], fail_at))) stop("SAEM failed")
    },
    getEstimatedPopulationParameters = function() c(k_pop = env$current$initialValue[1] + 0.1, omega_k = 0.5, a = 0.2),
    getEstimatedStandardErrors = function() {
      stats::setNames(list(data.frame(parameter = c("k_pop", "omega_k", "a"), rse = c(10, 20, 5))), se_method)
    },
    getEstimatedLogLikelihood = function() {
      stats::setNames(list(c(OFV = 1, AIC = 2, BIC = 3, BICc = env$current$initialValue[1])), ll_method)
    },
    getSAEMiterations = function() {
      list(estimates = data.frame(k_pop = c(1, 1.1), omega_k = c(1, 0.5), a = c(1, 0.2)), iterationNumbers = c(1, 2))
    },
    .package = "lixoftConnectors",
    .env = parent.frame()
  )
}

test_that("predefined initial values are each fitted once and failed runs are not replaced", {
  skip_if_not(monolix_available())
  env <- new.env()
  mock_monolix(env, fail_at = 3)
  predefined <- data.frame(parameter = "k_pop", run = c(5, 6, 7), initial_value = c(2, 3, 4))

  expect_warning(
    expect_warning(
      res <- initial_guess_sampling(NULL, nits = 4, predefined_initialvalues = predefined),
      "using only the 3 predefined sets"),
    "1 of 3 predefined runs failed \\(runs 6\\)")

  expect_equal(sort(unique(res$tabestimates$run)), c(5, 7))
  expect_equal(unique(res$initialvalues$run), c(5, 6, 7))
  expect_equal(res$initialvalues$initial_value[res$initialvalues$parameter == "k_pop"], c(2, 3, 4))
  # 3 predefined fits, then the refit of the best run (run 5, started at 2)
  expect_equal(unlist(env$set_values), c(2, 3, 4, 2))
  expect_true(all(c("rse", "bicc") %in% names(res$tabestimates)))
})

test_that("random sampling keeps drawing until nits runs succeed; error parameters keep their start", {
  skip_if_not(monolix_available())
  env <- new.env()
  mock_monolix(env)
  ranges <- data.frame(min = c(1, 0.1), max = c(10, 1), row.names = c("k_pop", "omega_k"))
  res <- initial_guess_sampling(ranges, log10 = "k_pop", nits = 3, seed = 1, set_error = c(a = 0.3))
  expect_equal(sort(unique(res$tabestimates$run)), 1:3)
  iv <- res$initialvalues
  ks <- iv$initial_value[iv$parameter == "k_pop"]
  expect_true(all(ks >= 1 & ks <= 10))
  expect_equal(iv$initial_value[iv$parameter == "a"], rep(0.3, 3))
  # started at 0.3 but still estimated (not FIXED)
  expect_equal(unique(iv$method[iv$parameter == "a"]), "MLE")
  res2 <- initial_guess_sampling(ranges, log10 = "k_pop", nits = 3, seed = 1, set_error = c(a = 0.3))
  expect_equal(res2$initialvalues, res$initialvalues)
})

test_that("with linearization, results are read from the linearization method", {
  skip_if_not(monolix_available())
  env <- new.env()
  mock_monolix(env, linearization = TRUE)
  ranges <- data.frame(min = c(1, 0.1), max = c(10, 1), row.names = c("k_pop", "omega_k"))
  res <- initial_guess_sampling(ranges, nits = 2, seed = 1, set_error = c(a = 1))
  expect_equal(res$tabestimates$rse[res$tabestimates$parameter == "a"], c(5, 5))
  expect_false(anyNA(res$tabestimates$bicc))
})

test_that("missing standard errors give NA RSEs instead of failing the run", {
  skip_if_not(monolix_available())
  local_mocked_bindings(getEstimatedStandardErrors = function() list(), .package = "lixoftConnectors")
  expect_equal(read_rses(result_methods(FALSE), c("k_pop", "a")), c(k_pop = NA_real_, a = NA_real_))
})

test_that("shrinkage is read for whichever individual estimates exist", {
  skip_if_not(monolix_available())
  local_mocked_bindings(
    getEtaShrinkage = function() list(
      conditionalMode = data.frame(parameters = c("k", "v"), shrinkage = c(10, 20))),
    .package = "lixoftConnectors")
  sh <- read_shrinkage(3)
  expect_equal(sh$conditionalMode, c(10, 20))
  expect_true(all(is.na(sh$conditionalMean)))
  expect_equal(sh$run, c(3, 3))
})

test_that("result_methods maps linearization to Monolix's method names", {
  expect_equal(result_methods(TRUE), list(se = "linearization", ll = "linearization"))
  expect_equal(result_methods(FALSE), list(se = "stochasticApproximation", ll = "importanceSampling"))
})

test_that("condition number is NA, not an error, when the correlation matrix has nan", {
  # Seen with real Monolix: m_pop's standard error was nan, so its correlation row was all nan
  ok <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(c("k_pop", "v_pop"), c("k_pop", "v_pop")))
  expect_equal(condition_number_of(ok), 1.5 / 0.5)
  bad <- ok; bad["v_pop", ] <- NaN; bad[, "v_pop"] <- NaN
  expect_output(cn <- condition_number_of(bad, run = 2), "run 2.*\\(v_pop\\)")
  expect_identical(cn, NA_real_)
  expect_identical(condition_number_of(NULL), NA_real_)
  expect_equal(condition_number_of(as.data.frame(ok)), 3)
})
