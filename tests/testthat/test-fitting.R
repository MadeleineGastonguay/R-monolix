# Fake Monolix: 2 parameters; runScenario() fails whenever k_pop starts at the `fail_at` value.
mock_monolix <- function(env, fail_at = NULL) {
  env$set_values <- list()
  env$current <- NULL
  local_mocked_bindings(
    getPopulationParameterInformation = function() {
      data.frame(name = c("k_pop", "omega_k"), initialValue = c(1, 1), method = "MLE")
    },
    setPopulationParameterInformation = function(x) {
      env$current <- x
      env$set_values[[length(env$set_values) + 1]] <- x$initialValue[x$name == "k_pop"]
    },
    runScenario = function() {
      if (!is.null(fail_at) && isTRUE(all.equal(env$current$initialValue[1], fail_at))) stop("SAEM failed")
    },
    getEstimatedPopulationParameters = function() c(k_pop = env$current$initialValue[1] + 0.1, omega_k = 0.5),
    getEstimatedStandardErrors = function() {
      list(stochasticApproximation = data.frame(parameter = c("k_pop", "omega_k"), rse = c(10, 20)))
    },
    getEstimatedLogLikelihood = function() {
      list(importanceSampling = c(OFV = 1, AIC = 2, BIC = 3, BICc = env$current$initialValue[1]))
    },
    getSAEMiterations = function() {
      list(estimates = data.frame(k_pop = c(1, 1.1), omega_k = c(1, 0.5)), iterationNumbers = c(1, 2))
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

test_that("random sampling keeps drawing until nits runs succeed", {
  skip_if_not(monolix_available())
  env <- new.env()
  mock_monolix(env)
  ranges <- data.frame(min = c(1, 0.1), max = c(10, 1), row.names = c("k_pop", "omega_k"))
  res <- initial_guess_sampling(ranges, log10 = "k_pop", nits = 3, seed = 1)
  expect_equal(sort(unique(res$tabestimates$run)), 1:3)
  ks <- res$initialvalues$initial_value[res$initialvalues$parameter == "k_pop"]
  expect_true(all(ks >= 1 & ks <= 10))
  res2 <- initial_guess_sampling(ranges, log10 = "k_pop", nits = 3, seed = 1)
  expect_equal(res2$initialvalues, res$initialvalues)
})
