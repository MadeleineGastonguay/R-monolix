# set_base() talks to Monolix only through two calls, so mock them.
fake_popparams <- function() {
  data.frame(name = c("lambda_pop", "p_pop", "omega_lambda", "beta_lambda_vac_1", "beta_p_vac_1", "a"),
             initialValue = 1, method = "MLE")
}

run_set_base <- function(initial_values, fixed = NULL, ...) {
  captured <- NULL
  local_mocked_bindings(
    getPopulationParameterInformation = function() fake_popparams(),
    setPopulationParameterInformation = function(x) captured <<- x,
    .package = "lixoftConnectors"
  )
  set_base(initial_values, fixed = fixed, ...)
  captured
}

test_that("set_base matches fixed values by name, not position", {
  skip_if_not(monolix_available())
  iv <- data.frame(parameter = fake_popparams()$name, estimate = c(10, 2, 0.5, 0.3, NA, 0.4))
  out <- run_set_base(iv, fixed = c(p_pop = 0.7, lambda_pop = 50))
  expect_equal(out$initialValue[out$name == "lambda_pop"], 50)
  expect_equal(out$initialValue[out$name == "p_pop"], 0.7)
  expect_equal(out$method[out$name %in% c("lambda_pop", "p_pop")], c("FIXED", "FIXED"))
})

test_that("set_base fills only missing treatment effects with 0", {
  skip_if_not(monolix_available())
  iv <- data.frame(parameter = fake_popparams()$name, estimate = c(10, 2, 0.5, 0.3, NA, 0.4))
  out <- run_set_base(iv)
  expect_equal(out$initialValue[out$name == "beta_lambda_vac_1"], 0.3)
  expect_equal(out$initialValue[out$name == "beta_p_vac_1"], 0)
})

test_that("set_base raises omega initial values below 1 unless told not to", {
  skip_if_not(monolix_available())
  iv <- data.frame(parameter = fake_popparams()$name, estimate = c(10, 2, 0.5, 0.3, 0.1, 0.4))
  expect_equal(run_set_base(iv)$initialValue[3], 1)
  expect_equal(run_set_base(iv, omega_under_1 = TRUE)$initialValue[3], 0.5)
})

test_that("is_error_param picks out error-model parameters", {
  expect_equal(is_error_param(c("a", "b", "a_y1", "k_pop", "omega_k", "beta_k_sex_1", "corr_k_v")),
               c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})
