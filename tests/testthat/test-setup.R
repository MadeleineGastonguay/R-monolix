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

test_that("set_base keeps the project's value for parameters it has no value for", {
  skip_if_not(monolix_available())
  # like the fixed c of a proportional error model: in the project, not in the estimates
  iv <- data.frame(parameter = c("lambda_pop", "p_pop", "omega_lambda"), estimate = c(10, 2, 0.5))
  out <- run_set_base(iv)
  expect_equal(out$initialValue[out$name == "a"], 1)
  expect_false(anyNA(out$initialValue))
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

test_that("create_project stops when Monolix rejects the model", {
  # Seen with real Monolix: a structural parameter named `a` makes newProject() print
  # "[ERROR] Name 'a' is already used as error parameter" and return FALSE, with no R error
  skip_if_not(monolix_available())
  local_mocked_bindings(newProject = function(...) FALSE, .package = "lixoftConnectors")
  expect_error(create_project("model_param_a.txt", "data.txt", headerTypes = c("id", "time", "observation"),
                              error_model = "constant", exclude_random_effects = NULL),
               "could not create the project from model_param_a.txt")
})

test_that("parameter_bounds follows each distribution, lognormal by default", {
  b <- parameter_bounds(list(p2 = "normal", p3 = list(distribution = "logitnormal", limits = c(5, 50)),
                             p4 = "probitnormal", p5 = "logitnormal"), c("p1", "p2", "p3", "p4", "p5"))
  expect_equal(b$parameter, paste0("p", 1:5, "_pop"))
  expect_equal(b$distribution, c("lognormal", "normal", "logitnormal", "probitnormal", "logitnormal"))
  expect_equal(b$lower, c(0, -Inf, 5, 0, 0))
  expect_equal(b$upper, c(Inf, Inf, 50, 1, 1))
  expect_equal(outside_bounds(c(0, 0.5, 1), 0, 1), c(TRUE, FALSE, TRUE))
})

test_that("create_project passes Monolix's distribution names and sets logit limits", {
  skip_if_not(monolix_available())
  dist <- NULL; lim <- list(); pop_set <- NULL
  state <- list(limits = list(p = c(0, 1)))
  local_mocked_bindings(
    newProject = function(...) TRUE, setErrorModel = function(...) TRUE,
    getIndividualParameterModel = function() list(distribution = c(k = "logNormal", p = "logNormal"),
                                                  limits = state$limits),
    setIndividualParameterDistribution = function(x) dist <<- x,
    # like Monolix when the starting value is outside the new limits: the first call is ignored
    setIndividualLogitLimits = function(...) { lim[[length(lim) + 1]] <<- list(...)
      if (length(lim) > 1) state$limits <<- list(...) ; TRUE },
    getPopulationParameterInformation = function() data.frame(name = c("k_pop", "p_pop"), initialValue = 1, method = "MLE"),
    setPopulationParameterInformation = function(x) pop_set <<- x,
    .package = "lixoftConnectors")
  local_mocked_bindings(set_project_settings = function(...) NULL)
  suppressMessages(capture.output(
    create_project("m.txt", "d.txt", headerTypes = "id", error_model = list(y = "constant"),
                   exclude_random_effects = NULL, distributions = c(k = "normal", p = "logitnormal"),
                   limits = list(p = c(5, 50)))))
  expect_equal(dist, c(k = "normal", p = "logitNormal"))
  expect_length(lim, 2)
  expect_equal(lim[[2]], list(p = c(5, 50)))
  expect_equal(pop_set$initialValue[pop_set$name == "p_pop"], 27.5)  # moved inside the limits
  expect_equal(pop_set$initialValue[pop_set$name == "k_pop"], 1)
})

test_that("create_project stops when Monolix doesn't take the limits", {
  skip_if_not(monolix_available())
  local_mocked_bindings(
    newProject = function(...) TRUE, setErrorModel = function(...) TRUE,
    getIndividualParameterModel = function() list(distribution = c(p = "logNormal"), limits = list()),
    setIndividualParameterDistribution = function(x) TRUE,
    setIndividualLogitLimits = function(...) FALSE,
    getPopulationParameterInformation = function() data.frame(name = "p_pop", initialValue = 1, method = "MLE"),
    setPopulationParameterInformation = function(x) TRUE,
    .package = "lixoftConnectors")
  expect_error(suppressMessages(create_project("m.txt", "d.txt", headerTypes = "id", error_model = list(y = "constant"),
                                               exclude_random_effects = NULL, distributions = c(p = "logitnormal"),
                                               limits = list(p = c(5, 50)))),
               "did not set the limits of p")
})
