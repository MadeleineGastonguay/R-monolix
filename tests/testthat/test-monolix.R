test_that("package loads and reports Monolix availability without initializing it", {
  expect_type(monolix_available(), "logical")
})

test_that("init_monolix asks for a path when none is set", {
  skip_if_not(monolix_available())
  withr_path <- Sys.getenv("MONOLIX_PATH")
  Sys.setenv(MONOLIX_PATH = "")
  on.exit(Sys.setenv(MONOLIX_PATH = withr_path))
  expect_error(init_monolix(), "MONOLIX_PATH")
})

test_that("init_monolix prints and returns the demo estimates", {
  skip_if_not(monolix_available())
  local_mocked_bindings(initializeLixoftConnectors = function(...) invisible(NULL), .package = "lixoftConnectors")
  local_mocked_bindings(check_installation = function() c(ka_pop = 1.5, V_pop = 0.5))
  expect_output(est <- init_monolix("/some/path"), "theophylline demo estimates")
  expect_equal(est, c(ka_pop = 1.5, V_pop = 0.5))
  expect_true(init_monolix("/some/path", check = FALSE))
})
