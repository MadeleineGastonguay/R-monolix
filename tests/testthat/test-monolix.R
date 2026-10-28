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
