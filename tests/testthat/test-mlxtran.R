test_that("parse_mlxtran reads parameters, outputs and rescalings of the VACTLR7 model", {
  p <- parse_mlxtran(fixture("2regime_model_delay.txt"))
  expect_equal(p$parameters, c("lambda", "betas", "ta", "NPs", "p", "m"))
  expect_equal(p$outputs, "tot")
  expect_equal(p$tables, c("T", "V", "E", "I", "P"))
  expect_equal(guess_report_transforms(p),
               list(NPs_pop = list(label = "NP_pop", expr = "x * 10^5"),
                    betas_pop = list(label = "beta_pop", expr = "x * 10^(-7)")))
})

test_that("parse_mlxtran skips regressors and reads an output without braces", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("[LONGITUDINAL]", "input = {ka, V, Cl, WT}", "WT = {use = regressor}",
               "EQUATION:", "k = Cl/V  ; not a rescaling", "Vs = 1000 * V", "Cc = amt/V",
               "OUTPUT:", "output = Cc"), f)
  p <- parse_mlxtran(f)
  expect_equal(p$parameters, c("ka", "V", "Cl"))
  expect_equal(p$regressors, "WT")
  expect_equal(p$outputs, "Cc")
  expect_equal(names(p$rescalings), "V")
  expect_equal(p$rescalings$V$expr, "x * 1000")
})

test_that("guess_header_types recognises the VACTLR7-style columns", {
  ht <- guess_header_types(read_monolix_data(fixture("synthetic_vl_data.txt")))
  expect_equal(setNames(ht$type, ht$column),
               c(id = "id", time = "time", vl = "observation", censor = "cens", tlr7 = "catcov",
                 vac = "catcov", study = "catcov", intTLR7VAC = "catcov", intTLR7STUD = "catcov"))
})

test_that("guess_header_types falls back on data patterns and flags regressors", {
  d <- data.frame(ID = rep(1:3, each = 4), TIME = rep(0:3, 3), conc_obs = c(1:4, 2:5, 3:6) + 0.5,
                  WT = rep(c(60, 70, 80), each = 4), dose_now = c(1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1) * 5 + rep(0:3, 3))
  ht <- setNames(guess_header_types(d)$type, names(d))
  expect_equal(ht[["ID"]], "id")
  expect_equal(ht[["TIME"]], "time")
  expect_equal(ht[["WT"]], "catcov")
  expect_true(sum(ht == "observation") == 1)
  expect_true("ignore" %in% ht)
})
