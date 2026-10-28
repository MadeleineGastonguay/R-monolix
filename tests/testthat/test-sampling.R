ranges_fixture <- function() {
  r <- data.frame(min = c(1, 0.01, 1), max = c(1000, 1, 5), row.names = c("k_pop", "v_pop", "omega_k"))
  attr(r, "log10") <- c("k_pop", "v_pop")
  r
}
params <- c("k_pop", "v_pop", "omega_k", "beta_k_sex_1", "a")

test_that("sample_initial_values follows the sampling rules", {
  iv <- sample_initial_values(ranges_fixture(), params, n = 200, fixed = list(v_pop = 0.5),
                              set_error = c(a = 1), set_treatment = 0, seed = 1)
  expect_equal(nrow(iv), 200 * length(params))
  expect_equal(sort(unique(iv$run)), 1:200)
  val <- function(p) iv$initial_value[iv$parameter == p]
  expect_true(all(val("v_pop") == 0.5))       # fixed
  expect_true(all(val("a") == 1))             # error start
  expect_true(all(val("beta_k_sex_1") == 0))  # treatment start
  expect_true(all(val("k_pop") >= 1 & val("k_pop") <= 1000))
  expect_true(all(val("omega_k") >= 1 & val("omega_k") <= 5))
  # log-uniform: about half the draws below the geometric midpoint (31.6), not ~3%
  expect_gt(mean(val("k_pop") < sqrt(1000)), 0.35)
})

test_that("set_variance replaces sampling of omegas", {
  iv <- sample_initial_values(ranges_fixture(), c("k_pop", "omega_k"), n = 3, set_variance = 1, seed = 2)
  expect_equal(iv$initial_value[iv$parameter == "omega_k"], rep(1, 3))
})

test_that("sample_initial_values is reproducible with a seed and errors without a range", {
  a <- sample_initial_values(ranges_fixture(), c("k_pop", "omega_k"), n = 5, seed = 3)
  b <- sample_initial_values(ranges_fixture(), c("k_pop", "omega_k"), n = 5, seed = 3)
  expect_identical(a, b)
  expect_error(sample_initial_values(ranges_fixture(), c("k_pop", "m_pop"), n = 1), "No sampling range for m_pop")
})

test_that("prepare_initial_values rounds _pop values and resets omegas and betas", {
  best <- data.frame(parameter = c("k_pop", "omega_k", "beta_k_sex_1", "a"),
                     estimate = c(123.4, 0.4, 0.8, 0.27))
  iv <- prepare_initial_values(best, sig_figs = 1)
  expect_equal(iv$initial_value, c(100, 1, 0, 0.27))
  expect_equal(prepare_initial_values(best, sig_figs = 2, omega_value = 0.5)$initial_value[1:2], c(120, 0.5))
})

test_that("treatment effects start at 0 without set_treatment or a sampling range", {
  ranges <- data.frame(min = c(1, 10), max = c(2, 20), row.names = c("k_pop", "v_pop"))
  iv <- sample_initial_values(ranges, c("k_pop", "v_pop", "beta_k_sex_1"), n = 3, seed = 1)
  expect_equal(iv$initial_value[iv$parameter == "beta_k_sex_1"], c(0, 0, 0))
  # the other parameters are drawn as before (a beta uses no random numbers)
  iv2 <- sample_initial_values(ranges, c("k_pop", "v_pop"), n = 3, seed = 1)
  expect_equal(iv$initial_value[iv$parameter != "beta_k_sex_1"], iv2$initial_value)
  # other parameters without a range still error
  expect_error(sample_initial_values(ranges, c("k_pop", "p_pop"), n = 1), "No sampling range for p_pop")
})
