test_that("calculate_BICw gives normalized weights ordered best first", {
  w <- calculate_BICw(fake_tabestimates())
  expect_equal(sum(w$w), 1)
  expect_equal(w$run[1], 1)
  expect_true(all(diff(w$cumprob) >= 0))
  expect_equal(w$w[1] / w$w[2], exp(0.5 / 2))
})

test_that("get_narrow_ranges keeps runs within the cumulative weight and widens to [best/2, 2*best]", {
  r <- get_narrow_ranges(fake_tabestimates(), p_threshold = 0.95)
  k <- r[r$parameter == "k_pop", ]
  # runs 1 and 2 carry ~95% of the weight; widening wins over their 1.0-1.2 span
  expect_equal(k$best, 1.0)
  expect_equal(k$min, 0.5)
  expect_equal(k$max, 2.0)
})

test_that("outlier_detector flags values beyond 1.5 IQR", {
  x <- c(1, 2, 3, 4, 100)
  expect_equal(outlier_detector(x), c(FALSE, FALSE, FALSE, FALSE, TRUE))
  expect_equal(outlier_detector(c(-100, 1, 2, 3, 4), upper_only = TRUE), rep(FALSE, 5))
})

test_that("filter_results drops BICc outlier runs and says so in the caption", {
  te <- fake_tabestimates()
  te$bicc[te$run == 4] <- 1e4
  out <- filter_results(te, remove_estimate_outliers = FALSE, remove_bicc_outliers = TRUE)
  expect_false(4 %in% out$tabestimates$run)
  expect_match(out$caption, "outlier BICcs removed")
  none <- filter_results(te, FALSE, FALSE)
  expect_equal(nrow(none$tabestimates), nrow(te))
  expect_equal(none$caption, "")
})
