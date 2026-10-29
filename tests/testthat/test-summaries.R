# Fake IGS: 6 runs. k_pop is well identified; v_pop follows its initial value and k_pop/v_pop
# aren't correlated otherwise.
fake_igs <- function() {
  runs <- 1:6
  iv_v <- c(1, 2, 4, 8, 16, 32)
  te <- rbind(
    data.frame(parameter = "k_pop", estimate = c(1.00, 1.01, 0.99, 1.02, 1.5, 3), run = runs),
    data.frame(parameter = "v_pop", estimate = iv_v * 1.1, run = runs),
    data.frame(parameter = "omega_k", estimate = 0.3, run = runs),
    data.frame(parameter = "lag_pop", estimate = 0.5, run = runs)
  )
  te$rse <- ifelse(te$parameter == "v_pop", 80, 10)
  te$bicc <- c(100, 100.2, 100.4, 100.6, 120, 150)[te$run]
  iv <- rbind(
    data.frame(parameter = "k_pop", initial_value = c(5, 0.1, 2, 0.5, 9, 0.2), run = runs, method = "MLE"),
    data.frame(parameter = "v_pop", initial_value = iv_v, run = runs, method = "MLE"),
    data.frame(parameter = "omega_k", initial_value = 1, run = runs, method = "MLE"),
    data.frame(parameter = "lag_pop", initial_value = 0.5, run = runs, method = "FIXED")
  )
  sh <- data.frame(run = runs, parameter = "k", conditionalMean = 40, conditionalMode = 45, conditionalDist = 40)
  list(tabestimates = te, initialvalues = iv, shrinkage = sh)
}

test_that("identifiability_summary flags weakly identified parameters", {
  s <- identifiability_summary(fake_igs())
  expect_setequal(s$parameter, c("k_pop", "v_pop", "omega_k"))  # lag_pop is fixed
  k <- s[s$parameter == "k_pop", ]
  v <- s[s$parameter == "v_pop", ]
  expect_equal(k$best, 1.00)
  expect_equal(k$n_top_runs, 3)          # the 4th run would push the cumulative weight past 0.95
  expect_lt(k$top_fold, 1.1)
  expect_equal(v$cor_initial, 1)
  expect_match(v$flags, "follows initial value")
  expect_match(v$flags, "RSE > 50%")
  expect_match(v$flags, "spread > 2-fold")
  expect_false(grepl("follows", k$flags))
  expect_equal(s$shrinkage[s$parameter == "omega_k"], 45)
  expect_match(s$flags[s$parameter == "omega_k"], "shrinkage")
})


test_that("identifiability_summary takes parameter correlations from the top runs", {
  # runs 5 and 6 (poor BICc) would make k_pop and v_pop look strongly correlated
  s <- identifiability_summary(fake_igs())
  expect_equal(s$max_abs_cor[s$parameter == "k_pop"], 0.5)
})

test_that("compare_models counts wins, ties within delta and pairwise differences", {
  igs <- list(
    m1 = fake_candidate(c(100, 100, 100, 100), "beta_k_sex_1", p = c(0.001, 0.02, 0.2, 0.0001), cn = c(10, 20, 30, 40)),
    m2 = fake_candidate(c(101, 99, 105, 103)),
    m3 = fake_candidate(c(110, 110, 110, 110), c("beta_k_sex_1", "beta_k_age"), p = c(0.5, 0.5, 0.5, 0.5))
  )
  cmp <- compare_models(igs, delta = 2, fdr_method = "none")
  s <- cmp$summary
  expect_equal(s$times_best, c(3, 1, 0))
  expect_equal(s$times_within_delta, c(4, 2, 0))
  expect_equal(s$n_treatment_effects, c(1, 0, 2))
  expect_equal(s$median_condition_number, c(25, NA, NA))
  pw <- cmp$pairwise
  expect_equal(pw$runs_better_by_delta[pw$model == "m1" & pw$other == "m2"], 2)   # runs 3 and 4
  expect_equal(pw$runs_better_by_delta[pw$model == "m2" & pw$other == "m1"], 0)
  expect_equal(pw$runs_better_by_delta[pw$model == "m1" & pw$other == "m3"], 4)
  ws <- cmp$wald_summary
  expect_equal(ws$share_p05[ws$model == "m1"], 0.75)
  expect_equal(ws$share_p01[ws$model == "m1"], 0.5)
  expect_equal(nrow(cmp$omegas), 3)
})

test_that("compare_models only uses runs that succeeded in every model", {
  igs <- list(a = fake_candidate(c(100, 100, 100, 100)),
              b = fake_candidate(c(101, 101, 101), runs = c(1, 2, 4)))
  expect_message(cmp <- compare_models(igs), "run\\(s\\) failed in at least one model and are left out: 3")
  expect_equal(cmp$paired_runs, c(1, 2, 4))
  expect_equal(cmp$summary$n_runs, c(3, 3))
  expect_error(compare_models(list(fake_candidate(100))), "named list")
})

test_that("fdr adjustment is applied over all models together", {
  igs <- list(a = fake_candidate(c(100, 100), "beta_k_x", p = c(0.01, 0.04), runs = 1:2),
              b = fake_candidate(c(100, 100), "beta_k_x", p = c(0.03, 0.02), runs = 1:2))
  w <- compare_models(igs, fdr_method = "fdr")$wald
  expect_equal(w$p_adjusted, stats::p.adjust(c(0.01, 0.04, 0.03, 0.02), "fdr"))
})

test_that("profile_summary finds the lowest OFV and whether the 95% interval is bounded", {
  pr <- list(lambda_pop = fake_profile("lambda_pop", c(10, 50, 100, 200, 400), c(120, 104, 100, 103, 110)),
             p_pop = fake_profile("p_pop", c(1, 2, 3), c(100, 100.5, 101)))
  s <- profile_summary(pr, n_values = c(lambda_pop = 6, p_pop = 3))
  expect_equal(s$profiled_parameter, c("lambda_pop", "p_pop"))
  expect_equal(s$profiled_value, c(100, 1))
  expect_equal(s$ci_lower, c(100, 1)); expect_equal(s$ci_upper, c(200, 3))
  expect_equal(s$lower_bounded, c(TRUE, FALSE)); expect_equal(s$upper_bounded, c(TRUE, FALSE))
  expect_equal(s$n_failed, c(1, 0))
})

test_that("profile_grid uses round numbers, at least n of them", {
  expect_equal(profile_grid(c(0.5, 20), 10, log10 = TRUE), c(0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20))
  expect_equal(profile_grid(c(0.5, 20), 6, log10 = TRUE), c(0.5, 1, 2, 5, 10, 20))
  expect_equal(profile_grid(c(0.3, 25), 8, log10 = TRUE), c(0.3, 0.5, 1, 2, 5, 10, 20, 25))  # ends kept
  expect_equal(profile_grid(c(0.5, 20), 10), c(0.5, seq(2, 20, 2)))
  expect_equal(profile_grid(c(1, 9), 3), c(1, 5, 9))
  expect_equal(profile_grid(c(90, 110), 10, log10 = TRUE), seq(90, 110, 2))  # too narrow for log steps
  for (r in list(c(0.0013, 0.27), c(20, 600), c(1e-4, 1e4))) for (n in c(2, 7, 20)) for (lg in c(TRUE, FALSE)) {
    g <- profile_grid(r, n, lg)
    expect_gte(length(g), n); expect_equal(range(g), r); expect_false(is.unsorted(g))
  }
})

test_that("wald_note names treatment effects with p >= 0.05 or no p-value", {
  w <- data.frame(parameter = c("beta_ta_vac_1", "beta_m_vac_1", "beta_NPs_tlr7_1"), `p-value(sa)` = c(1e-4, 0.051, NA),
                  check.names = FALSE)
  expect_equal(wald_note(w), "not significant at 5%: beta_m_vac_1 (p = 0.051), beta_NPs_tlr7_1 (no p-value)")
  expect_equal(wald_note(w[1, ]), "")
  expect_equal(wald_note(NULL), "")
  expect_equal(wald_note(w[0, ]), "")
  expect_equal(wald_note(w, alpha = 0.1), "not significant at 10%: beta_NPs_tlr7_1 (no p-value)")
  # one p-value column per standard-error method: se_method picks one
  two <- data.frame(parameter = "beta_ta_vac_1", `p-value(sa)` = 0.2, `p-value(lin)` = 0.01, check.names = FALSE)
  expect_equal(wald_note(two, se_method = "linearization"), "")
  expect_match(wald_note(two, se_method = "stochasticApproximation"), "p = 0.2")
  expect_equal(wald_p_column(data.frame(parameter = "x", p.value = 1)), "p.value")
  expect_true(is.na(wald_p_column(data.frame(parameter = "x", statistics = 1))))
})

test_that("compare_models returns estimates relative to each model's best, and condition numbers", {
  igs <- list(m1 = fake_candidate(c(100, 102, 105, 101), cn = c(10, 20, 30, 40)),
              m2 = fake_candidate(c(90, 91, 99, 95), cn = c(5, NA, 7, 8)))
  cmp <- compare_models(igs)
  e <- cmp$estimates
  expect_equal(e$delta_bicc[e$model == "m1" & e$parameter == "k_pop"], c(0, 2, 5, 1))
  expect_equal(e$delta_bicc[e$model == "m2" & e$parameter == "k_pop"], c(0, 1, 9, 5))
  expect_equal(unique(e$run[e$model_best & e$model == "m1"]), 1)
  expect_equal(nrow(cmp$condition_numbers), 8)
  expect_null(compare_models(list(a = fake_candidate(1:4), b = fake_candidate(4:1)))$condition_numbers)
})
