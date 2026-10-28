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

fake_candidate <- function(bicc, betas = character(), p = NULL, cn = NULL, runs = 1:4) {
  te <- rbind(
    data.frame(parameter = "k_pop", estimate = 1, run = runs),
    data.frame(parameter = "omega_k", estimate = 0.3 + bicc / 1000, run = runs),
    if (length(betas)) data.frame(parameter = rep(betas, each = length(runs)), estimate = 0.1, run = runs)
  )
  te$rse <- 10
  te$bicc <- bicc[match(te$run, runs)]
  out <- list(tabestimates = te)
  if (!is.null(p)) out$pvals <- data.frame(parameter = betas[1], p.value = p, run = runs)
  if (!is.null(cn)) out$condition_number <- data.frame(run = runs, condition_number = cn)
  out
}

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
