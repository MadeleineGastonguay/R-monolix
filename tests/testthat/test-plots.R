test_that("save_igs_diagnostics writes every output from saved results, without Monolix", {
  igs <- fake_igs_full()
  out <- file.path(withr_tempdir(), "wide")
  res <- suppressWarnings(save_igs_diagnostics(igs, out))
  expect_true(all(res$ok), info = paste(res$output[!res$ok], res$message[!res$ok], collapse = "; "))
  files <- list.files(out)
  expect_true(all(c("estimates.csv", "identifiability_summary.csv", "individual_fits.pdf",
                    "SAEM_trajectories.pdf", "cor_iv_estimates.pdf") %in% files))
})

test_that("a failing plot is reported and the others are still written", {
  igs <- fake_igs_full()
  igs$tabiters <- igs$tabiters[, c("iteration", "run")]  # breaks the SAEM plot only
  res <- suppressWarnings(save_igs_diagnostics(igs, file.path(withr_tempdir(), "x")))
  expect_false(res$ok[res$output == "saem_trajectories"])
  expect_true(res$ok[res$output == "estimates.csv"])
})

test_that("parameter_info marks parameters fixed in every run", {
  te <- data.frame(parameter = rep(c("a", "k_pop", "lag_pop", "omega_k"), 2), run = rep(1:2, each = 4),
                   estimate = c(0.1, 1, 0.5, 0.3, 0.2, 1.1, 0.5, 0.4), rse = c(5, 10, NA, 20, 5, 10, NA, 20))
  info <- parameter_info(te)
  expect_equal(info$name, c("k_pop", "lag_pop", "omega_k", "a"))
  expect_equal(info$method, c("MLE", "FIXED", "MLE", "MLE"))
})

test_that("plot_parameter_ranges applies report transforms", {
  igs <- fake_igs_full()
  wide <- data.frame(min = c(0.1, 1), max = c(10, 100), row.names = c("k_pop", "v_pop"))
  p <- plot_parameter_ranges(wide, igs$tabestimates, narrow_estimates = igs$tabestimates,
                             report_transforms = list(v_pop = list(label = "V (L)", expr = "x * 10")))
  d <- p$data
  expect_setequal(unique(d$label), c("k_pop", "V (L)"))
  expect_equal(d$max[d$label == "V (L)" & d$range == "Wide: initial values"], 1000)
  expect_setequal(levels(d$range)[levels(d$range) %in% d$range],
                  c("Wide: initial values", "Wide: estimates", "Narrow: initial values", "Narrow: estimates",
                    "Narrow: best-supported runs"))
})

test_that("fit plots use the best run and find the observation column", {
  igs <- fake_igs_full()
  p <- plot_individual_fits(igs$indivFit, igs$tabestimates)
  expect_equal(unique(p$data$run), igs$tabestimates$run[which.min(igs$tabestimates$bicc)])
  expect_equal(observation_column(igs$indivFit), "y")
  data <- data.frame(id = 1:3, grp = c("A", "B", "B"))
  pp <- plot_population_fits(igs$indivFit, igs$tabestimates, data = data, stratify_by = "grp")
  expect_setequal(as.character(unique(pp$data$group)), c("grp=A", "grp=B"))
  expect_error(plot_population_fits(igs$indivFit, stratify_by = "grp"), "`data` is needed")
})

test_that("model comparison plots build from compare_models output", {
  mk <- function(b) { x <- fake_igs_full(4); x$tabestimates$bicc <- x$tabestimates$bicc + b; x }
  cmp <- compare_models(list(m1 = mk(0), m2 = mk(1)))
  expect_s3_class(plot_model_comparison(cmp), "ggplot")
  expect_null(plot_treatment_pvalues(cmp))
})
