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

test_that("BICc-vs-estimate overlays every model; condition numbers get violins", {
  igs <- list(m1 = fake_candidate(c(100, 102, 105, 101), "beta_k_sex_1", cn = c(10, 20, 30, 40)),
              m2 = fake_candidate(c(90, 91, 99, 95), cn = c(5, NA, 7, 8)))
  igs$m1$tabestimates$estimate[igs$m1$tabestimates$parameter == "beta_k_sex_1"] <- c(-0.2, 0.1, 0.3, 0)
  cmp <- compare_models(igs)
  p <- plot_bicc_vs_estimate_models(cmp)
  expect_s3_class(p, "patchwork")
  # k_pop is the same in every run (fixed): left out; omega_k and the treatment effect stay
  titles <- vapply(p$patches$plots, function(g) g$labels$title, "")
  expect_equal(c(titles, p$labels$title %||% character()), c("omega_k", "beta_k_sex_1"))
  last <- p[[length(p)]]
  expect_equal(sort(unique(as.character(last$data$model))), "m1")      # only m1 has the effect
  expect_false(any(vapply(last$scales$scales, function(s) "log10" %in% s$trans$name %||% "", logical(1))))
  b <- ggplot2::ggplot_build(p[[1]])
  expect_equal(sort(unique(b$data[[1]]$y)), sort(unique(c(0, 2, 5, 1, 0, 1, 9, 5))))
  expect_equal(nrow(b$data[[2]]), 2)                                     # one best-run triangle per model
  capped <- plot_bicc_vs_estimate_models(cmp, max_delta_bicc = 4)
  expect_match(capped$patches$annotation$caption, "3 model run\\(s\\) more than 4 above their model's best left out")

  cnp <- plot_condition_numbers(cmp)
  expect_s3_class(cnp, "ggplot")
  expect_equal(nrow(cnp$data), 7)                                        # the NA run is left out
  expect_match(cnp$labels$caption, "1 run\\(s\\) without a condition number left out")
  expect_null(plot_condition_numbers(compare_models(list(a = fake_candidate(1:4), b = fake_candidate(4:1)))))

  dir <- withr_tempdir()
  files <- suppressMessages(save_comparison_plots(cmp, dir))
  expect_setequal(basename(files), c("model_comparison.pdf", "condition_numbers.pdf", "bicc_vs_estimates.pdf"))
})

test_that("plot_individual_fits shows every run and censored observations as open circles", {
  fit <- data.frame(id = rep(1:2, each = 3), time = rep(c(0, 4, 8), 2), y = c(1, 2, 3, 1, 2, 3),
                    popPred = 1, indivPred_mode = 1.5, censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  indivFit <- rbind(transform(fit, run = 1), transform(fit, run = 2, indivPred_mode = 2))
  te <- data.frame(run = c(1, 2), bicc = c(10, 5))
  p <- plot_individual_fits(indivFit, te)
  b <- ggplot2::ggplot_build(p)
  layers <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_equal(unname(layers), c("GeomLine", "GeomLine", "GeomPoint"))
  expect_equal(unique(b$data[[1]]$colour), "grey70")       # run 1 drawn as another run
  expect_equal(sort(unique(b$data[[3]]$shape)), c(1, 16))  # open circles for censored
  expect_match(p$labels$caption, "Best run: 2")
  # without the censored column: plain points; without other runs: no grey layer
  p2 <- plot_individual_fits(indivFit[indivFit$run == 2, names(indivFit) != "censored"], te)
  expect_equal(length(p2$layers), 2)
})

test_that("package plots use theme_report", {
  p <- ggplot_report(data.frame(x = 1), ggplot2::aes(x, x))
  expect_equal(p$theme$panel.grid.minor, ggplot2::element_blank())
})

test_that("plot_profile_likelihood draws one panel per parameter with the 95% line", {
  pr <- list(lambda_pop = fake_profile("lambda_pop", c(10, 100, 1000), c(110, 100, 108)),
             p_pop = fake_profile("p_pop", c(1, 2, 3), c(101, 100, 102)))
  g <- plot_profile_likelihood(pr, log10 = "lambda_pop", best = c(lambda_pop = 90))
  expect_s3_class(g, "patchwork")
  expect_length(g$patches$plots, 1)   # lambda_pop; p_pop is g itself
  b <- ggplot2::ggplot_build(g$patches$plots[[1]])$data
  expect_equal(b[[1]]$yintercept, stats::qchisq(0.95, 1))
  expect_equal(b[[2]]$y, c(10, 0, 8))   # OFV relative to the lowest
  expect_equal(b[[2]]$x, c(1, 2, 3))    # log10 axis
  expect_equal(b[[4]]$xintercept, log10(90))
  expect_equal(ggplot2::ggplot_build(g)$data[[2]]$x, c(1, 2, 3))   # p_pop: linear
  expect_equal(g$patches$plots[[1]]$theme$panel.grid.minor, ggplot2::element_blank())
})
