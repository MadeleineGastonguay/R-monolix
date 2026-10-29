# The IGS figures: which runs each one shows, their captions, and the one-page HTML version.

test_that("top_runs keeps the lowest-BICc runs up to the BICc weight, and always the best run", {
  te <- data.frame(run = 1:4, parameter = "k_pop", estimate = 1, bicc = c(100, 100, 100, 130))
  expect_equal(top_runs(te, 0.95), c(1, 2))   # weights 1/3 each: the third would pass 95%
  expect_equal(top_runs(te, 0.1), 1)   # the best run alone carries more than 10%
})

test_that("captions say which runs are shown, and wrap long lines", {
  te <- data.frame(run = 1:20, bicc = 1:20)
  expect_equal(runs_note(te), "All 20 runs.")
  expect_match(runs_note(te, 1:18, 0.95, sets_ranges = TRUE), "^Top 18 of 20 runs: .*95% of the BICc weight \\(these set the narrow ranges\\)\\.$")
  expect_false(grepl("narrow ranges", runs_note(te, 1:18, 0.95)))
  expect_null(join_caption(NULL, "", NA))
  expect_null(join_caption(NULL, NULL))
  long <- join_caption(strrep("word ", 40), "second")
  expect_true(all(nchar(strsplit(long, "\n")[[1]]) <= 100))
  expect_equal(tail(strsplit(long, "\n")[[1]], 1), "second")
})

test_that("consistency plots show only the top runs; the others show every run", {
  igs <- fake_igs_full(12)
  top <- top_runs(igs$tabestimates, 0.95)
  skip_if(length(top) < 3 || length(top) == 12, "fixture needs between 3 and 11 top runs")
  figs <- igs_figures(igs, 0.95, sets_ranges = TRUE)
  expect_true(all(c("bicc_vs_estimates", "convergence_estimates", "estimate_boxplots", "paired_estimates",
                    "estimate_correlations", "coefficient_of_variation", "saem_trajectories", "individual_fits") %in% names(figs)))
  box <- figs$estimate_boxplots$draw()
  expect_setequal(unique(box$data$run), top)
  expect_match(ggplot2::get_labs(box)$caption, paste0("Top ", length(top), " of 12 runs"))
  expect_match(ggplot2::get_labs(box)$caption, "these set the narrow\\s+ranges")
  bv <- figs$bicc_vs_estimates$draw()
  expect_equal(length(unique(bv$data$run)), 12)
  expect_match(ggplot2::get_labs(bv)$caption, "All 12 runs.*Dashed line")
  conv <- figs$convergence_estimates$draw()
  expect_s3_class(conv, "patchwork")
  expect_match(conv$patches$annotation$caption, "Solid line: estimate of the best run")
  saem <- figs$saem_trajectories$draw()
  expect_match(saem$patches$annotation$caption, "0 \\(dashed line\\) is the\\s+switch")
})

test_that("with fewer than 3 top runs the consistency plots show every run and say why", {
  igs <- fake_igs_full(6)
  igs$tabestimates$bicc <- igs$tabestimates$bicc + ifelse(igs$tabestimates$run == 1, -50, 0)  # one run dominates
  figs <- igs_figures(igs, 0.95)
  box <- figs$estimate_boxplots$draw()
  expect_equal(length(unique(box$data$run)), 6)
  expect_match(ggplot2::get_labs(box)$caption, "only 1 run\\(s\\) make up 95%")
})

test_that("normalized boxplots keep every error-model parameter, in Monolix's order", {
  te <- expand.grid(run = 1:4, parameter = c("b", "omega_k", "k_pop", "c1"), stringsAsFactors = FALSE)
  te$estimate <- 1 + te$run / 10
  te$bicc <- 100 + te$run
  p <- estimate_boxplots(te)
  expect_false(anyNA(p$data$parameter))
  expect_equal(levels(p$data$parameter), c("k_pop", "omega_k", "b", "c1"))
})

test_that("the figures page has every figure of the wide and narrow IGS", {
  skip_if_not(rmarkdown::pandoc_available())
  igs <- fake_igs_full(8)
  root <- withr_tempdir()
  dir.create(file.path(root, "step1_igs"))
  saveRDS(list(wide = igs, narrow = igs, narrow_ranges = get_narrow_ranges(igs$tabestimates, p_threshold = 0.95),
               sampling_ranges = data.frame(min = c(0.1, 1), max = c(10, 100), row.names = c("k_pop", "v_pop"))),
          file.path(root, "step1_igs", "results.rds"))
  cfg <- list(analysis = list(name = "test", out_folder = root), step1_igs = list(narrow = list(p_threshold = 0.95)))
  page <- suppressWarnings(render_igs_figures(cfg, 1))
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  n_figs <- 2 * length(igs_figures(igs)) + 1
  expect_equal(lengths(regmatches(html, gregexpr("<img role=\"img\" src=\"data:image/png", html, fixed = TRUE))), n_figs)
  expect_false(grepl("Could not draw", html))
  expect_match(html, "Wide IGS[\\s\\S]*Narrow ranges[\\s\\S]*Narrow IGS", perl = TRUE)
})
