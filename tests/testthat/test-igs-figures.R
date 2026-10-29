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

# --- step 3 and 4 pages -------------------------------------------------------------------

# A config and saved step 3 results (the SCM fixture, two candidates) in a temporary folder
saved_step3 <- function(candidates = TRUE) {
  cfg <- read_analysis_config(make_config())
  paths <- pipeline_paths(cfg)
  dir.create(paths$step3_search, recursive = TRUE)
  mb <- readRDS(testthat::test_path("fixtures", "scm_results_synthetic.rds"))
  saveRDS(list(SCM = mb), file.path(paths$step3_search, "results.rds"))
  if (candidates) {
    cdir <- paths$step3_candidates
    dir.create(cdir, recursive = TRUE)
    cands <- structure(list(scm_best = list(ta = "vac"), model1 = list(NPs = "tlr7")),
                       duplicates = data.frame(dropped = "scm_delta2", same_as = "scm_best"))
    sf <- data.frame(model = c("scm_best", "model1"), covariate_model = c("ta: vac", "NPs: tlr7"), bicc = c(100, 104),
                     delta_bicc = c(0, 4), condition_number = c(10, 30), wald_note = c("", "beta_NPs_tlr7_1 (p = 0.2)"),
                     source = c("copied from the search", "fitted from step 2's best values"))
    saveRDS(list(candidates = cands, single_fits = sf), file.path(cdir, "results.rds"))
    for (m in names(cands)) {
      dir.create(file.path(cdir, m))
      utils::write.csv(data.frame(parameter = c("ta_pop", "omega_ta", "a"), value = c(1.5, 0.3, 0.2), rse_sa = c(10, 25, NA)),
                       file.path(cdir, m, "populationParameters.txt"), row.names = FALSE)
    }
  }
  cfg
}

test_that("each candidate's single-fit estimates go in one table, with RSEs", {
  cfg <- saved_step3()
  est <- candidate_estimates(pipeline_paths(cfg)$step3_candidates, c("scm_best", "model1", "missing"))
  expect_equal(names(est), c("parameter", "scm_best", "model1"))
  expect_equal(est$parameter, c("ta_pop", "omega_ta", "a"))
  expect_equal(est$scm_best, c("1.5 (10%)", "0.3 (25%)", "0.2"))   # no RSE: no brackets
})

test_that("the search settings are described, including LRT thresholds and fixed relationships", {
  cfg <- read_analysis_config(make_config(list(
    "  criterion: BIC" = "  criterion: LRT",
    "  relationships: ~" = "  relationships: {include: {ta: [vac]}, exclude: {m: [tlr7]}}")))
  txt <- search_settings_text(cfg, "SCM")
  expect_match(txt, "likelihood-ratio test \\(add a covariate at p < 0.05, remove it at p > 0.01\\)")
  expect_match(txt, "In every model: ta: vac")
  expect_match(txt, "In no model: m: tlr7")
  expect_match(search_settings_text(cfg, character()), "No covariate search")
})

test_that("the step 3 page has the search, scm_delta2's path and the candidate fits", {
  skip_if_not(rmarkdown::pandoc_available())
  cfg <- saved_step3()
  page <- suppressWarnings(render_model_building_page(cfg))
  expect_equal(basename(page), "model_building.html")
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  for (x in c("Search settings", "SCM search", "Tested models", "Forward path for scm_delta2", "Single fits",
              "Same covariate model, fitted once: scm_delta2 = scm_best", "beta_NPs_tlr7_1")) {
    expect_true(grepl(x, html, fixed = TRUE), info = x)
  }
  expect_false(grepl("Could not draw", html))
})

test_that("before checkpoint 2, the step 3 page says the candidates aren't fitted yet", {
  skip_if_not(rmarkdown::pandoc_available())
  html <- paste(readLines(suppressWarnings(render_model_building_page(saved_step3(candidates = FALSE))), warn = FALSE), collapse = "\n")
  expect_match(html, "Not fitted yet")
})

test_that("the step 4 page has the comparison tables, plots and a section per model", {
  skip_if_not(rmarkdown::pandoc_available())
  cfg <- read_analysis_config(make_config())
  dir <- pipeline_paths(cfg)$step4
  dir.create(dir, recursive = TRUE)
  a <- fake_igs_full(6)
  b <- fake_igs_full(6)
  b$tabestimates$bicc <- b$tabestimates$bicc + 3
  igs <- list(m1 = a, m2 = b)
  saveRDS(list(igs = igs, comparison = suppressMessages(compare_models(igs)), initial_values = data.frame(run = 1:6)),
          file.path(dir, "results.rds"))
  page <- suppressWarnings(render_comparison_page(cfg))
  expect_equal(basename(page), "comparison.html")
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  for (x in c("2 models (m1, m2)", "Summary", "Pairwise", "BICc per run", "Model m1", "Model m2", "Normalized estimates")) {
    expect_true(grepl(x, html, fixed = TRUE), info = x)
  }
  expect_false(grepl("Could not draw", html))
})
