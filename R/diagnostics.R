#' Figures of an IGS run
#'
#' The standard diagnostic figures of one initial guess sampling run, as a
#' list of figures to draw: used by [save_igs_diagnostics()] (one PDF each)
#' and by [render_igs_figures()] (one HTML page per step). Plots about the
#' consistency of the estimates (convergence, normalized boxplots, CV, paired
#' estimates, correlations between estimates) show only the best-supported
#' runs ([top_runs()]); the others show every run, coloured by BICc. Each
#' caption says which runs it shows and what the marks mean. With fewer than
#' 3 best-supported runs, the consistency plots show every run instead.
#'
#' @param igs Output of [initial_guess_sampling()].
#' @param p_threshold Cumulative BICc weight defining the best-supported runs.
#' @param remove_bicc_outliers,remove_estimate_outliers Passed to the plots
#'   that support them (`step1_igs.plots` in the config); the caption says
#'   what was removed.
#' @param sets_ranges `TRUE` for step 1's wide IGS, whose best-supported runs
#'   set the narrow ranges (said in the captions).
#' @return Named list; each figure has `title`, `file`, `width`, `height` and
#'   `draw` (a function returning the plot). Figures that need data the IGS
#'   didn't record (SAEM iterations, individual fits) are left out.
#' @export
igs_figures <- function(igs, p_threshold = 0.95, remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE,
                        sets_ranges = FALSE) {
  te <- igs$tabestimates
  iv <- igs$initialvalues
  rb <- remove_bicc_outliers
  re <- remove_estimate_outliers
  top <- top_runs(te, p_threshold)
  all_note <- runs_note(te)
  if (length(top) >= 3) {
    te_top <- te[te$run %in% top, ]
    iv_top <- iv[iv$run %in% top, ]
    top_note <- runs_note(te, top, p_threshold, sets_ranges)
  } else {
    te_top <- te
    iv_top <- iv
    top_note <- paste0(all_note, " (only ", length(top), " run(s) make up ", signif(100 * p_threshold, 3),
                       "% of the BICc weight, too few to show alone).")
  }
  threshold <- max(te$bicc[te$run %in% top])
  fig <- function(title, file, width, height, draw) list(title = title, file = file, width = width, height = height, draw = draw)

  figs <- list(
    bicc_vs_estimates = fig("BICc vs estimates", "BIC_vs_estimates_treatment_effects.pdf", 10, 7, function()
      plot_bicc_v_estimate(te, narrow = FALSE, remove_bicc_outliers = rb, remove_estimate_outliers = re,
                           BICc_threshold = threshold, note = all_note)),
    distributions_pop = fig("Initial values and estimates: population parameters", "iv_estimate_dist_pop.pdf", 10, 7, function()
      plot_distributions(iv, te, parameter_set = "pop", save = FALSE, remove_estimate_outliers = re,
                         remove_bicc_outliers = rb, note = all_note)),
    distributions_omega = fig("Initial values and estimates: random-effect SDs", "iv_estimate_dist_omega.pdf", 10, 7, function()
      plot_distributions(iv, te, parameter_set = "omega", save = FALSE, remove_estimate_outliers = re,
                         remove_bicc_outliers = rb, note = all_note)),
    estimates_vs_initial_values = fig("Estimates vs initial values", "estimates_vs_ivs.pdf", 12, 10, function()
      plot_iv_vs_estimates(te, iv, remove_estimate_outliers = re, remove_bicc_outliers = rb, note = all_note)),
    cor_initial_values = fig("Correlation of initial values and estimates", "cor_iv_estimates.pdf", 7, 8, function()
      cor_iv_estimates(te, iv, remove_estimate_outliers = re, remove_bicc_outliers = rb, note = all_note)),
    saem_trajectories = if (!is.null(igs$tabiters)) fig("SAEM trajectories", "SAEM_trajectories.pdf", 15, 7, function()
      plot_SAEM(igs$tabiters, igs$phase, te, remove_bicc_outliers = rb, note = all_note)),
    convergence_estimates = fig("Estimates and RSEs per run", "Convergence_estimates.pdf", 15, 7, function()
      plot_estimates(te_top, remove_bicc_outliers = rb, note = top_note)),
    estimate_boxplots = fig("Normalized estimates", "normalized_estimate_boxplots.pdf", 12, 7, function()
      estimate_boxplots(te_top, remove_bicc_outliers = rb, remove_estimate_outliers = re, note = top_note)),
    coefficient_of_variation = fig("CV of initial values and estimates", "coefficient_of_variation.pdf", 8.5, 10, function()
      calc_CV(te_top, iv_top, remove_estimate_outliers = re, remove_bicc_outliers = rb, note = top_note)),
    paired_estimates = fig("Paired population estimates", "paired_SAEM_population_estimates.pdf", 7, 7, function()
      compare_estimates(te_top, remove_estimate_outliers = re, remove_bicc_outliers = rb, note = top_note)),
    estimate_correlations = fig("Correlations between estimates", "SAEM_correlations.pdf", 10, 10, function()
      correlate_estimates(te_top, remove_estimate_outliers = re, remove_bicc_outliers = rb, note = top_note)),
    rses_vs_bicc = fig("RSE vs BICc", "rses_vs_bicc.pdf", 10, 7, function() plot_vs_rses(te, note = all_note)),
    individual_fits = if (!is.null(igs$indivFit)) fig("Individual fits", "individual_fits.pdf", 12, 10, function()
      plot_individual_fits(igs$indivFit, te))
  )
  figs[!vapply(figs, is.null, logical(1))]
}

#' Save the standard diagnostic plots and tables of an IGS run
#'
#' One call per IGS step in the driver. Writes one PDF per figure of
#' [igs_figures()], plus `identifiability_summary.csv`, `estimates.csv` and
#' `condition_number.csv`. A plot that fails is reported and skipped; the
#' others are still written.
#'
#' @inheritParams igs_figures
#' @param out_dir Folder to write to (created if needed).
#' @return Invisibly, a data frame with one row per output: `output`, `ok`,
#'   `message`.
#' @export
save_igs_diagnostics <- function(igs, out_dir, p_threshold = 0.95,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE, sets_ranges = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  te <- igs$tabestimates
  tables <- list(
    estimates.csv = function() utils::write.csv(te, file.path(out_dir, "estimates.csv"), row.names = FALSE),
    identifiability_summary.csv = function() utils::write.csv(identifiability_summary(igs, p_threshold),
                                                              file.path(out_dir, "identifiability_summary.csv"), row.names = FALSE),
    condition_number.csv = function() if (!is.null(igs$condition_number))
      utils::write.csv(igs$condition_number, file.path(out_dir, "condition_number.csv"), row.names = FALSE)
  )
  figs <- igs_figures(igs, p_threshold, remove_bicc_outliers, remove_estimate_outliers, sets_ranges)
  steps <- c(tables, lapply(figs, function(f) function() ggsave(file.path(out_dir, f$file), f$draw(),
                                                                 width = f$width, height = f$height)))

  res <- lapply(names(steps), function(nm) {
    msg <- tryCatch({
      grDevices::pdf(NULL)  # keep stray plots from opening a device
      on.exit(grDevices::dev.off(), add = TRUE)
      suppressMessages(steps[[nm]]())
      ""
    }, error = function(e) conditionMessage(e))
    if (nzchar(msg)) cat("Could not write ", nm, ": ", msg, "\n", sep = "")
    data.frame(output = nm, ok = !nzchar(msg), message = msg, stringsAsFactors = FALSE)
  })
  invisible(do.call(rbind, res))
}

#' Plots of the final model's step 4 runs
#'
#' The paired-estimates plot ([compare_estimates()]: population estimates
#' against each other across runs, colored by BICc), the normalized boxplots
#' ([estimate_boxplots()]: each estimate divided by its median, with the CV)
#' and BICc against RSE ([plot_vs_rses()]), from the final model's runs in the
#' step 4 comparison. Fixed
#' parameters (one value in every run) are left out. Doesn't need Monolix.
#'
#' @param tabestimates The final model's `tabestimates` from step 4
#'   (`run_step4()$igs[[final_model]]$tabestimates`).
#' @param dir Folder to write `paired_SAEM_population_estimates.pdf`,
#'   `normalized_estimate_boxplots.pdf` and `rses_vs_bicc.pdf` to.
#' @inheritParams save_igs_diagnostics
#' @return Invisibly, a data frame with one row per plot: `output`, `ok`,
#'   `message`.
#' @export
save_final_model_plots <- function(tabestimates, dir, remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  varies <- tapply(tabestimates$estimate, tabestimates$parameter, function(x) length(unique(signif(x, 10))) > 1)
  te <- tabestimates[tabestimates$parameter %in% names(varies)[varies], ]
  of <- dirname(dir); pn <- basename(dir)
  rb <- remove_bicc_outliers; re <- remove_estimate_outliers
  steps <- list(
    paired_estimates = function() compare_estimates(te, of, pn, remove_estimate_outliers = re, remove_bicc_outliers = rb),
    estimate_boxplots = function() estimate_boxplots(te, remove_bicc_outliers = rb, remove_estimate_outliers = re,
                                                     out_folder = of, project_name = pn),
    rses_vs_bicc = function() plot_vs_rses(te, of, pn))
  res <- lapply(names(steps), function(nm) {
    msg <- tryCatch({
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
      suppressMessages(steps[[nm]]())
      ""
    }, error = function(e) conditionMessage(e))
    if (nzchar(msg)) cat("Could not write ", nm, ": ", msg, "\n", sep = "")
    data.frame(output = nm, ok = !nzchar(msg), message = msg, stringsAsFactors = FALSE)
  })
  invisible(do.call(rbind, res))
}

#' One HTML page with every IGS figure of a step
#'
#' Redraws the figures of [igs_figures()] from a step's saved results (no
#' Monolix needed) and writes them to one page, `igs_figures.html` in the
#' step's folder. Step 1's page has the wide IGS, the narrow ranges chosen from
#' it, the narrow IGS and any profile likelihoods; step 2's page has its IGS.
#' The pipeline writes the page at the end of step 1, after profile
#' likelihoods and at the end of step 2; call this to remake it, e.g. for
#' results made before the page existed. The separate PDFs are unchanged.
#'
#' @param config Output of [read_analysis_config()].
#' @param step 1 or 2.
#' @return Path to the page, invisibly.
#' @export
render_igs_figures <- function(config, step = 1) {
  stopifnot(step %in% 1:2)
  paths <- pipeline_paths(config)
  dir <- if (step == 1) paths$step1 else paths$step2
  rds <- file.path(dir, "results.rds")
  if (!file.exists(rds)) stop("No saved results for step ", step, " in ", dir, call. = FALSE)
  res <- readRDS(rds)
  p <- config$step1_igs$narrow$p_threshold %||% 0.95
  pl <- config$step1_igs$plots
  figs <- function(igs, sets_ranges = FALSE) unname(igs_figures(igs, p, isTRUE(pl$remove_bicc_outliers),
                                                               isTRUE(pl$remove_estimate_outliers), sets_ranges))
  block <- function(title, figures) list(title = title, figures = figures)

  if (step == 1) {
    ranges <- list(title = "Initial values and estimates of the wide and narrow IGS", width = 12, height = 6, draw = function()
      plot_parameter_ranges(res$sampling_ranges, res$wide$tabestimates, res$narrow_ranges, res$narrow$tabestimates,
                            p_threshold = p, report_transforms = config$model$report_transforms))
    blocks <- list(block("Wide IGS", figs(res$wide, sets_ranges = TRUE)),
                   block("Narrow ranges", list(ranges)),
                   block("Narrow IGS", figs(res$narrow)))
    prof <- saved_profiles(dir, res$narrow$tabestimates)
    if (!is.null(prof)) blocks <- c(blocks, list(block("Profile likelihoods", list(prof))))
  } else {
    blocks <- list(block("Step 2 IGS", figs(res)))
  }

  fig_dir <- tempfile("igs_figures_")
  dir.create(fig_dir)
  on.exit(unlink(fig_dir, recursive = TRUE), add = TRUE)
  template <- system.file("rmarkdown", "templates", "igs_figures", "skeleton", "skeleton.Rmd", package = "monolixPipelineR")
  out <- rmarkdown::render(
    template,
    params = list(title = paste0(config$analysis$name, ": step ", step, " IGS figures"), blocks = blocks, fig_dir = fig_dir),
    output_file = "igs_figures.html", output_dir = normalizePath(dir),
    intermediates_dir = tempdir(), envir = new.env(parent = asNamespace("monolixPipelineR")), quiet = TRUE)
  invisible(out)
}

# The profile likelihoods saved after step 1, as one figure (NULL if there are none)
saved_profiles <- function(step1_dir, narrow_estimates) {
  dir <- file.path(step1_dir, "profile_likelihood")
  log_file <- file.path(dir, "run_log.yaml")
  if (!file.exists(log_file)) return(NULL)
  done <- yaml::read_yaml(log_file, handlers = yaml_handlers())$extra[["profile_likelihood"]]
  done <- done[vapply(names(done), function(p) file.exists(file.path(dir, paste0(p, ".csv"))), logical(1))]
  if (!length(done)) return(NULL)
  profiles <- lapply(stats::setNames(names(done), names(done)), function(p) utils::read.csv(file.path(dir, paste0(p, ".csv"))))
  log10 <- names(done)[vapply(done, function(pl) isTRUE(pl$log10), logical(1))]
  best <- best_run_estimates(narrow_estimates)
  list(title = "Profile likelihoods", width = 5 * min(length(profiles), 3), height = 4 * ceiling(length(profiles) / 3),
       draw = function() plot_profile_likelihood(profiles, log10 = log10, best = stats::setNames(best$estimate, best$parameter)))
}

# Write a step's figures page; a failure (e.g. no pandoc) is reported, never stops the pipeline
figures_page <- function(config, step) {
  tryCatch({
    page <- render_igs_figures(config, step)
    message("Step ", step, ": all IGS figures on one page: ", page)
  }, error = function(e) message("Could not write the step ", step, " figures page: ", conditionMessage(e)))
  invisible(NULL)
}
