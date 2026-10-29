# One HTML page per pipeline step, from its saved results (no Monolix needed): every IGS figure
# (steps 1 and 2), the covariate search and candidate fits (step 3), the model comparison
# (step 4). Pages are blocks (sections) of items: figures, tables and notes. The template is
# inst/rmarkdown/templates/results_page.

page_figure <- function(title, width, height, draw) list(type = "figure", title = title, width = width, height = height, draw = draw)
page_table <- function(title, data, caption = NULL) list(type = "table", title = title, data = data, caption = caption)
page_text <- function(text, title = NULL) list(type = "text", title = title, text = paste(text, collapse = "\n\n"))
page_block <- function(title, items) list(title = title, items = Filter(Negate(is.null), items))

# Numbers to `digits` significant figures, for tables on a page
round_table <- function(d, digits = 4) {
  if (is.null(d)) return(NULL)
  d <- as.data.frame(d)
  d[] <- lapply(d, function(x) if (is.numeric(x)) signif(x, digits) else x)
  d
}

# Render blocks to <dir>/<file>
render_page <- function(blocks, title, dir, file) {
  fig_dir <- tempfile("results_page_")
  dir.create(fig_dir)
  on.exit(unlink(fig_dir, recursive = TRUE), add = TRUE)
  template <- system.file("rmarkdown", "templates", "results_page", "skeleton", "skeleton.Rmd", package = "monolixPipelineR")
  out <- rmarkdown::render(
    template, params = list(title = title, blocks = blocks, fig_dir = fig_dir),
    output_file = file, output_dir = normalizePath(dir),
    intermediates_dir = tempdir(), envir = new.env(parent = asNamespace("monolixPipelineR")), quiet = TRUE)
  invisible(out)
}

read_step_results <- function(dir, step) {
  rds <- file.path(dir, "results.rds")
  if (!file.exists(rds)) stop("No saved results for step ", step, " in ", dir, call. = FALSE)
  readRDS(rds)
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
  res <- read_step_results(dir, step)
  p <- config$step1_igs$narrow$p_threshold %||% 0.95
  pl <- config$step1_igs$plots
  figs <- function(igs, sets_ranges = FALSE) unname(igs_figures(igs, p, isTRUE(pl$remove_bicc_outliers),
                                                               isTRUE(pl$remove_estimate_outliers), sets_ranges))
  if (step == 1) {
    ranges <- page_figure("Initial values and estimates of the wide and narrow IGS", 12, 6, function()
      plot_parameter_ranges(res$sampling_ranges, res$wide$tabestimates, res$narrow_ranges, res$narrow$tabestimates,
                            p_threshold = p, report_transforms = config$model$report_transforms))
    blocks <- list(page_block("Wide IGS", figs(res$wide, sets_ranges = TRUE)),
                   page_block("Narrow ranges", list(ranges)),
                   page_block("Narrow IGS", figs(res$narrow)))
    prof <- saved_profiles(dir, res$narrow$tabestimates)
    if (!is.null(prof)) blocks <- c(blocks, list(page_block("Profile likelihoods", list(prof))))
  } else {
    blocks <- list(page_block("Step 2 IGS", figs(res)))
  }
  render_page(blocks, paste0(config$analysis$name, ": step ", step, " IGS figures"), dir, "igs_figures.html")
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
  page_figure("Profile likelihoods", 5 * min(length(profiles), 3), 4 * ceiling(length(profiles) / 3), function()
    plot_profile_likelihood(profiles, log10 = log10, best = stats::setNames(best$estimate, best$parameter)))
}

#' One HTML page with step 3's covariate search and candidate fits
#'
#' Written to `step3_model_building/model_building.html` from step 3's saved
#' results (no Monolix needed), for checkpoints 2 and 2b. The first part has
#' the search settings and, for each search method, its path, every model it
#' tested and the model it chose (with `scm_delta2` and its forward path when
#' that candidate is used); with SCM and COSSAC, whether they agree. Once the
#' candidates are fitted, the second part has the single-fit table and each
#' candidate's estimates and RSEs. The pipeline writes the page after the
#' search and again after the candidate fits.
#'
#' @inheritParams render_igs_figures
#' @return Path to the page, invisibly.
#' @export
render_model_building_page <- function(config) {
  paths <- pipeline_paths(config)
  searches <- read_step_results(paths$step3_search, 3)
  s3 <- config$step3_model_building

  blocks <- list(page_block("Search settings", list(page_text(search_settings_text(config, names(searches))))))
  for (m in names(searches)) {
    blocks <- c(blocks, list(page_block(paste(m, "search"), search_items(m, searches[[m]], config, paths$step3_search))))
  }
  if (length(searches) > 1) {
    best <- vapply(searches, function(mb) format_covariate_model(scm_best_model(mb)), character(1))
    note <- if (length(unique(best)) == 1) {
      paste0("All searches chose the same model: ", best[[1]], ". It is fitted once, as `",
             tolower(names(searches)[1]), "_best`.")
    } else {
      paste0("The searches chose different models: ", paste0(names(best), ": ", best, collapse = "; "),
             ". Each is a candidate.")
    }
    blocks <- c(blocks, list(page_block(paste(names(searches), collapse = " and "), list(page_text(note)))))
  }

  cdir <- paths$step3_candidates
  if (file.exists(file.path(cdir, "results.rds"))) {
    blocks <- c(blocks, list(page_block("Candidate models", candidate_items(readRDS(file.path(cdir, "results.rds")), cdir))))
  } else {
    blocks <- c(blocks, list(page_block("Candidate models", list(page_text(
      "Not fitted yet: record `step3_model_building.candidate_models` (checkpoint 2), then run step 3 again.")))))
  }
  render_page(blocks, paste0(config$analysis$name, ": step 3, covariate search and candidates"),
              dirname(paths$step3_search), "model_building.html")
}

# Bulleted settings of the covariate search
search_settings_text <- function(config, methods) {
  s3 <- config$step3_model_building
  if (!length(methods)) return("No covariate search (`method: none`): the candidates are the models in `candidate_models`.")
  crit <- if (identical(s3$criterion, "LRT")) {
    paste0("likelihood-ratio test (add a covariate at p < ", s3$lrt_thresholds$forward, ", remove it at p > ",
           s3$lrt_thresholds$backward, ")")
  } else "BICc"
  pairs <- function(x) if (!length(x)) "none" else paste0(names(x), ": ", vapply(x, function(v) paste(unlist(v), collapse = ", "), ""),
                                                           collapse = "; ")
  rel <- s3$relationships
  c(paste0("- Method: ", paste(methods, collapse = ", ")),
    paste0("- Criterion: ", crit),
    paste0("- Covariates searched: ", paste(unlist(s3$covariates), collapse = ", ")),
    paste0("- Parameters searched: ", paste(unlist(s3$parameters), collapse = ", ")),
    paste0("- In every model: ", pairs(rel$include)),
    paste0("- In no model: ", pairs(rel$exclude))) %>% paste(collapse = "\n")
}

# Items for one search method: its path, the models it tested, the model it chose
search_items <- function(m, mb, config, sdir) {
  s3 <- config$step3_model_building
  lrt <- identical(s3$criterion, "LRT")
  tab <- scm_results_table(mb)
  tab$BICc <- round(tab$BICc, 2)
  tab$LL <- round(tab$LL, 2)
  best <- format_covariate_model(scm_best_model(mb))
  items <- list(
    page_figure("Search path", 10, 7, function() plot_model_building(mb)$plot),
    page_table("Tested models", tab, caption = paste0(
      "Every model the search fitted, in order. best_at_step: the model kept at that step",
      if (lrt) " (by the likelihood-ratio test)" else " (lowest BICc)", "; best_model: the model the search chose.")),
    page_text(paste0("**Chosen model (", tolower(m), "_best):** ", best)))
  if (m == "SCM" && "scm_delta2" %in% auto_candidates(config)) {
    delta <- s3$scm_stop_rule_delta_bicc %||% 2
    d2 <- scm_truncated_model(mb, delta)
    path <- attr(d2, "path")
    path$BICc <- round(path$BICc, 2)
    path$improvement <- ifelse(is.na(path$improvement), "", format(round(path$improvement, 2), nsmall = 2))
    items <- c(items, list(
      page_text(paste0("**scm_delta2:** ", format_covariate_model(d2), " (the forward path stopped at the first step ",
                       "improving BICc by less than ", delta, ")")),
      page_table("Forward path for scm_delta2", path[, c("n_covariates", "parameter", "covariate", "BICc", "improvement", "kept")],
                 caption = "The best model at each forward step; improvement: BICc drop from the step before."),
      if (!is.null(note <- scm_delta2_backward_note(mb, delta))) page_text(paste0("**Note:** ", note))))
  }
  report <- file.path(sdir, tolower(m), "final_model", paste0(m, "_final_model_report.html"))
  if (file.exists(report)) {
    items <- c(items, list(page_text(paste0("[Monolix report on the chosen model](search/", tolower(m), "/final_model/",
                                            basename(report), ")"))))
  }
  items
}

# Items for the fitted candidates: the single-fit table and each candidate's estimates
candidate_items <- function(res, cdir) {
  sf <- round_table(res$single_fits)
  dup <- attr(res$candidates, "duplicates")
  list(
    if (!is.null(dup) && nrow(dup)) page_text(paste0("Same covariate model, fitted once: ",
                                                     paste0(dup$dropped, " = ", dup$same_as, collapse = "; "), ".")),
    page_table("Single fits", sf, caption = paste(
      "One fit per candidate. delta_bicc: BICc minus the lowest; condition_number: of the correlation matrix of the",
      "estimates (large means poorly determined); wald_note: covariate effects with Wald p >= 0.05; source: copied",
      "from the search's chosen model, or fitted from step 2's best run.")),
    page_table("Estimates (RSE %) of each candidate's single fit", candidate_estimates(cdir, names(res$candidates))))
}

# Estimate (RSE %) per parameter and candidate, from each candidate's populationParameters.txt
candidate_estimates <- function(cdir, models) {
  rows <- lapply(models, function(m) {
    f <- file.path(cdir, m, "populationParameters.txt")
    if (!file.exists(f)) return(NULL)
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    rse_col <- intersect(c("rse_sa", "rse_lin"), names(d))
    rse <- if (length(rse_col)) suppressWarnings(as.numeric(d[[rse_col[1]]])) else rep(NA_real_, nrow(d))
    data.frame(parameter = d$parameter, model = m,
               value = paste0(signif(d$value, 3), ifelse(is.finite(rse), paste0(" (", round(rse), "%)"), "")),
               stringsAsFactors = FALSE)
  })
  long <- do.call(rbind, rows)
  if (is.null(long)) return(NULL)
  wide <- as.data.frame(tidyr::pivot_wider(long, names_from = model, values_from = value, values_fill = ""))
  wide[order(match(wide$parameter, parameter_levels(long))), , drop = FALSE]
}

#' One HTML page with step 4's model comparison
#'
#' Written to `step4_model_comparison/comparison.html` from step 4's saved
#' results (no Monolix needed), for checkpoint 3: the summary of how often each
#' model fits best, the pairwise counts, the Wald tests of the covariate
#' effects and the comparison plots, then for each model its covariate model
#' and the consistency of its estimates across its best runs (convergence,
#' normalized estimates, paired estimates, correlations). The pipeline writes
#' it at the end of step 4.
#'
#' @inheritParams render_igs_figures
#' @return Path to the page, invisibly.
#' @export
render_comparison_page <- function(config) {
  paths <- pipeline_paths(config)
  res <- read_step_results(paths$step4, 4)
  cmp <- res$comparison
  models <- names(res$igs)
  p <- config$step1_igs$narrow$p_threshold %||% 0.95
  covs <- tryCatch(readRDS(file.path(paths$step3_candidates, "results.rds"))$candidates, error = function(e) NULL)

  n_start <- length(unique(res$initial_values$run))
  overview <- c(
    if (length(models) == 1) {
      paste0("Only one model (", models, "): this is its IGS, with nothing to compare it with.")
    } else {
      paste0(length(models), " models (", paste(models, collapse = ", "), "), each fitted from the same ", n_start,
             " sets of starting values, so they're compared run by run.")
    },
    paste0(length(cmp$paired_runs), " runs succeeded in every model and are compared",
           if (length(cmp$dropped_runs)) paste0("; left out: run ", paste(cmp$dropped_runs, collapse = ", ")) else "", "."))
  summary_items <- list(
    page_text(overview),
    page_table("Summary", round_table(cmp$summary), caption = paste0(
      "times_best: compared runs where the model has the lowest BICc; times_within_delta: runs within ", cmp$delta,
      " BICc of that run's best model; min_bicc, median_bicc: over the compared runs; median_condition_number: of the ",
      "correlation matrix of the estimates.")),
    if (length(models) > 1) page_table("Pairwise", cmp$pairwise, caption = paste0(
      "runs_better_by_delta: compared runs where `model`'s BICc is more than ", cmp$delta, " lower than `other`'s.")),
    if (!is.null(cmp$wald_summary)) page_table("Wald tests of the covariate effects", round_table(cmp$wald_summary), caption =
      "share_p05, share_p01: share of runs where the effect's adjusted p-value is below 0.05 or 0.01; median_p_adjusted: over runs."))
  blocks <- list(page_block("Comparison", summary_items),
                 page_block("Comparison plots", unname(comparison_figures(cmp))))
  for (m in models) {
    figs <- igs_figures(res$igs[[m]], p)
    figs <- unname(figs[intersect(c("convergence_estimates", "estimate_boxplots", "paired_estimates", "estimate_correlations"),
                                  names(figs))])
    cm <- if (!is.null(covs[[m]])) format_covariate_model(covs[[m]]) else NULL
    blocks <- c(blocks, list(page_block(paste("Model", m), c(
      if (!is.null(cm)) list(page_text(paste0("**Covariate model:** ", cm))), figs))))
  }
  render_page(blocks, paste0(config$analysis$name, ": step 4, model comparison"), paths$step4, "comparison.html")
}

# Write a step's results page; a failure (e.g. no pandoc) is reported, never stops the pipeline
figures_page <- function(config, step) {
  what <- c("all IGS figures", "all IGS figures", "the covariate search and candidates", "the model comparison")[step]
  tryCatch({
    page <- switch(step, render_igs_figures(config, 1), render_igs_figures(config, 2),
                   render_model_building_page(config), render_comparison_page(config))
    message("Step ", step, ": ", what, " on one page: ", page)
  }, error = function(e) message("Could not write the step ", step, " results page: ", conditionMessage(e)))
  invisible(NULL)
}
