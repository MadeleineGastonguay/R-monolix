# Plots for the step 4 comparison (from compare_models()), generalized from the
# model-comparison script of the SIV analysis.

#' BICc of each candidate model across paired runs
#'
#' @param cmp Output of [compare_models()].
#' @return A ggplot: one point per run and model, colored by whether the model
#'   was the best in that run, within `delta` of it, or neither.
#' @export
plot_model_comparison <- function(cmp) {
  run_violins(run_status(cmp), "bicc", cmp) + labs(y = "BICc")
}

#' FIM condition number of each candidate model across paired runs
#'
#' Same layout as [plot_model_comparison()]: one point per run and model,
#' colored by whether the model was the best (lowest BICc) in that run. The
#' condition number is the ratio of the largest to the smallest eigenvalue of
#' the correlation matrix of the estimates, on a log10 axis; runs without one
#' (a standard error not computed) are left out.
#'
#' @inheritParams plot_model_comparison
#' @return A ggplot, or `NULL` when no condition numbers were recorded.
#' @export
plot_condition_numbers <- function(cmp) {
  cn <- cmp$condition_numbers
  if (is.null(cn) || !any(is.finite(cn$condition_number))) return(NULL)
  d <- merge(run_status(cmp), cn[is.finite(cn$condition_number), c("model", "run", "condition_number")], by = c("model", "run"))
  d$model <- factor(d$model, levels = levels(cmp$per_run$model))
  n_missing <- nrow(cmp$per_run) - nrow(d)
  run_violins(d, "condition_number", cmp) +
    scale_y_log10() +
    labs(y = "Condition number (log10 scale)",
         caption = paste0(length(cmp$paired_runs), " runs from the same initial values in every model",
                          if (n_missing) paste0("; ", n_missing, " run(s) without a condition number left out") else ""))
}

#' BICc against each estimate, all candidate models overlaid
#'
#' One panel per estimated parameter, with every model's runs in its own
#' color. The y axis is each model's BICc minus that model's lowest BICc, so
#' models whose BICc differ overall still share an axis: a model with a
#' well-defined optimum shows its runs piling up at one estimate near 0.
#' Triangles mark each model's best run. Fixed parameters (one value in every
#' run) are left out.
#'
#' @inheritParams plot_model_comparison
#' @param log10_x Log10 x axis for parameters whose estimates are all positive
#'   (never for treatment effects, which can be negative).
#' @param max_delta_bicc Leave out runs more than this above their model's
#'   best (e.g. failed fits that squash the axis); `NULL` keeps all.
#' @return A patchwork of ggplots, one per parameter, with one legend.
#' @export
plot_bicc_vs_estimate_models <- function(cmp, log10_x = TRUE, max_delta_bicc = NULL) {
  d <- cmp$estimates
  if (is.null(d) || !nrow(d)) stop("cmp has no estimates; it was made by an older compare_models()", call. = FALSE)
  models <- levels(cmp$per_run$model) %||% unique(as.character(d$model))
  d$model <- factor(as.character(d$model), levels = models)
  n_runs <- length(cmp$paired_runs)
  dropped <- 0
  if (!is.null(max_delta_bicc)) {
    keep <- d$delta_bicc <= max_delta_bicc
    dropped <- length(unique(paste(d$model, d$run)[!keep]))
    d <- d[keep, ]
  }
  varies <- tapply(d$estimate, d$parameter, function(x) length(unique(signif(x, 10))) > 1)
  params <- intersect(parameter_levels(d), names(varies)[varies])
  panels <- lapply(params, function(p) {
    dp <- d[d$parameter == p, ]
    g <- ggplot_report(dp, aes(estimate, delta_bicc, fill = model)) +
      geom_point(shape = 21, alpha = 0.7, size = 2, color = "grey30") +
      geom_point(data = dp[dp$model_best, ], shape = 24, size = 3.5, color = "black", show.legend = FALSE) +
      scale_fill_brewer(palette = "Dark2", limits = models, drop = FALSE) +  # same colors in every panel
      guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, alpha = 1))) +
      labs(x = NULL, y = "BICc - model's lowest", fill = "Model", title = p) +
      theme(plot.title = element_text(size = 10, face = "bold"))
    # one legend: patchwork merges only identical ones, so panels missing a model don't show theirs
    g <- g + theme(legend.position = if (all(models %in% dp$model)) "bottom" else "none")
    # narrow ranges get crowded tick labels: drop the ones that would overlap
    g + if (log10_x && !grepl("^beta_", p) && all(dp$estimate > 0)) {
      scale_x_log10(n.breaks = 4, guide = guide_axis(check.overlap = TRUE))
    } else {
      scale_x_continuous(n.breaks = 4, guide = guide_axis(check.overlap = TRUE))
    }
  })
  patchwork::wrap_plots(panels) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "BICc against estimates, all models", theme = theme(legend.position = "bottom"),
      caption = paste0(n_runs, " runs from the same initial values in every model; y = BICc minus that model's lowest BICc; ",
                       "triangles: each model's best run",
                       if (dropped) paste0("\n", dropped, " model run(s) more than ", max_delta_bicc, " above their model's best left out") else ""))
}

# per_run with the status of each model in each run, for the violin plots
run_status <- function(cmp) {
  cmp$per_run %>%
    mutate(status = ifelse(is_best, "Best model", ifelse(within_delta, paste0("Within ", cmp$delta, " of best"), "Not the best")))
}

# Violin and beeswarm of `y` per model, points colored by run status, models on the y axis
run_violins <- function(d, y, cmp) {
  ggplot_report(d, aes(model, .data[[y]])) +
    geom_violin(color = "gray70") +
    ggbeeswarm::geom_beeswarm(aes(color = status), alpha = 0.75, size = 2) +
    scale_color_manual(values = stats::setNames(c("blue", "darkorchid", "darkgray"),
                                                c("Best model", paste0("Within ", cmp$delta, " of best"), "Not the best"))) +
    coord_flip() +
    labs(x = "Model", color = NULL,
         caption = paste0(length(cmp$paired_runs), " runs from the same initial values in every model")) +
    theme(legend.position = "bottom")
}

#' Adjusted Wald p-values of treatment effects across runs
#'
#' @inheritParams plot_model_comparison
#' @return A ggplot, or `NULL` when no Wald tests were recorded.
#' @export
plot_treatment_pvalues <- function(cmp) {
  if (is.null(cmp$wald) || !nrow(cmp$wald)) return(NULL)
  ggplot_report(cmp$wald, aes(model, log10(p_adjusted))) +
    ggbeeswarm::geom_beeswarm(color = "darkgray") +
    geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05")) +
    geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01")) +
    facet_wrap(~parameter) +
    labs(x = "Model", y = "log10 adjusted p-value", linetype = NULL)
}

#' Save the step 4 comparison plots
#'
#' Writes `model_comparison.pdf` (BICc per run), `condition_numbers.pdf`,
#' `bicc_vs_estimates.pdf` (all models overlaid) and `treatment_pvalues.pdf`
#' to `dir`. A plot that can't be drawn is skipped with a message. Doesn't need
#' Monolix, so it also works on saved step 4 results.
#'
#' @inheritParams plot_model_comparison
#' @param dir Folder to write to.
#' @return Paths of the files written, invisibly.
#' @export
save_comparison_plots <- function(cmp, dir) {
  n_par <- if (is.null(cmp$estimates)) 0 else length(unique(cmp$estimates$parameter))
  side <- max(2, ceiling(sqrt(n_par)))
  plots <- list(
    model_comparison = list(function() plot_model_comparison(cmp), 8, 5),
    condition_numbers = list(function() plot_condition_numbers(cmp), 8, 5),
    bicc_vs_estimates = list(function() plot_bicc_vs_estimate_models(cmp), 3.2 * side, 2.8 * side + 1),
    treatment_pvalues = list(function() plot_treatment_pvalues(cmp), 10, 7))
  written <- character()
  for (nm in names(plots)) {
    f <- file.path(dir, paste0(nm, ".pdf"))
    tryCatch({
      g <- plots[[nm]][[1]]()
      if (!is.null(g)) {
        ggsave(f, g, width = plots[[nm]][[2]], height = plots[[nm]][[3]])
        written <- c(written, f)
      }
    }, error = function(e) message("Could not draw ", nm, ": ", conditionMessage(e)))
  }
  invisible(written)
}
