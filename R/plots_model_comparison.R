# Plots for the step 4 comparison (from compare_models()), generalized from the
# model-comparison script of the SIV analysis.

#' BICc of each candidate model across paired runs
#'
#' @param cmp Output of [compare_models()].
#' @return A ggplot: one point per run and model, colored by whether the model
#'   was the best in that run, within `delta` of it, or neither.
#' @export
plot_model_comparison <- function(cmp) {
  d <- cmp$per_run %>%
    mutate(status = ifelse(is_best, "Best model", ifelse(within_delta, paste0("Within ", cmp$delta, " of best"), "Not the best")))
  ggplot(d, aes(model, bicc)) +
    geom_violin(color = "gray70") +
    ggbeeswarm::geom_beeswarm(aes(color = status), alpha = 0.75, size = 2) +
    scale_color_manual(values = stats::setNames(c("blue", "darkorchid", "darkgray"),
                                                c("Best model", paste0("Within ", cmp$delta, " of best"), "Not the best"))) +
    coord_flip() +
    labs(x = "Model", y = "BICc", color = NULL,
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
  ggplot(cmp$wald, aes(model, log10(p_adjusted))) +
    ggbeeswarm::geom_beeswarm(color = "darkgray") +
    geom_hline(aes(yintercept = log10(0.05), linetype = "p = 0.05")) +
    geom_hline(aes(yintercept = log10(0.01), linetype = "p = 0.01")) +
    facet_wrap(~parameter) +
    labs(x = "Model", y = "log10 adjusted p-value", linetype = NULL)
}
