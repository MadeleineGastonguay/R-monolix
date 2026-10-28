# Replaces the SIV plot_parameter_ranges() and plot_parameter_ranges_v2(): one function, with
# labels and scaling from model.report_transforms instead of hardcoded NPs/betas.

#' Sampling and estimate ranges across the two IGS phases
#'
#' For each `_pop` parameter, one bar per row: the wide sampling range, the
#' estimates of the wide run, the narrow sampling range, the estimates of the
#' narrow run and those of its best-supported runs. A point marks the best
#' estimate.
#'
#' @param wide_ranges Wide sampling ranges (row names = parameters, columns
#'   `min`, `max`, as from [sampling_ranges()]).
#' @param wide_estimates `tabestimates` of the wide run.
#' @param narrow_ranges Narrow sampling ranges ([get_narrow_ranges()] output);
#'   computed from `wide_estimates` with `p_threshold` when `NULL`.
#' @param narrow_estimates `tabestimates` of the narrow run, optional.
#' @param p_threshold Cumulative BICc weight for the best-supported runs.
#' @param report_transforms `model.report_transforms` from the config: per
#'   parameter a `label` and an `expr` in `x` used to rescale values.
#' @param log_scale Log10 x axis (only used when every value is positive).
#' @return A ggplot.
#' @export
plot_parameter_ranges <- function(wide_ranges, wide_estimates, narrow_ranges = NULL, narrow_estimates = NULL,
                                  p_threshold = 0.95, report_transforms = NULL, log_scale = TRUE) {
  span <- function(te, label) {
    te %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate), .groups = "drop") %>%
      mutate(range = label)
  }
  wr <- data.frame(parameter = rownames(wide_ranges), min = wide_ranges[, 1], max = wide_ranges[, 2],
                   range = "Wide: initial values")
  if (is.null(narrow_ranges)) narrow_ranges <- get_narrow_ranges(wide_estimates, p_threshold = p_threshold)
  df <- rbind(wr, span(wide_estimates, "Wide: estimates"),
              data.frame(parameter = narrow_ranges$parameter, min = narrow_ranges$min, max = narrow_ranges$max,
                         range = "Narrow: initial values"))
  final <- wide_estimates
  if (!is.null(narrow_estimates)) {
    top <- calculate_BICw(narrow_estimates)
    top_runs <- top$run[top$cumprob <= p_threshold | seq_len(nrow(top)) == 1]
    df <- rbind(df, span(narrow_estimates, "Narrow: estimates"),
                span(narrow_estimates[narrow_estimates$run %in% top_runs, ], "Narrow: best-supported runs"))
    final <- narrow_estimates
  }
  best <- best_run_estimates(final)[, c("parameter", "estimate")]
  df <- df[grepl("_pop$", df$parameter), ]
  best <- best[best$parameter %in% df$parameter, ]

  df$label <- df$parameter
  best$label <- best$parameter
  for (p in intersect(names(report_transforms), df$parameter)) {
    rt <- report_transforms[[p]]
    f <- function(x) eval_transform(rt$expr, x)
    i <- df$parameter == p
    df$min[i] <- f(df$min[i]); df$max[i] <- f(df$max[i]); df$label[i] <- rt$label %||% p
    j <- best$parameter == p
    best$estimate[j] <- f(best$estimate[j]); best$label[j] <- rt$label %||% p
  }
  lev <- c("Wide: initial values", "Wide: estimates", "Narrow: initial values", "Narrow: estimates",
           "Narrow: best-supported runs")
  df$range <- factor(df$range, levels = rev(lev))
  best$range <- factor(if (is.null(narrow_estimates)) "Wide: estimates" else "Narrow: best-supported runs", levels = rev(lev))

  plot <- ggplot_report(df, aes(y = range, xmin = min, xmax = max, color = range)) +
    geom_linerange(linewidth = 2) +
    geom_point(data = best, aes(x = estimate, y = range), inherit.aes = FALSE, size = 2) +
    facet_wrap(~label, scales = "free_x") +
    scale_color_manual(values = c("Wide: initial values" = "gray60", "Wide: estimates" = "gray30",
                                  "Narrow: initial values" = "darkred", "Narrow: estimates" = "steelblue",
                                  "Narrow: best-supported runs" = "darkgreen"), guide = "none") +
    labs(x = NULL, y = NULL, caption = paste0("Point: best estimate. Best-supported runs: cumulative BICc weight <= ",
                                              p_threshold, "."))
  if (log_scale && all(c(df$min, df$max) > 0)) plot <- plot + scale_x_log10()
  plot
}
