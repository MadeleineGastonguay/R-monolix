# Profile likelihood plots: change in OFV against the profiled value, one panel per parameter.

#' Plot profile likelihoods
#'
#' One panel per profiled parameter: OFV minus its lowest value on the grid
#' against the profiled value, with a dashed line at `threshold` (3.84, the
#' 95% cutoff of a chi-squared test with one degree of freedom). Values below
#' the line form the approximate 95% confidence interval. A dotted line marks
#' the starting estimate, when given.
#'
#' @param profiles Named list of [profile_likelihood()] tables (one per
#'   parameter), or one table with all of them.
#' @param log10 Names of the parameters whose x axis is on a log10 scale.
#' @param best Optional named numeric vector of starting estimates, e.g. the
#'   best narrow IGS run.
#' @param threshold Change in OFV marking the confidence interval.
#' @return A patchwork of ggplots.
#' @export
plot_profile_likelihood <- function(profiles, log10 = character(), best = NULL, threshold = stats::qchisq(0.95, 1)) {
  pts <- profile_points(profiles)
  panels <- lapply(unique(pts$profiled_parameter), function(p) {
    d <- pts[pts$profiled_parameter == p, ]
    d$delta_ofv <- d$ofv - min(d$ofv)
    g <- ggplot_report(d, aes(profiled_value, delta_ofv)) +
      geom_hline(yintercept = threshold, linetype = "dashed", color = "grey50") +
      geom_line(color = "#1B9E77") +
      geom_point(color = "#1B9E77") +
      labs(x = p, y = "OFV - lowest OFV", title = p)
    if (!is.null(best) && p %in% names(best)) g <- g + geom_vline(xintercept = best[[p]], linetype = "dotted", color = "#D95F02")
    if (p %in% log10) g <- g + scale_x_log10()
    g
  })
  patchwork::wrap_plots(panels) +
    patchwork::plot_annotation(caption = paste0("Dashed: 95% cutoff (", round(threshold, 2), ")",
                                                if (!is.null(best)) "; dotted: starting estimate" else ""))
}

# One row per grid value (the profile tables repeat it for each estimated parameter)
profile_points <- function(profiles) {
  pr <- if (is.data.frame(profiles)) profiles else do.call(rbind, unname(profiles))
  pr <- pr[!is.na(pr$ofv), c("profiled_parameter", "profiled_value", "ofv")]
  pr <- pr[!duplicated(pr[, c("profiled_parameter", "profiled_value")]), ]
  pr[order(pr$profiled_parameter, pr$profiled_value), ]
}
