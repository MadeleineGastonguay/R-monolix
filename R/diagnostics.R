#' Save the standard diagnostic plots and tables of an IGS run
#'
#' One call per IGS step in the driver. Writes PDFs with the ported plotting
#' functions, plus `identifiability_summary.csv` and `estimates.csv`. A plot that
#' fails is reported and skipped; the others are still written.
#'
#' @param igs Output of [initial_guess_sampling()].
#' @param out_dir Folder to write to (created if needed).
#' @param p_threshold Cumulative BICc weight for the identifiability summary.
#' @param remove_bicc_outliers,remove_estimate_outliers Passed to the plots
#'   that support them (`step1_igs.plots` in the config); the plot caption says
#'   what was removed.
#' @return Invisibly, a data frame with one row per output: `output`, `ok`,
#'   `message`.
#' @export
save_igs_diagnostics <- function(igs, out_dir, p_threshold = 0.95,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # the ported plots save to file.path(out_folder, project_name)
  of <- dirname(out_dir)
  pn <- basename(out_dir)
  te <- igs$tabestimates
  iv <- igs$initialvalues
  rb <- remove_bicc_outliers
  re <- remove_estimate_outliers

  steps <- list(
    estimates.csv = function() utils::write.csv(te, file.path(out_dir, "estimates.csv"), row.names = FALSE),
    identifiability_summary.csv = function() utils::write.csv(identifiability_summary(igs, p_threshold),
                                                              file.path(out_dir, "identifiability_summary.csv"), row.names = FALSE),
    condition_number.csv = function() if (!is.null(igs$condition_number))
      utils::write.csv(igs$condition_number, file.path(out_dir, "condition_number.csv"), row.names = FALSE),
    bicc_vs_estimates = function() plot_bicc_v_estimate(te, of, pn, narrow = FALSE, remove_bicc_outliers = rb,
                                                        remove_estimate_outliers = re),
    distributions = function() for (set in c("pop", "omega")) plot_distributions(iv, te, of, pn, parameter_set = set,
                                                                                 remove_estimate_outliers = re,
                                                                                 remove_bicc_outliers = rb),
    estimates_vs_initial_values = function() plot_iv_vs_estimates(te, iv, of, pn, remove_estimate_outliers = re,
                                                                  remove_bicc_outliers = rb),
    cor_initial_values = function() cor_iv_estimates(te, iv, of, pn, remove_estimate_outliers = re, remove_bicc_outliers = rb),
    estimate_correlations = function() correlate_estimates(te, out_folder = of, project_name = pn,
                                                           remove_estimate_outliers = re, remove_bicc_outliers = rb),
    paired_estimates = function() compare_estimates(te, of, pn, remove_estimate_outliers = re, remove_bicc_outliers = rb),
    coefficient_of_variation = function() calc_CV(te, iv, of, pn, remove_estimate_outliers = re, remove_bicc_outliers = rb),
    saem_trajectories = function() if (!is.null(igs$tabiters)) plot_SAEM(igs$tabiters, igs$phase, te, of, pn,
                                                                         remove_bicc_outliers = rb),
    convergence_estimates = function() plot_estimates(te, of, pn, remove_bicc_outliers = rb),
    estimate_boxplots = function() estimate_boxplots(te, remove_bicc_outliers = rb, remove_estimate_outliers = re,
                                                     out_folder = of, project_name = pn),
    rses_vs_bicc = function() plot_vs_rses(te, of, pn),
    individual_fits.pdf = function() if (!is.null(igs$indivFit))
      ggsave(file.path(out_dir, "individual_fits.pdf"), plot_individual_fits(igs$indivFit, te), width = 12, height = 10)
  )

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
