# Drawing initial values. Taken out of initial_guess_sampling() so step 4 can draw one shared
# set of initial values for every candidate model with the same code.

# One set of initial values, one per population parameter. Same rules and RNG use as the
# SIV initial_guess_sampling() loop:
# fixed value > error-parameter start > log-uniform (if in log10) > set_variance for omegas and
# error parameters > set_treatment for betas > uniform. A beta without set_treatment or a
# sampling range starts at 0 (no effect), as in set_base(); sampling_ranges() gives betas none.
draw_initial_values <- function(param_names, ranges, log10 = NULL, fixed = NULL, set_error = NULL,
                                set_variance = NULL, set_treatment = NULL) {
  vapply(param_names, function(j) {
    if (j %in% names(fixed)) {
      unlist(fixed)[[j]]
    } else if (is_error_param(j) && j %in% names(set_error)) {
      set_error[[j]]  # same start every run; still estimated
    } else if (j %in% log10) {
      10^(runif(1, log10(ranges[j, 1]), log10(ranges[j, 2])))
    } else if ((grepl("omega", j) || is_error_param(j)) && !is.null(set_variance)) {
      set_variance
    } else if (grepl("beta_", j) && !is.null(set_treatment)) {
      set_treatment
    } else if (grepl("beta_", j) && !j %in% rownames(ranges)) {
      0
    } else {
      if (!j %in% rownames(ranges)) stop("No sampling range for ", j, call. = FALSE)
      runif(1, ranges[j, 1], ranges[j, 2])
    }
  }, numeric(1))
}

#' Draw sets of initial values without fitting
#'
#' Uses the same rules as [initial_guess_sampling()]. In step 4 the result is
#' passed to every candidate model as `predefined_initialvalues`, so the models
#' are compared from the same starting points (a paired comparison).
#'
#' @param ranges Sampling ranges: data frame with row names = parameters and
#'   columns `min`, `max` (as from [sampling_ranges()]).
#' @param param_names Parameters to draw values for. Parameters missing from
#'   `ranges` need a value from `fixed`, `set_error`, `set_variance` or
#'   `set_treatment`.
#' @param n Number of sets.
#' @param log10 Parameters sampled log-uniformly; defaults to the `log10`
#'   attribute of `ranges`.
#' @param fixed Named values of fixed parameters.
#' @param set_error Named starting values of error-model parameters
#'   ([error_starts()]).
#' @param set_variance Starting value for every `omega_*` (and error parameters
#'   without a `set_error` value), instead of sampling.
#' @param set_treatment Starting value for every `beta_*`. Without it, a
#'   `beta_*` with no sampling range starts at 0.
#' @param seed Seed for the draws (`NULL` leaves the RNG alone).
#' @return Data frame with `parameter`, `run`, `initial_value`.
#' @export
sample_initial_values <- function(ranges, param_names, n, log10 = attr(ranges, "log10"), fixed = NULL,
                                  set_error = NULL, set_variance = NULL, set_treatment = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  do.call(rbind, lapply(seq_len(n), function(i) {
    v <- draw_initial_values(param_names, ranges, log10, fixed, set_error, set_variance, set_treatment)
    data.frame(parameter = param_names, run = i, initial_value = unname(v), stringsAsFactors = FALSE)
  }))
}

#' Starting values for a model built on an earlier fit
#'
#' Used before model building (step 3) and for candidate models: population
#' parameters start at the earlier best estimate rounded to `sig_figs`
#' significant figures, random-effect SDs at `omega_value`, treatment effects at
#' `treatment_value`, and error parameters at their earlier estimate.
#'
#' @param best Estimates of the best run: data frame with `parameter` and
#'   `estimate` (e.g. the best run of a `tabestimates` table).
#' @param sig_figs Significant figures for `_pop` values.
#' @param omega_value Starting value for every `omega_*`.
#' @param treatment_value Starting value for every `beta_*`.
#' @return Data frame with `parameter` and `initial_value`.
#' @export
prepare_initial_values <- function(best, sig_figs = 1, omega_value = 1, treatment_value = 0) {
  p <- best$parameter
  v <- best$estimate
  v[grepl("_pop$", p)] <- signif(v[grepl("_pop$", p)], sig_figs)
  v[grepl("^omega", p)] <- omega_value
  v[grepl("^beta_", p)] <- treatment_value
  data.frame(parameter = p, initial_value = v, stringsAsFactors = FALSE)
}

# Best run of a tabestimates table (lowest BICc; first run on ties)
best_run_estimates <- function(tabestimates) {
  run <- tabestimates$run[which.min(tabestimates$bicc)]
  tabestimates[tabestimates$run == run, , drop = FALSE]
}
