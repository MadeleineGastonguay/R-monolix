# Fits of the best run, from the predictions.txt tables initial_guess_sampling() records
# (`indivFit`, with record_indiv_params = TRUE).

# Prediction columns Monolix writes to predictions.txt
PREDICTION_COLUMN <- "^(popPred|indivPred|indWRes|popWRes|iWRes|npde|pwRes|predictions)"

# The observation column of a predictions.txt table: the one column that isn't id, time, run,
# censored or a prediction
observation_column <- function(pred, observation = NULL) {
  if (!is.null(observation)) {
    if (!observation %in% names(pred)) stop("No column '", observation, "' in the predictions", call. = FALSE)
    return(observation)
  }
  cand <- setdiff(names(pred), c("id", "time", "run", "model", "censored"))
  cand <- cand[!grepl(PREDICTION_COLUMN, cand)]
  if (length(cand) != 1) stop("Can't tell the observation column apart; pass `observation` (candidates: ",
                              paste(cand, collapse = ", "), ")", call. = FALSE)
  cand
}

# The chosen run (default: the run with the lowest BICc, or the first run)
best_run <- function(pred, run = NULL, tabestimates = NULL) {
  if (!is.null(run)) return(run)
  if (!is.null(tabestimates)) return(tabestimates$run[which.min(tabestimates$bicc)])
  if ("run" %in% names(pred)) min(pred$run) else NULL
}

# Rows of the chosen run
pick_run <- function(pred, run = NULL, tabestimates = NULL) {
  run <- best_run(pred, run, tabestimates)
  if ("run" %in% names(pred)) pred <- pred[pred$run == run, ]
  pred
}

# Observation points, open circles for censored observations (the `censored` column that
# initial_guess_sampling() adds from the data)
observation_points <- function(pred, obs, ...) {
  if ("censored" %in% names(pred)) {
    pred$observation <- ifelse(pred$censored %in% TRUE, "Censored", "Observed")
    list(geom_point(data = pred, aes(time, .data[[obs]], shape = observation), ...),
         scale_shape_manual(values = c(Observed = 16, Censored = 1), name = NULL))
  } else {
    list(geom_point(data = pred, aes(time, .data[[obs]]), ...))
  }
}

#' Individual fits of an IGS
#'
#' One panel per individual: the observations (open circles when censored),
#' the individual (conditional mode) and population predictions of the best
#' run, and, with `all_runs`, the individual predictions of every other run as
#' thin grey lines, to show how much the fits vary between initial values.
#' Predictions are at the observation times (predictions.txt).
#'
#' @param indivFit Predictions table (`indivFit` from [initial_guess_sampling()]).
#' @param tabestimates Optional `tabestimates`, used to pick the best run.
#' @param run Run to highlight (default: best run, or the first run).
#' @param observation Observation column (guessed when `NULL`).
#' @param ids Individuals to show (default: all).
#' @param log_y Log10 y axis.
#' @param all_runs Also draw the other runs' individual predictions.
#' @return A ggplot.
#' @export
plot_individual_fits <- function(indivFit, tabestimates = NULL, run = NULL, observation = NULL, ids = NULL,
                                 log_y = FALSE, all_runs = TRUE) {
  all <- as.data.frame(indivFit)
  if (!is.null(ids)) all <- all[all$id %in% ids, ]
  best <- best_run(all, run, tabestimates)
  pred <- pick_run(all, best)
  obs <- observation_column(pred, observation)
  ind <- intersect(c("indivPred_mode", "indivPred_mean"), names(pred))[1]
  lines <- pred %>% select(id, time, any_of(c(ind, "popPred"))) %>%
    pivot_longer(-c(id, time), names_to = "prediction", values_to = "value") %>%
    mutate(prediction = ifelse(prediction == "popPred", "Population (best run)", "Individual (best run)"))
  others <- if (all_runs && "run" %in% names(all)) all[all$run != best, ] else all[0, ]
  styles <- c("Individual (best run)" = "#1B9E77", "Population (best run)" = "#D95F02", "Other runs" = "grey70")
  types <- c("Individual (best run)" = "solid", "Population (best run)" = "dashed", "Other runs" = "solid")
  p <- ggplot_report(pred, aes(time, .data[[obs]]))
  if (nrow(others)) {
    p <- p + geom_line(data = others, aes(time, .data[[ind]], group = run, color = "Other runs",
                                          linetype = "Other runs"), linewidth = 0.3, alpha = 0.6)
  }
  p <- p +
    geom_line(data = lines, aes(time, value, color = prediction, linetype = prediction)) +
    observation_points(pred, obs, size = 1.2) +
    facet_wrap(~id) +
    scale_color_manual(values = styles, name = "Prediction") +
    scale_linetype_manual(values = types, name = "Prediction") +
    labs(y = obs, caption = paste0("Best run: ", best, if (nrow(others)) paste0("; ", length(unique(others$run)), " other runs in grey")))
  if (log_y) p <- p + scale_y_log10()
  p
}

#' Population fit of one run, optionally by group
#'
#' Observations per individual (thin lines) with the population prediction,
#' faceted by the covariates in `stratify_by`.
#'
#' @inheritParams plot_individual_fits
#' @param data The analysis data set (for the covariates in `stratify_by`).
#' @param id_column Name of the id column in `data`.
#' @param stratify_by Covariate columns of `data` to facet by.
#' @return A ggplot.
#' @export
plot_population_fits <- function(indivFit, tabestimates = NULL, run = NULL, observation = NULL,
                                 data = NULL, id_column = "id", stratify_by = NULL, log_y = FALSE) {
  pred <- pick_run(as.data.frame(indivFit), run, tabestimates)
  obs <- observation_column(pred, observation)
  if (length(stratify_by)) {
    if (is.null(data)) stop("`data` is needed to stratify by ", paste(stratify_by, collapse = ", "), call. = FALSE)
    covs <- unique(data[, c(id_column, stratify_by), drop = FALSE])
    names(covs)[1] <- "id"
    pred <- merge(pred, covs, by = "id")
    pred$group <- interaction(lapply(stratify_by, function(s) paste0(s, "=", pred[[s]])), sep = ", ")
  } else {
    pred$group <- "All"
  }
  pop <- pred %>% group_by(group, time) %>% summarise(popPred = stats::median(popPred), .groups = "drop")
  p <- ggplot_report(pred, aes(time, .data[[obs]])) +
    geom_line(aes(group = id), alpha = 0.3) +
    observation_points(pred, obs, size = 0.8, alpha = 0.6) +
    geom_line(data = pop, aes(time, popPred), color = "darkred", linewidth = 1) +
    facet_wrap(~group) +
    labs(y = obs, caption = "Red: population prediction (median over individuals at each time)")
  if (log_y) p <- p + scale_y_log10()
  p
}
