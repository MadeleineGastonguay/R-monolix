# Fits of the best run, from the predictions.txt tables initial_guess_sampling() records
# (`indivFit`, with record_indiv_params = TRUE).

# Prediction columns Monolix writes to predictions.txt
PREDICTION_COLUMN <- "^(popPred|indivPred|indWRes|popWRes|iWRes|npde|pwRes|predictions)"

# The observation column of a predictions.txt table: the one column that isn't id, time, run or
# a prediction
observation_column <- function(pred, observation = NULL) {
  if (!is.null(observation)) {
    if (!observation %in% names(pred)) stop("No column '", observation, "' in the predictions", call. = FALSE)
    return(observation)
  }
  cand <- setdiff(names(pred), c("id", "time", "run", "model"))
  cand <- cand[!grepl(PREDICTION_COLUMN, cand)]
  if (length(cand) != 1) stop("Can't tell the observation column apart; pass `observation` (candidates: ",
                              paste(cand, collapse = ", "), ")", call. = FALSE)
  cand
}

# Rows of the chosen run (default: the run with the lowest BICc)
pick_run <- function(pred, run = NULL, tabestimates = NULL) {
  if (is.null(run)) {
    run <- if (!is.null(tabestimates)) tabestimates$run[which.min(tabestimates$bicc)] else min(pred$run)
  }
  if ("run" %in% names(pred)) pred <- pred[pred$run == run, ]
  pred
}

#' Individual fits of one run
#'
#' Observations with the individual (conditional mode) and population
#' predictions, one panel per individual.
#'
#' @param indivFit Predictions table (`indivFit` from [initial_guess_sampling()]).
#' @param tabestimates Optional `tabestimates`, used to pick the best run.
#' @param run Run to plot (default: best run, or the first run).
#' @param observation Observation column (guessed when `NULL`).
#' @param ids Individuals to show (default: all).
#' @param log_y Log10 y axis.
#' @return A ggplot.
#' @export
plot_individual_fits <- function(indivFit, tabestimates = NULL, run = NULL, observation = NULL, ids = NULL,
                                 log_y = FALSE) {
  pred <- pick_run(as.data.frame(indivFit), run, tabestimates)
  obs <- observation_column(pred, observation)
  if (!is.null(ids)) pred <- pred[pred$id %in% ids, ]
  ind <- intersect(c("indivPred_mode", "indivPred_mean"), names(pred))[1]
  lines <- pred %>% select(id, time, any_of(c(ind, "popPred"))) %>%
    pivot_longer(-c(id, time), names_to = "prediction", values_to = "value") %>%
    mutate(prediction = ifelse(prediction == "popPred", "Population", "Individual"))
  p <- ggplot(pred, aes(time, .data[[obs]])) +
    geom_line(data = lines, aes(time, value, linetype = prediction), color = "steelblue") +
    geom_point(size = 1) +
    facet_wrap(~id) +
    scale_linetype_manual(values = c(Individual = "solid", Population = "dashed")) +
    labs(y = obs, linetype = "Prediction")
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
  p <- ggplot(pred, aes(time, .data[[obs]])) +
    geom_line(aes(group = id), alpha = 0.3) +
    geom_point(size = 0.6, alpha = 0.5) +
    geom_line(data = pop, aes(time, popPred), color = "darkred", linewidth = 1) +
    facet_wrap(~group) +
    labs(y = obs, caption = "Red: population prediction (median over individuals at each time)")
  if (log_y) p <- p + scale_y_log10()
  p
}
