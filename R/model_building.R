# Reading the results of a Monolix covariate search (getModelBuildingResults()).
# Each entry is one tested model with `BICc`, `LL`, `tested` (the parameter-covariate pair
# added at that step; c("", "") for the base model), `individualModels` (logical table,
# rows = parameters, columns = covariates: the full covariate model) and `bestModel`.

#' Covariate model of one model-building entry
#'
#' @param individual_models The `individualModels` table of one entry.
#' @return Named list `parameter = covariates`, the same shape as a
#'   `candidate_models` entry in analysis.yaml. Parameters without covariates are
#'   left out; the base model gives an empty list.
#' @export
covariate_model_of <- function(individual_models) {
  m <- as.matrix(individual_models)
  out <- lapply(rownames(m), function(p) colnames(m)[m[p, ] %in% TRUE])
  names(out) <- rownames(m)
  out[lengths(out) > 0]
}

#' Every model tested in a covariate search, as a table
#'
#' @param mb Output of `lixoftConnectors::getModelBuildingResults()`.
#' @return Data frame with one row per tested model: `entry` (position in
#'   `mb`), `parameter` and `covariate` added at that step (empty for the base
#'   model), `BICc`, `LL`, `n_covariates`, `direction` (`base`, `forward` or
#'   `backward`), `best_at_step` (lowest BICc among forward models with the same
#'   number of covariates), `best_model` (Monolix's final model) and
#'   `covariate_model` (the full model, as text).
#' @export
scm_results_table <- function(mb) {
  tab <- data.frame(
    entry = seq_along(mb),
    parameter = vapply(mb, function(x) x$tested[1] %||% "", character(1)),
    covariate = vapply(mb, function(x) x$tested[2] %||% "", character(1)),
    BICc = vapply(mb, function(x) as.numeric(x$BICc), numeric(1)),
    LL = vapply(mb, function(x) as.numeric(x$LL), numeric(1)),
    n_covariates = vapply(mb, function(x) sum(as.matrix(x$individualModels) %in% TRUE), numeric(1)),
    best_model = vapply(mb, function(x) isTRUE(x$bestModel), logical(1)),
    covariate_model = vapply(mb, function(x) format_covariate_model(covariate_model_of(x$individualModels)), character(1)),
    stringsAsFactors = FALSE
  )
  # Forward selection runs until the last model with the most covariates; anything after
  # that removes covariates (backward elimination)
  last_forward <- max(which(tab$n_covariates == max(tab$n_covariates)))
  tab$direction <- ifelse(tab$entry > last_forward, "backward",
                          ifelse(tab$n_covariates == 0 & tab$parameter == "", "base", "forward"))
  fwd <- tab$direction != "backward"
  tab$best_at_step <- FALSE
  for (k in unique(tab$n_covariates[fwd])) {
    rows <- which(fwd & tab$n_covariates == k)
    tab$best_at_step[rows[which.min(tab$BICc[rows])]] <- TRUE
  }
  tab[, c("entry", "direction", "n_covariates", "parameter", "covariate", "BICc", "LL",
          "best_at_step", "best_model", "covariate_model")]
}

#' The final model of a covariate search
#'
#' @inheritParams scm_results_table
#' @return Named list `parameter = covariates` (see [covariate_model_of()]),
#'   with attributes `entry` and `BICc`.
#' @export
scm_best_model <- function(mb) {
  i <- which(vapply(mb, function(x) isTRUE(x$bestModel), logical(1)))
  if (length(i) != 1) stop("Expected exactly one model with bestModel = TRUE, found ", length(i), call. = FALSE)
  structure(covariate_model_of(mb[[i]]$individualModels), entry = i, BICc = as.numeric(mb[[i]]$BICc))
}

#' Stop the forward selection once BICc stops improving enough
#'
#' Walks the forward-selection path of a stepwise covariate search and stops
#' at the first step where the best added covariate improves BICc by less than
#' `delta_bicc` over the previous step's model. Returns that previous model:
#' the `scm_delta2` candidate when `delta_bicc = 2`. Backward-elimination
#' entries are ignored.
#'
#' @inheritParams scm_results_table
#' @param delta_bicc Smallest BICc improvement worth adding a covariate for.
#' @return Named list `parameter = covariates` (see [covariate_model_of()]),
#'   with attributes `entry`, `BICc` and `path` (the best model at each step,
#'   with its improvement over the step before and whether it was kept).
#' @export
scm_truncated_model <- function(mb, delta_bicc = 2) {
  tab <- scm_results_table(mb)
  path <- tab[tab$best_at_step, c("entry", "n_covariates", "parameter", "covariate", "BICc")]
  path <- path[order(path$n_covariates), ]
  path$improvement <- c(NA, -diff(path$BICc))
  stop_at <- which(!is.na(path$improvement) & path$improvement < delta_bicc)[1]
  keep <- if (is.na(stop_at)) nrow(path) else stop_at - 1
  path$kept <- seq_len(nrow(path)) <= keep
  entry <- path$entry[keep]
  structure(covariate_model_of(mb[[entry]]$individualModels), entry = entry,
            BICc = as.numeric(mb[[entry]]$BICc), path = path)
}

# "ta: vac; NPs: tlr7, vac", or "(none)"
format_covariate_model <- function(cm) {
  if (!length(cm)) return("(none)")
  paste(names(cm), vapply(cm, paste, character(1), collapse = ", "), sep = ": ", collapse = "; ")
}
