# Helpers for the HTML report (inst/rmarkdown/templates/monolix_report). They take plain
# tables, so they're tested without Monolix; the template reads the tables from Monolix.

#' Group of each individual, from the report's stratifying covariates
#'
#' Replaces the SIV-specific treatment groups (tlr7/vac/ab/nnrti) in the
#' report. A group key is the individual's values of `stratify_by` joined
#' with "." in `stratify_by` order, e.g. `"1.0"` for `tlr7 = 1, vac = 0`.
#' `group_labels` (`report.group_labels` in the config) names keys; unnamed
#' keys are shown as `tlr7=1, vac=0`.
#'
#' @param data The analysis data set.
#' @param id_column Name of the id column in `data`.
#' @param stratify_by Covariate columns of `data` (`report.stratify_by`).
#'   `NULL` puts everyone in one group, "All".
#' @param group_labels Named list or vector, group key -> label.
#' @return Data frame with `id`, `group_key` and `group` (a factor, levels in
#'   key order).
#' @export
report_groups <- function(data, id_column = "id", stratify_by = NULL, group_labels = NULL) {
  stratify_by <- unlist(stratify_by)
  if (!id_column %in% names(data)) stop("No id column '", id_column, "' in the data", call. = FALSE)
  missing <- setdiff(stratify_by, names(data))
  if (length(missing)) stop("report.stratify_by column(s) not in the data: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!length(stratify_by)) {
    ids <- unique(data[[id_column]])
    return(data.frame(id = ids, group_key = "all", group = factor("All")))
  }
  covs <- unique(data[, c(id_column, stratify_by), drop = FALSE])
  if (anyDuplicated(covs[[id_column]])) {
    stop("report.stratify_by column(s) change within an individual: ", paste(stratify_by, collapse = ", "), call. = FALSE)
  }
  key <- do.call(paste, c(unname(as.list(covs[stratify_by])), sep = "."))
  default <- do.call(paste, c(lapply(stratify_by, function(s) paste0(s, "=", covs[[s]])), sep = ", "))
  labels <- unlist(group_labels)
  label <- ifelse(key %in% names(labels), labels[key], default)
  ord <- order(key)
  data.frame(id = covs[[id_column]], group_key = key,
             group = factor(label, levels = unique(label[ord])))
}

#' Rescale population estimates for reporting
#'
#' Applies `model.report_transforms` (per parameter a `label` and an `expr` in
#' `x`, e.g. `{NPs_pop: {label: NP_pop, expr: "x * 10^5"}}`) to the estimate
#' and any confidence limits. Replaces the hard-coded `betas_pop` and
#' `NPs_pop` scaling in the SIV report.
#'
#' @param estimates Data frame with `parameter` and `estimate`, and optionally
#'   `P2.5` and `P97.5`.
#' @param report_transforms Named list from the config.
#' @return `estimates` with `reported_parameter` and `reported_estimate` (and
#'   `reported_P2.5`, `reported_P97.5` when the limits are there). Parameters
#'   without a transform are reported as they are.
#' @export
transform_estimates <- function(estimates, report_transforms = NULL) {
  out <- estimates
  out$reported_parameter <- out$parameter
  cols <- intersect(c("estimate", "P2.5", "P97.5"), names(out))
  for (col in cols) out[[paste0("reported_", col)]] <- out[[col]]
  for (p in intersect(names(report_transforms), out$parameter)) {
    rt <- report_transforms[[p]]
    rows <- out$parameter == p
    if (!is.null(rt$label)) out$reported_parameter[rows] <- rt$label
    for (col in cols) out[[paste0("reported_", col)]][rows] <- eval_transform(rt$expr, out[[col]][rows])
  }
  out
}

# Error-model parameters by name: a, b, c with one observation; a1, b2, ... with several.
# (is_error_param() can't be used here: by-group values have no _pop suffix either.)
is_error_param_name <- function(names) grepl("^[abc][0-9]*$", names)

#' Fine-grained individual and population fits from Monolix's charts data
#'
#' Reads `ChartsData/IndividualFits/<obs>_fits.txt` from a results folder, the
#' curves Monolix draws in its individual fits plot (on the `grid` project
#' setting's time points, not only the observation times). The file is written
#' by `lixoftConnectors::computeChartsData(plot = "indfits")`, with columns
#' `ID`, `time`, `pop`, `popPred`, `indivPredMean` and `indivPredMode` (Monolix 2024R1).
#'
#' @param results_dir The project's results folder.
#' @param obs Observation name; picks `<obs>_fits.txt` when there are several.
#' @return Data frame with `id`, `time`, `popPred` and `indivPred`, and the
#'   attribute `indiv_column` (the Monolix column used for `indivPred`), or
#'   `NULL` when there is no fits file or its columns aren't recognized.
#' @export
read_indfits <- function(results_dir, obs = NULL) {
  files <- list.files(file.path(results_dir, "ChartsData", "IndividualFits"),
                      pattern = "_fits\\.txt$", full.names = TRUE)
  if (!is.null(obs) && any(basename(files) == paste0(obs, "_fits.txt"))) {
    files <- files[basename(files) == paste0(obs, "_fits.txt")]
  }
  if (length(files) != 1) return(NULL)
  fits <- utils::read.csv(files, check.names = FALSE)
  nm <- names(fits)
  id <- nm[tolower(nm) == "id"][1]
  time <- nm[tolower(nm) == "time"][1]
  # popPred: population parameters with each individual's covariates. (Monolix's `pop` column
  # uses the median covariates, like popPred_medianCOV in predictions.txt.)
  pop <- if ("popPred" %in% nm) "popPred" else NA
  indiv <- c(grep("^indiv.*mode", nm, ignore.case = TRUE, value = TRUE),
             grep("^indiv.*mean", nm, ignore.case = TRUE, value = TRUE))[1]
  if (anyNA(c(id, time, pop, indiv))) return(NULL)
  out <- data.frame(id = fits[[id]], time = fits[[time]], popPred = fits[[pop]], indivPred = fits[[indiv]])
  attr(out, "indiv_column") <- indiv
  out
}

#' Individual predictions of all model variables, from tables.txt
#'
#' Monolix writes the variables listed in the model's `table = {...}` (OUTPUT
#' section) to `tables.txt`, at the observation times, as `<variable>_mode`,
#' `<variable>_mean`, `<variable>_popPred`, ... columns.
#'
#' @param results_dir The project's results folder.
#' @param prediction Which prediction to read: `"mode"` (default), `"mean"`,
#'   `"popPred"` or `"SAEM"`. Falls back to `"mean"` when there are no
#'   `_mode` columns.
#' @return Long data frame with `id`, `time`, `variable` (a factor in the
#'   model's order) and `value`, and the attribute `prediction`; `NULL` when
#'   there is no tables.txt or no matching columns.
#' @export
read_model_variables <- function(results_dir, prediction = "mode") {
  f <- file.path(results_dir, "tables.txt")
  if (!file.exists(f)) return(NULL)
  tab <- utils::read.csv(f, check.names = FALSE)
  pick <- function(suffix) grep(paste0("_", suffix, "$"), names(tab), value = TRUE)
  cols <- pick(prediction)
  if (!length(cols) && prediction == "mode") { prediction <- "mean"; cols <- pick(prediction) }
  id <- names(tab)[tolower(names(tab)) == "id"][1]
  if (!length(cols) || is.na(id) || !"time" %in% names(tab)) return(NULL)
  vars <- sub(paste0("_", prediction, "$"), "", cols)
  out <- data.frame(id = rep(tab[[id]], length(cols)), time = rep(tab$time, length(cols)),
                    variable = factor(rep(vars, each = nrow(tab)), levels = vars),
                    value = unlist(tab[cols], use.names = FALSE))
  attr(out, "prediction") <- prediction
  out
}

#' Simulation times for each individual
#'
#' `n_points` equally spaced times from each individual's first to last time
#' in the data, written to the subject-specific output file Simulx reads in
#' [simulate_model_variables()].
#'
#' @param data The analysis data set.
#' @param id_column,time_column Names of the id and time columns in `data`.
#' @param n_points Number of times per individual.
#' @return Data frame with `id` and `time`.
#' @export
simulation_grid <- function(data, id_column, time_column, n_points = 250) {
  ranges <- lapply(split(data[[time_column]], data[[id_column]]), range, na.rm = TRUE)
  ids <- unique(data[[id_column]])
  do.call(rbind, lapply(ids, function(i) {
    r <- ranges[[as.character(i)]]
    data.frame(id = i, time = unique(seq(r[1], r[2], length.out = n_points)))
  }))
}

# Where simulate_model_variables() writes its curves, in the project's results folder
simulated_variables_file <- function(results_dir) {
  file.path(results_dir, "ModelVariables", "simulated_model_variables.txt")
}

#' Model variables simulated with Simulx
#'
#' Reads the curves written by [simulate_model_variables()].
#'
#' @param results_dir The project's results folder.
#' @return Long data frame with `id`, `time`, `variable` (factor, in the
#'   simulated order) and `value`, with the attribute `prediction` ("mode"
#'   for the EBEs, "mean" for conditional means); `NULL` when there is no file.
#' @export
read_simulated_variables <- function(results_dir) {
  f <- simulated_variables_file(results_dir)
  if (!file.exists(f)) return(NULL)
  sims <- utils::read.csv(f)
  if (!nrow(sims)) return(NULL)
  out <- data.frame(id = sims$id, time = sims$time,
                    variable = factor(sims$variable, levels = unique(sims$variable)), value = sims$value)
  attr(out, "prediction") <- c(mlx_EBEs = "mode", mlx_CondMean = "mean")[[sims$parameters[1]]]
  out
}
