# Turning the config into the settings for one Monolix project.

#' Settings for one Monolix project
#'
#' Resolves the config into everything needed to build one project: the base
#' model, or the base model plus a candidate's treatment effects. Every
#' candidate shares the data, structural model, error model, variability and
#' fixed parameters; only the covariate model changes.
#'
#' @param config Output of [read_analysis_config()], checked with
#'   [validate_config()].
#' @param treatment_effects Named list of covariates per parameter, e.g.
#'   `list(param1 = "cov1", param2 = c("cov1", "cov2"))`. `NULL` for no covariates.
#' @param use_fixed Apply the checkpoint 1 decisions (`decisions.fixed` and
#'   `decisions.random_effects_removed`)? `FALSE` for step 1.
#' @return A list with `model_file`, `data_file`, `header_types` (in data column
#'   order), `error_model`, `exclude_random_effects` (`model.random_effects_excluded`,
#'   plus `decisions.random_effects_removed` with `use_fixed`), `distributions`,
#'   `treatment_effects`, `fixed`, `estimation`, `scenario`, `initial_ranges`,
#'   `omega_range` and `error_initial` (overrides of Monolix's default
#'   starting values for error-model parameters).
#' @export
model_spec <- function(config, treatment_effects = NULL, use_fixed = TRUE) {
  data_cols <- names(read_monolix_data(config$data$file))
  ht <- unlist(config$data$header_types)
  missing <- setdiff(data_cols, names(ht))
  if (length(missing)) stop("No header type for column(s): ", paste(missing, collapse = ", "), call. = FALSE)

  fixed <- config$decisions$fixed
  if (use_fixed && is.null(fixed)) {
    stop("decisions.fixed is not decided yet (checkpoint 1); use use_fixed = FALSE for step 1.", call. = FALSE)
  }

  nonempty <- function(x) if (length(x)) x else NULL
  list(
    model_file = config$model$structural_model,
    data_file = config$data$file,
    header_types = unname(ht[data_cols]),
    error_model = as.list(unlist(config$model$error_model)),
    error_initial = unlist(config$model$error_initial),
    exclude_random_effects = nonempty(unique(c(unlist(config$model$random_effects_excluded),
                                               if (use_fixed) removed_random_effects(config)))),
    distributions = nonempty(unlist(config$model$distributions)),
    treatment_effects = nonempty(lapply(treatment_effects, unlist)),
    fixed = if (use_fixed) nonempty(unlist(fixed)) else NULL,
    estimation = config$estimation,
    scenario = config$scenario,
    initial_ranges = initial_ranges_table(config$step1_igs$initial_ranges,
                                          parse_mlxtran(config$model$structural_model)$parameters),
    omega_range = config$step1_igs$omega_range
  )
}

# initial_ranges from the YAML as a table, with _pop names
initial_ranges_table <- function(ir, model_parameters) {
  ir <- ir[names(ir) %in% model_parameters]
  if (!length(ir)) return(data.frame(parameter = character(), min = numeric(), max = numeric(), log10 = logical()))
  nm <- names(ir)
  data.frame(
    parameter = paste0(nm, "_pop"),
    min = vapply(ir, function(r) as.numeric(r$min %||% NA), numeric(1)),
    max = vapply(ir, function(r) as.numeric(r$max %||% NA), numeric(1)),
    log10 = vapply(ir, function(r) isTRUE(r$log10), logical(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' Sampling ranges for every population parameter of a project
#'
#' Builds the `initialValue_ranges` table [initial_guess_sampling()] expects:
#' `_pop` parameters from `initial_ranges` and random-effect SDs from
#' `omega_range`. Error-model parameters ([error_starts()]) and treatment
#' effects (`beta_*`) get no range; they always start at the same value.
#'
#' @param spec Output of [model_spec()].
#' @param param_names Population parameter names of the project
#'   (`getPopulationParameterInformation()$name`).
#' @return Data frame with row names = parameters and columns `min`, `max`;
#'   attribute `log10` lists the parameters sampled on a log scale.
#' @export
sampling_ranges <- function(spec, param_names) {
  ir <- spec$initial_ranges
  rows <- param_names[!grepl("^beta_", param_names) & !is_error_param(param_names)]
  out <- data.frame(min = rep(NA_real_, length(rows)), max = NA_real_, row.names = rows)
  have <- intersect(rows, ir$parameter)
  out[have, "min"] <- ir$min[match(have, ir$parameter)]
  out[have, "max"] <- ir$max[match(have, ir$parameter)]
  fill <- rows[is.na(out$min) & grepl("^omega", rows)]
  out[fill, "min"] <- spec$omega_range[1]
  out[fill, "max"] <- spec$omega_range[2]
  still <- rows[is.na(out$min)]
  if (length(still)) stop("No sampling range for: ", paste(still, collapse = ", "), call. = FALSE)
  attr(out, "log10") <- intersect(ir$parameter[ir$log10], rows)
  out
}

#' Starting values for error-model parameters
#'
#' Error-model parameters are estimated by SAEM like any other parameter, but
#' their starting value isn't sampled: every run starts them at Monolix's
#' default (1 for `a` and `c`, 0.3 for `b`) unless `model.error_initial` in
#' the config sets another value.
#'
#' @param spec Output of [model_spec()].
#' @param param_names Population parameter names of the project.
#' @return Named numeric vector, one value per error-model parameter.
#' @export
error_starts <- function(spec, param_names) {
  err <- param_names[is_error_param(param_names)]
  defaults <- c(a = 1, b = 0.3, c = 1)
  out <- vapply(err, function(p) {
    if (p %in% names(spec$error_initial)) return(unname(spec$error_initial[[p]]))
    letter <- substr(p, 1, 1)
    if (letter %in% names(defaults)) unname(defaults[[letter]]) else 1
  }, numeric(1))
  stats::setNames(out, err)
}

#' Create a Monolix project from a model spec
#'
#' @param spec Output of [model_spec()].
#' @export
setup_project <- function(spec) {
  require_monolix()
  create_project(spec$model_file, spec$data_file, spec$header_types, spec$error_model,
                 spec$exclude_random_effects, treatment_effects = spec$treatment_effects,
                 fixed = spec$fixed, distributions = spec$distributions, estimation = spec$estimation,
                 scenario = spec$scenario)
}
