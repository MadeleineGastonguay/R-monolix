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
#'   `limits` (bounds of logit-normal parameters, when not (0, 1)),
#'   `bounds` ([parameter_bounds()]), `treatment_effects`, `fixed`,
#'   `estimation`, `scenario`, `initial_ranges`, `omega_range` and
#'   `error_initial` (overrides of Monolix's default starting values for
#'   error-model parameters).
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
  params <- parse_mlxtran(config$model$structural_model)$parameters
  dists <- parse_distributions(config$model$distributions)
  limits <- Filter(function(l) !identical(l, c(0, 1)), dists$limits)
  list(
    model_file = config$model$structural_model,
    data_file = config$data$file,
    header_types = unname(ht[data_cols]),
    error_model = as.list(unlist(config$model$error_model)),
    error_initial = unlist(config$model$error_initial),
    exclude_random_effects = nonempty(unique(c(unlist(config$model$random_effects_excluded),
                                               if (use_fixed) removed_random_effects(config)))),
    distributions = nonempty(dists$type),
    limits = nonempty(limits),
    bounds = parameter_bounds(config$model$distributions, params),
    treatment_effects = nonempty(lapply(treatment_effects, unlist)),
    fixed = if (use_fixed) nonempty(unlist(fixed)) else NULL,
    estimation = config$estimation,
    scenario = config$scenario,
    initial_ranges = initial_ranges_table(config$step1_igs$initial_ranges, params),
    omega_range = config$step1_igs$omega_range
  )
}

# Distribution names in analysis.yaml and Monolix's spelling of them
DISTRIBUTIONS <- c(lognormal = "logNormal", normal = "normal", logitnormal = "logitNormal",
                   probitnormal = "probitNormal")

# model.distributions as a distribution per parameter (`type`, analysis.yaml names) and the
# bounds of each logit-normal parameter (`limits`, (0, 1) unless given). An entry is either a
# name (`param2: logitnormal`) or a map (`param2: {distribution: logitnormal, limits: [0, 100]}`).
parse_distributions <- function(x) {
  x <- x %||% list()
  type <- vapply(x, function(d) {
    d <- if (is.list(d)) d$distribution else d
    if (is.character(d) && length(d) == 1) tolower(d) else NA_character_
  }, character(1))
  limits <- lapply(x[which(type == "logitnormal")], function(d) {
    if (is.list(d) && !is.null(d$limits)) as.numeric(unlist(d$limits)) else c(0, 1)
  })
  list(type = type, limits = limits)
}

#' Allowed values of each structural parameter
#'
#' From `model.distributions`: lognormal parameters (Monolix's default) must
#' be above 0, logit-normal ones between their limits ((0, 1) unless
#' `limits` is given), probit-normal ones between 0 and 1, and normal ones
#' can take any value. Sampling ranges, fixed values and profile grids must
#' lie strictly inside these bounds.
#'
#' @param distributions `model.distributions` from the config, e.g.
#'   `list(param1 = "normal", param2 = list(distribution = "logitnormal", limits = c(0, 100)))`.
#' @param parameters Structural parameter names.
#' @return Data frame with `parameter` (the `_pop` name), `distribution`,
#'   `lower` and `upper`.
#' @export
parameter_bounds <- function(distributions, parameters) {
  d <- parse_distributions(distributions)
  type <- stats::setNames(rep("lognormal", length(parameters)), parameters)
  known <- intersect(names(d$type)[!is.na(d$type)], parameters)
  type[known] <- d$type[known]
  lower <- stats::setNames(ifelse(type == "normal", -Inf, 0), parameters)
  upper <- stats::setNames(ifelse(type %in% c("logitnormal", "probitnormal"), 1, Inf), parameters)
  for (p in intersect(names(d$limits), parameters)) {
    l <- d$limits[[p]]
    if (length(l) == 2 && all(is.finite(l))) { lower[p] <- l[1]; upper[p] <- l[2] }
  }
  data.frame(parameter = paste0(parameters, "_pop"), distribution = unname(type),
             lower = unname(lower), upper = unname(upper), stringsAsFactors = FALSE)
}

# Which of `values` lie outside the open interval (lower, upper)
outside_bounds <- function(values, lower, upper) !(values > lower & values < upper)

# "above 0", "between 0 and 100" or "" (any value): what a parameter's bounds allow
bounds_text <- function(lower, upper) {
  if (is.finite(lower) && is.finite(upper)) paste0("between ", format(lower), " and ", format(upper))
  else if (is.finite(lower)) paste0("above ", format(lower))
  else if (is.finite(upper)) paste0("below ", format(upper))
  else ""
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
                 fixed = spec$fixed, distributions = spec$distributions, limits = spec$limits,
                 estimation = spec$estimation,
                 scenario = spec$scenario)
}
