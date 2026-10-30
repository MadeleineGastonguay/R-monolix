# The analysis config (analysis.yaml): defaults, reading and validation.
#
# Conventions:
#   `~` (NULL) for a decision means "not decided yet" -> the driver stops at that checkpoint.
#   `{}` (empty map) means "decided: none", e.g. `fixed: {}` fixes no parameters.

#' Default settings for every key in analysis.yaml
#'
#' These are the values used for any key you leave out (or set to `~`, except
#' for the checkpoint decisions).
#'
#' @return Nested list with the same layout as analysis.yaml.
#' @export
config_defaults <- function() {
  list(
    analysis = list(name = NULL, out_folder = NULL, r_seed = 1),
    monolix = list(version = NULL),
    data = list(file = NULL, header_types = NULL),
    model = list(structural_model = NULL, error_model = NULL, error_initial = list(),
                 random_effects_excluded = list(), distributions = list(), report_transforms = list()),
    estimation = list(
      saem = list(burn_in = 50, exploratory = 1000, smoothing = 500, smoothing_alpha = 0.7,
                  exploratory_autostop = TRUE, smoothing_autostop = TRUE)
    ),
    scenario = list(population_parameters = TRUE, conditional_mode = TRUE, conditional_distribution = TRUE,
                    standard_errors = TRUE, log_likelihood = TRUE, linearization = FALSE, plots = FALSE),
    step1_igs = list(
      initial_ranges = NULL, omega_range = c(1, 5), wide_nits = 100,
      narrow = list(nits = 100, p_threshold = 0.95, min_fold_width = NULL, omega_initial = 1),
      profile_likelihood = list(),
      plots = list(remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE)
    ),
    decisions = list(fixed = NULL, fixed_rationale = NULL, random_effects_removed = list(), fixed_confirmed = NULL),
    step2_igs = list(nits = 100),
    step3_model_building = list(
      method = "SCM", criterion = "BIC", lrt_thresholds = list(forward = 0.05, backward = 0.01),
      covariates = "all", parameters = "all", relationships = NULL, use_linearization = FALSE,
      initial_value_rounding = 1, omega_initial = 1, scm_stop_rule_delta_bicc = 2,
      candidate_models = NULL
    ),
    step4_model_comparison = list(nits = 50, compare = NULL, omega_initial = 1, treatment_initial = 0,
                                  fdr_method = "fdr", final_model = NULL, final_rationale = NULL),
    report = list(stratify_by = NULL, group_labels = NULL)
  )
}

# Random effects removed at checkpoint 1 (applied from step 2 on)
removed_random_effects <- function(config) unlist(config$decisions$random_effects_removed)

#' Does step 2 run?
#'
#' Step 2 refits the base model after checkpoint 1 when anything was fixed
#' (`decisions.fixed`) or any random effect removed
#' (`decisions.random_effects_removed`); otherwise step 1's narrow IGS is used.
#'
#' @param config Output of [read_analysis_config()].
#' @return `TRUE` or `FALSE`.
#' @export
step2_runs <- function(config) {
  length(config$decisions$fixed) > 0 || length(removed_random_effects(config)) > 0
}

# Keys whose values are user-named maps (any keys allowed), not fixed settings
FREE_MAPS <- c("header_types", "error_model", "error_initial", "distributions", "report_transforms", "initial_ranges",
               "profile_likelihood", "fixed", "candidate_models", "group_labels", "relationships")

# Keys holding text that must never be converted to numbers
TEXT_KEYS <- c("name", "out_folder", "file", "structural_model", "version", "header_types", "error_model",
               "distributions", "random_effects_excluded", "random_effects_removed", "method", "criterion", "covariates",
               "parameters", "fdr_method", "final_model", "fixed_rationale", "final_rationale", "label",
               "expr", "stratify_by", "group_labels", "compare", "candidate_models", "relationships")

# Covariate search methods for step3_model_building.method, and the candidate models each one
# adds automatically. scm_delta2 is added only with criterion BIC (see auto_candidates()); COSSAC
# has no stop-rule variant: its path isn't one covariate per step.
MODEL_BUILDING_METHODS <- c("SCM", "COSSAC")
METHOD_CANDIDATES <- list(SCM = c("scm_best", "scm_delta2"), COSSAC = "cossac_best")
# Every name reserved for an automatic candidate
AUTO_CANDIDATES <- unlist(METHOD_CANDIDATES, use.names = FALSE)

#' Covariate search methods to run in step 3
#'
#' `step3_model_building.method` is `SCM`, `COSSAC`, a list of both (each search
#' runs and their results are compared in step 4), or `none` (no search: only
#' the hand-written `candidate_models` are compared). Case doesn't matter.
#'
#' @param config Output of [read_analysis_config()].
#' @return Character vector of methods in upper case (`character()` for none).
#'   Errors on anything else.
#' @export
model_building_methods <- function(config) {
  m <- toupper(unlist(config$step3_model_building$method))
  if (!length(m)) stop("step3_model_building.method is not set", call. = FALSE)
  if (identical(m, "NONE")) return(character())
  bad <- setdiff(m, MODEL_BUILDING_METHODS)
  if (length(bad) || anyDuplicated(m)) {
    stop("step3_model_building.method must be SCM, COSSAC, [SCM, COSSAC] or none; got ",
         paste(unlist(config$step3_model_building$method), collapse = ", "), call. = FALSE)
  }
  m
}

#' Names of the candidate models for step 4
#'
#' The automatic candidates of each covariate search method (`scm_best`, and
#' `scm_delta2` when the criterion is BIC, for SCM; `cossac_best` for COSSAC),
#' then the hand-written `candidate_models`.
#'
#' @inheritParams model_building_methods
#' @return Character vector.
#' @export
candidate_names <- function(config) {
  c(auto_candidates(config), names(config$step3_model_building$candidate_models))
}

# Automatic candidates of the covariate searches that run. scm_delta2 is a BICc stop rule, so
# it's added only with criterion BIC; with LRT the p-value thresholds already decide when to stop.
auto_candidates <- function(config, methods = model_building_methods(config)) {
  out <- unlist(METHOD_CANDIDATES[methods], use.names = FALSE)
  if (!identical(toupper(config$step3_model_building$criterion %||% ""), "BIC")) out <- setdiff(out, "scm_delta2")
  out
}

# Keep y/n/yes/no/on/off as text; only true/false are logicals
yaml_handlers <- function() {
  list("bool#yes" = function(x) if (tolower(x) == "true") TRUE else x,
       "bool#no" = function(x) if (tolower(x) == "false") FALSE else x)
}

#' Read analysis.yaml
#'
#' Fills in defaults from [config_defaults()], converts numeric settings
#' written in forms YAML reads as text (e.g. `1e-7`), and turns the data, model
#' and output paths into absolute paths relative to the YAML file's folder.
#' Use [validate_config()] to check the result.
#'
#' @param path Path to analysis.yaml.
#' @return The config as a nested list.
#' @export
read_analysis_config <- function(path) {
  if (!file.exists(path)) stop("Config file not found: ", path, call. = FALSE)
  user <- yaml::read_yaml(path, handlers = yaml_handlers())
  if (is.null(user)) user <- list()
  config <- fill_defaults(user, config_defaults())
  config <- coerce_numbers(config)

  dir <- dirname(normalizePath(path))
  resolve <- function(p) {
    if (is.null(p) || !nzchar(p)) return(p)
    if (grepl("^(/|~|[A-Za-z]:)", p)) return(normalizePath(p, mustWork = FALSE))
    normalizePath(file.path(dir, p), mustWork = FALSE)
  }
  config$data$file <- resolve(config$data$file)
  config$model$structural_model <- resolve(config$model$structural_model)
  if (is.null(config$analysis$out_folder) && !is.null(config$analysis$name)) {
    config$analysis$out_folder <- file.path("results", config$analysis$name)
  }
  config$analysis$out_folder <- resolve(config$analysis$out_folder)
  attr(config, "path") <- normalizePath(path)
  attr(config, "unknown_keys") <- unknown_keys(user, config_defaults())
  config
}

# Fill keys missing from `user` (or set to ~) with defaults, recursing into fixed-layout sections
fill_defaults <- function(user, defaults) {
  out <- user
  for (k in names(defaults)) {
    d <- defaults[[k]]
    u <- user[[k]]
    if (is.null(u)) {
      out[k] <- list(d)
    } else if (is.list(d) && length(d) > 0 && !is.null(names(d)) && is.list(u) && !k %in% FREE_MAPS) {
      out[[k]] <- fill_defaults(u, d)
    }
  }
  out
}

unknown_keys <- function(user, defaults, prefix = "") {
  out <- character()
  for (k in names(user)) {
    if (!k %in% names(defaults)) {
      out <- c(out, paste0(prefix, k))
    } else if (is.list(user[[k]]) && is.list(defaults[[k]]) && length(defaults[[k]]) > 0 &&
               !is.null(names(defaults[[k]])) && !k %in% FREE_MAPS) {
      out <- c(out, unknown_keys(user[[k]], defaults[[k]], paste0(prefix, k, ".")))
    }
  }
  out
}

coerce_numbers <- function(x, key = "") {
  if (key %in% TEXT_KEYS) return(x)
  if (is.list(x)) {
    nm <- names(x)
    for (i in seq_along(x)) {
      k <- if (is.null(nm)) "" else nm[i]
      if (!is.null(x[[i]])) x[[i]] <- coerce_numbers(x[[i]], k)
    }
    return(x)
  }
  if (is.character(x) && all(grepl("^[-+]?([0-9]*\\.)?[0-9]+([eE][-+]?[0-9]+)?$", x))) {
    return(as.numeric(x))
  }
  x
}

#' Check analysis.yaml for problems
#'
#' Checks only what is needed up to `up_to_step`, so the file can be filled in
#' one step at a time. Every problem found is reported together. If the only
#' thing missing is a checkpoint decision (`decisions.fixed` before step 2,
#' `decisions.fixed_confirmed` before step 3 when step 2 ran ([step2_runs()]),
#' `step3_model_building.candidate_models` before step 3's candidate fits,
#' `step4_model_comparison.compare` before step 4,
#' `step4_model_comparison.final_model` before the final model), the error has
#' class `pipeline_checkpoint` so the driver can stop cleanly.
#'
#' @param config Output of [read_analysis_config()].
#' @param up_to_step 1 to 4.
#' @param decisions_needed If `FALSE`, don't require the decision made partway
#'   through step `up_to_step` (candidate models in step 3, final model in step
#'   4), so the part of the step before that checkpoint can run.
#' @return `config`, invisibly, if everything is fine.
#' @export
validate_config <- function(config, up_to_step = 1, decisions_needed = TRUE) {
  stopifnot(up_to_step %in% 1:4)
  p <- character()        # problems
  todo <- character()     # checkpoint decisions still to make
  add <- function(...) p <<- c(p, paste0(...))

  unk <- attr(config, "unknown_keys")
  if (length(unk)) add("Unknown keys (typo?): ", paste(unk, collapse = ", "))

  # analysis / paths
  if (is_blank(config$analysis$name)) add("REQUIRED: analysis.name")
  if (!is_count(config$analysis$r_seed, zero_ok = TRUE)) add("analysis.r_seed must be a whole number")
  data_ok <- !is_blank(config$data$file) && file.exists(config$data$file)
  if (!data_ok) add("data.file not found: ", config$data$file %||% "(not set)")
  model_ok <- !is_blank(config$model$structural_model) && file.exists(config$model$structural_model)
  if (!model_ok) add("model.structural_model not found: ", config$model$structural_model %||% "(not set)")

  parsed <- if (model_ok) tryCatch(parse_mlxtran(config$model$structural_model),
                                   error = function(e) { add(conditionMessage(e)); NULL }) else NULL
  data <- if (data_ok) tryCatch(read_monolix_data(config$data$file),
                                error = function(e) { add(conditionMessage(e)); NULL }) else NULL
  params <- parsed$parameters
  bounds <- NULL  # allowed values of each parameter, from model.distributions (set below)
  # values of a _pop parameter (a sampling range, fixed value or profile grid) must lie inside its bounds
  check_bounds <- function(key, p, values) {
    b <- if (!is.null(bounds)) bounds[bounds$parameter == p, ]
    if (is.null(b) || !nrow(b) || !is.numeric(values)) return(invisible())
    if (any(outside_bounds(values, b$lower, b$upper), na.rm = TRUE)) {
      add(key, ": ", sub("_pop$", "", p), " has a ", b$distribution, " distribution, so its values must be ",
          bounds_text(b$lower, b$upper), " (see model.distributions)")
    }
  }

  # header types
  ht <- unlist(config$data$header_types)
  obs_cols <- character(); cov_cols <- character()
  if (is.null(ht)) {
    add("REQUIRED: data.header_types")
  } else {
    bad <- setdiff(ht, MONOLIX_HEADER_TYPES)
    if (length(bad)) add("Unknown header types: ", paste(unique(bad), collapse = ", "),
                         " (allowed: ", paste(MONOLIX_HEADER_TYPES, collapse = ", "), ")")
    for (t in c("id", "time")) if (sum(ht == t) != 1) add("data.header_types needs exactly one '", t, "' column")
    obs_cols <- names(ht)[ht == "observation"]
    cov_cols <- names(ht)[ht %in% c("catcov", "contcov")]
    if (length(obs_cols) == 0) add("data.header_types needs at least one 'observation' column")
    if (!is.null(data)) {
      missing_cols <- setdiff(names(data), names(ht))
      extra_cols <- setdiff(names(ht), names(data))
      if (length(missing_cols)) add("data.header_types has no type for column(s): ", paste(missing_cols, collapse = ", "))
      if (length(extra_cols)) add("data.header_types lists column(s) not in the data: ", paste(extra_cols, collapse = ", "))
      id_col <- names(ht)[ht == "id"]
      if (length(id_col) == 1 && id_col %in% names(data)) {
        for (cc in intersect(cov_cols, names(data))) {
          per_id <- tapply(data[[cc]], data[[id_col]], function(v) length(unique(v)))
          if (any(per_id > 1)) add("Covariate '", cc, "' changes within an individual; Monolix needs it constant per id")
        }
      }
      for (cc in intersect(names(ht)[ht == "cens"], names(data))) {
        if (!all(data[[cc]][!is.na(data[[cc]])] %in% c(-1, 0, 1))) add("Censoring column '", cc, "' must be 0, 1 or -1")
      }
    }
  }

  # model
  if (!is.null(parsed) && length(obs_cols) && length(parsed$outputs) != length(obs_cols)) {
    add("The model has ", length(parsed$outputs), " output(s) (", paste(parsed$outputs, collapse = ", "),
        ") but the data has ", length(obs_cols), " observation column(s)")
  }
  em <- unlist(config$model$error_model)
  if (is.null(em)) {
    add("REQUIRED: model.error_model (one entry per observation column)")
  } else {
    bad <- setdiff(em, c("constant", "proportional", "combined1", "combined2"))
    if (length(bad)) add("model.error_model values must be constant, proportional, combined1 or combined2; got ",
                         paste(bad, collapse = ", "))
    if (length(obs_cols) && !setequal(names(em), obs_cols)) {
      add("model.error_model keys (", paste(names(em), collapse = ", "), ") must match the observation column(s) (",
          paste(obs_cols, collapse = ", "), ")")
    }
  }
  if (!is.null(params)) {
    # Monolix refuses a structural parameter named like an error parameter
    n_obs <- max(length(obs_cols), length(em), 1)
    clash <- intersect(params, reserved_error_names(n_obs))
    if (length(clash)) add("Structural parameter(s) ", paste(clash, collapse = ", "),
                           " clash with Monolix's error-model parameter names; rename them in the model file")
    re <- unlist(config$model$random_effects_excluded)
    if (length(setdiff(re, params))) add("model.random_effects_excluded has unknown parameter(s): ", paste(setdiff(re, params), collapse = ", "))
    dl <- config$model$distributions %||% list()
    if (length(setdiff(names(dl), params))) add("model.distributions has unknown parameter(s): ", paste(setdiff(names(dl), params), collapse = ", "))
    dt <- parse_distributions(dl)$type
    if (any(dt %in% "probitnormal", na.rm = TRUE)) {
      add("model.distributions: probitnormal isn't supported (Monolix 2024R1 leaves the parameter lognormal ",
          "when lixoftConnectors asks for it); use logitnormal, whose default limits are also [0, 1]")
    }
    if (any(is.na(dt)) || length(setdiff(dt, c(names(DISTRIBUTIONS), "probitnormal")))) {
      add("model.distributions values must be lognormal, normal or logitnormal, ",
          "or {distribution: logitnormal, limits: [min, max]}")
    }
    for (k in names(dl)) {
      d <- dl[[k]]
      if (!is.list(d)) next
      extra <- setdiff(names(d), c("distribution", "limits"))
      if (length(extra)) add("model.distributions.", k, ": unknown key(s) ", paste(extra, collapse = ", "), " (use distribution and limits)")
      if (is.null(d$limits)) next
      l <- unlist(d$limits)
      if (!identical(tolower(as.character(d$distribution %||% "")), "logitnormal")) {
        add("model.distributions.", k, ": limits apply only to logitnormal parameters")
      } else if (!is.numeric(l) || length(l) != 2 || !all(is.finite(l)) || l[1] >= l[2]) {
        add("model.distributions.", k, ".limits must be [min, max] with min < max")
      }
    }
    bounds <- parameter_bounds(dl, params)
    for (k in names(config$model$report_transforms)) {
      rt <- config$model$report_transforms[[k]]
      if (!k %in% paste0(params, "_pop")) add("model.report_transforms: '", k, "' is not a population parameter")
      if (is.null(rt$expr) || inherits(tryCatch(eval_transform(rt$expr, 1), error = function(e) e), "error")) {
        add("model.report_transforms.", k, ".expr must be an expression in x, e.g. \"x * 10^5\"")
      }
    }
  }

  # error-model starting values (estimated, not fixed; only the start is set)
  ei <- config$model$error_initial
  if (length(ei)) {
    eiv <- unlist(ei)
    if (!is.numeric(eiv) || any(eiv <= 0)) add("model.error_initial values must be positive numbers")
    if (length(em)) {
      bad <- setdiff(names(eiv), error_param_names(em))
      if (length(bad)) add("model.error_initial has parameter(s) not in the error model (",
                           paste(error_param_names(em), collapse = ", "), "): ", paste(bad, collapse = ", "))
    }
  }

  # estimation
  es <- config$estimation
  for (k in c("burn_in", "exploratory", "smoothing")) {
    if (!is_count(es$saem[[k]])) add("estimation.saem.", k, " must be a positive whole number")
  }
  if (!is_number(es$saem$smoothing_alpha) || es$saem$smoothing_alpha <= 0 || es$saem$smoothing_alpha > 1) {
    add("estimation.saem.smoothing_alpha must be in (0, 1]")
  }
  # scenario
  sc <- config$scenario
  for (k in c(names(SCENARIO_TASKS), "linearization")) {
    if (!is.logical(sc[[k]]) || length(sc[[k]]) != 1) add("scenario.", k, " must be true or false")
  }
  off <- REQUIRED_TASKS[!vapply(REQUIRED_TASKS, function(k) isTRUE(sc[[k]]), logical(1))]
  if (length(off)) {
    add("scenario: the pipeline needs ", paste(off, collapse = ", "),
        " (estimates, individual predictions and BICc to rank runs)")
  }

  # step 1
  s1 <- config$step1_igs
  error_params <- if (length(em)) error_param_names(em) else character()
  ir <- s1$initial_ranges
  if (is.null(ir)) {
    add("REQUIRED: step1_igs.initial_ranges")
  } else if (!is.null(params)) {
    miss <- setdiff(params, names(ir))
    if (length(miss)) add("REQUIRED: step1_igs.initial_ranges for parameter(s): ", paste(miss, collapse = ", "))
    err_keys <- intersect(names(ir), setdiff(error_params, params))
    if (length(err_keys)) {
      add("step1_igs.initial_ranges: error-model parameter(s) ", paste(err_keys, collapse = ", "),
          " aren't sampled; set their starting value in model.error_initial (default: Monolix's)")
    }
    extra <- setdiff(names(ir), c(params, err_keys))
    if (length(extra)) add("step1_igs.initial_ranges has unknown parameter(s): ", paste(extra, collapse = ", "))
    for (k in intersect(names(ir), params)) {
      r <- ir[[k]]
      if (is.null(r$min) || is.null(r$max)) { add("REQUIRED: step1_igs.initial_ranges.", k, " min and max"); next }
      if (!is_number(r$min) || !is_number(r$max) || r$min >= r$max) { add("step1_igs.initial_ranges.", k, ": min must be below max"); next }
      if (isTRUE(r$log10) && r$min <= 0) add("step1_igs.initial_ranges.", k, ": log10 sampling needs min > 0")
      check_bounds(paste0("step1_igs.initial_ranges.", k), paste0(k, "_pop"), c(r$min, r$max))
    }
  }
  if (!is.numeric(s1$omega_range) || length(s1$omega_range) != 2 || s1$omega_range[1] <= 0 || s1$omega_range[1] >= s1$omega_range[2]) {
    add("step1_igs.omega_range must be [min, max] with 0 < min < max")
  }
  if (!is_count(s1$wide_nits)) add("step1_igs.wide_nits must be a positive whole number")
  if (!is_count(s1$narrow$nits)) add("step1_igs.narrow.nits must be a positive whole number")
  if (!is_number(s1$narrow$p_threshold) || s1$narrow$p_threshold <= 0 || s1$narrow$p_threshold > 1) add("step1_igs.narrow.p_threshold must be in (0, 1]")
  if (!is.null(s1$narrow$min_fold_width) && (!is_number(s1$narrow$min_fold_width) || s1$narrow$min_fold_width <= 1)) {
    add("step1_igs.narrow.min_fold_width must be ~ or a number above 1")
  }
  if (!is_number(s1$narrow$omega_initial) || s1$narrow$omega_initial <= 0) add("step1_igs.narrow.omega_initial must be positive")
  for (k in names(s1$profile_likelihood)) {
    pl <- s1$profile_likelihood[[k]]
    pl$range <- unlist(pl$range)  # YAML reads [0.5, 30] (decimal and whole number) as a list
    pl$values <- unlist(pl$values)
    key <- paste0("step1_igs.profile_likelihood.", k)
    if (!is.null(params) && !k %in% c(paste0(params, "_pop"), error_params)) add("step1_igs.profile_likelihood: '", k, "' is not a population parameter")
    if (!is.null(pl$log10) && !(is.logical(pl$log10) && length(pl$log10) == 1 && !is.na(pl$log10))) {
      add(key, ".log10 must be true or false")
    }
    if (!is.null(pl$values)) {
      if (!is.null(pl$range) || !is.null(pl$n)) add(key, ": give either values or range and n, not both")
      if (!is.numeric(pl$values) || length(unique(pl$values)) < 2 || !all(is.finite(pl$values))) {
        add(key, ".values must be a list of at least 2 different numbers")
      } else {
        if (isTRUE(pl$log10) && any(pl$values <= 0)) add(key, ": a log10 axis needs values > 0")
        check_bounds(key, k, pl$values)
      }
      next
    }
    if (!is.numeric(pl$range) || length(pl$range) != 2 || pl$range[1] >= pl$range[2]) {
      add(key, ".range must be [min, max]")
    } else {
      if (isTRUE(pl$log10) && pl$range[1] <= 0) add(key, ": a log10 grid needs range min > 0")
      check_bounds(key, k, pl$range)
    }
    if (!is_count(pl$n) || pl$n < 2) add(key, ".n must be a whole number of at least 2")
  }

  # step 2: fixed parameters (checkpoint 1)
  # (needed before step 2 starts, so `decisions_needed` doesn't apply)
  if (up_to_step >= 2) {
    fx <- config$decisions$fixed
    if (is.null(fx)) {
      todo <- c(todo, paste0("decisions.fixed (checkpoint 1): parameters to fix after step 1, as `name_pop: value`, or {} to fix none ",
                             "(and decisions.random_effects_removed for random effects to drop)"))
    } else if (length(fx)) {
      fxv <- unlist(fx)
      if (!is.numeric(fxv)) add("decisions.fixed values must be numbers")
      if (!is.null(params)) {
        bad <- setdiff(names(fxv), c(paste0(params, "_pop"), error_params))
        if (length(bad)) add("decisions.fixed has unknown parameter(s): ", paste(bad, collapse = ", "), " (use the _pop names)")
        if (is.numeric(fxv)) for (k in names(fxv)) check_bounds(paste0("decisions.fixed.", k), k, fxv[[k]])
      }
    }
    rm <- removed_random_effects(config)
    if (length(rm) && !is.character(rm)) add("decisions.random_effects_removed must be a list of parameter names")
    if (length(rm) && is.character(rm) && !is.null(params)) {
      bad <- setdiff(rm, params)
      if (length(bad)) add("decisions.random_effects_removed has unknown parameter(s): ", paste(bad, collapse = ", "),
                           " (use the parameter names, e.g. k, not omega_k)")
      already <- intersect(rm, unlist(config$model$random_effects_excluded))
      if (length(already)) add("decisions.random_effects_removed: ", paste(already, collapse = ", "),
                               " already ha", if (length(already) == 1) "s" else "ve", " no random effect (model.random_effects_excluded)")
      if (!length(setdiff(params, c(rm, unlist(config$model$random_effects_excluded))))) {
        add("decisions.random_effects_removed leaves no parameter with a random effect; Monolix needs at least one")
      }
    }
    if (!is_count(config$step2_igs$nits)) add("step2_igs.nits must be a positive whole number")
  }

  # step 3: model building and candidates (checkpoint 2)
  if (up_to_step >= 3) {
    # checkpoint 1b: after step 2, confirm nothing else needs fixing (needed before step 3 starts)
    if (step2_runs(config)) {
      fc <- config$decisions$fixed_confirmed
      if (is.null(fc)) {
        todo <- c(todo, paste0("decisions.fixed_confirmed (checkpoint 1b): look at step2_igs/identifiability_summary.csv; ",
                               "true if nothing else needs fixing or removing, otherwise change decisions.fixed or ",
                               "decisions.random_effects_removed and rerun step 2"))
      } else if (!isTRUE(fc)) {
        add("decisions.fixed_confirmed must be true, or ~ until you've checked step 2")
      }
    }
    s3 <- config$step3_model_building
    methods <- tryCatch(model_building_methods(config), error = function(e) { add(conditionMessage(e)); NULL })
    searching <- length(methods) > 0
    if (searching) {
      # search settings matter only when a covariate search runs
      if (!isTRUE(s3$criterion %in% c("BIC", "LRT"))) add("step3_model_building.criterion must be BIC (Monolix's name for BICc) or LRT")
      covs <- unlist(s3$covariates)
      if (!identical(covs, "all") && length(cov_cols) && length(setdiff(covs, cov_cols))) {
        add("step3_model_building.covariates has non-covariate column(s): ", paste(setdiff(covs, cov_cols), collapse = ", "))
      }
      prm <- unlist(s3$parameters)
      if (!identical(prm, "all") && !is.null(params) && length(setdiff(prm, params))) {
        add("step3_model_building.parameters has unknown parameter(s): ", paste(setdiff(prm, params), collapse = ", "))
      }
      check_relationships(s3$relationships, params, cov_cols, prm, covs, add)
      if ("scm_delta2" %in% auto_candidates(config, methods) &&
          (!is_number(s3$scm_stop_rule_delta_bicc) || s3$scm_stop_rule_delta_bicc <= 0)) {
        add("step3_model_building.scm_stop_rule_delta_bicc must be positive")
      }
    }
    if (!is_count(s3$initial_value_rounding)) add("step3_model_building.initial_value_rounding must be a positive whole number")
    cm <- s3$candidate_models
    auto <- if (is.null(methods)) character() else auto_candidates(config, methods)
    if (is.null(cm)) {
      if (up_to_step > 3 || decisions_needed) {
        todo <- c(todo, if (searching) {
          paste0("step3_model_building.candidate_models (checkpoint 2): named treatment-effect models, or {} to compare only ",
                 paste(auto, collapse = " and "))
        } else {
          "step3_model_building.candidate_models: named treatment-effect models for step 4 (method is none, so there are no automatic ones)"
        })
      }
    } else {
      if (!is.null(methods) && !searching && length(cm) < 1) {
        add("step3_model_building.candidate_models needs at least one model when method is none (no covariate search adds any)")
      }
      clash <- intersect(names(cm), AUTO_CANDIDATES)
      if (length(clash)) add("step3_model_building.candidate_models: ", paste(clash, collapse = ", "), " are reserved for the automatic candidates of the covariate searches; pick other names")
      for (nm in names(cm)) {
        te <- cm[[nm]]
        if (!is.null(params) && length(setdiff(names(te), params))) add("candidate '", nm, "' has unknown parameter(s): ", paste(setdiff(names(te), params), collapse = ", "))
        bad_cov <- setdiff(unlist(te), cov_cols)
        if (length(cov_cols) && length(bad_cov)) add("candidate '", nm, "' uses non-covariate column(s): ", paste(bad_cov, collapse = ", "))
      }
    }
  }

  # step 4: comparison and final model (checkpoint 3)
  if (up_to_step >= 4) {
    s4 <- config$step4_model_comparison
    if (!is_count(s4$nits)) add("step4_model_comparison.nits must be a positive whole number")
    candidates <- tryCatch(candidate_names(config), error = function(e) character())
    # checkpoint 2b: compare is decided after step 3's single fits, before step 4 starts
    cmp <- unlist(s4$compare)
    if (is.null(cmp)) {
      todo <- c(todo, paste0("step4_model_comparison.compare (checkpoint 2b): look at step3_model_building/candidates/single_fits.csv; ",
                             "all, or a list of candidates (e.g. [",
                             paste(unique(c(utils::head(candidates, 1), utils::tail(candidates, 1))), collapse = ", "), "])"))
    } else if (!identical(cmp, "all")) {
      if (length(setdiff(cmp, candidates))) add("step4_model_comparison.compare has unknown model(s): ", paste(setdiff(cmp, candidates), collapse = ", "),
                                                " (candidates: ", paste(candidates, collapse = ", "), ")")
    }
    compared <- if (identical(cmp, "all") || is.null(cmp)) candidates else cmp
    if (!isTRUE(s4$fdr_method %in% stats::p.adjust.methods)) add("step4_model_comparison.fdr_method must be one of ", paste(stats::p.adjust.methods, collapse = ", "))
    if (is.null(s4$final_model)) {
      if (decisions_needed) todo <- c(todo, "step4_model_comparison.final_model (checkpoint 3): the candidate to keep as the final model")
    } else if (!s4$final_model %in% compared) {
      add("step4_model_comparison.final_model '", s4$final_model, "' is not one of the compared models")
    }
  }

  # report
  sb <- unlist(config$report$stratify_by)
  if (length(sb) && length(cov_cols) && length(setdiff(sb, cov_cols))) {
    add("report.stratify_by has non-covariate column(s): ", paste(setdiff(sb, cov_cols), collapse = ", "))
  }

  if (length(p)) {
    stop("analysis.yaml has ", length(p), " problem(s):\n", paste0("  - ", p, collapse = "\n"),
         if (length(todo)) paste0("\nAlso still to decide:\n", paste0("  - ", todo, collapse = "\n")) else "",
         call. = FALSE)
  }
  if (length(todo)) {
    msg <- paste0("Checkpoint: record your decision in ", attr(config, "path") %||% "analysis.yaml", ", then rerun:\n",
                  paste0("  - ", todo, collapse = "\n"))
    stop(structure(class = c("pipeline_checkpoint", "error", "condition"),
                   list(message = msg, call = NULL, todo = todo)))
  }
  invisible(config)
}

eval_transform <- function(expr, x) {
  f <- parse(text = expr)[[1]]
  vars <- all.vars(f)
  if (!identical(vars, "x") && length(vars) > 0) stop("only x allowed")
  eval(f, list(x = x), baseenv())
}

`%||%` <- function(a, b) if (is.null(a)) b else a
# step3_model_building.relationships: `include` and `exclude`, each `parameter: [covariates]`,
# within the parameters and covariates the search uses. `add` records a problem.
check_relationships <- function(rel, params, cov_cols, prm, covs, add) {
  if (!length(rel)) return(invisible())
  key <- "step3_model_building.relationships"
  if (!is.list(rel) || is.null(names(rel)) || length(setdiff(names(rel), c("include", "exclude")))) {
    add(key, " takes `include` and/or `exclude`, each as `parameter: [covariates]`, e.g. {include: {param1: [cov1]}}")
    return(invisible())
  }
  all_pairs <- list()
  for (k in names(rel)) {
    x <- rel[[k]]
    if (!length(x)) next
    if (!is.list(x) || is.null(names(x))) { add(key, ".", k, " must be `parameter: [covariates]`, e.g. {param1: [cov1]}"); next }
    p_bad <- if (!is.null(params)) setdiff(names(x), params) else character()
    if (length(p_bad)) add(key, ".", k, " has unknown parameter(s): ", paste(p_bad, collapse = ", "))
    p_out <- if (!identical(prm, "all")) setdiff(intersect(names(x), params %||% names(x)), prm) else character()
    if (length(p_out)) add(key, ".", k, " has parameter(s) not in step3_model_building.parameters: ", paste(p_out, collapse = ", "))
    c_bad <- if (length(cov_cols)) setdiff(unlist(x), cov_cols) else character()
    if (length(c_bad)) add(key, ".", k, " has non-covariate column(s): ", paste(c_bad, collapse = ", "))
    c_out <- if (!identical(covs, "all")) setdiff(intersect(unlist(x), if (length(cov_cols)) cov_cols else unlist(x)), covs) else character()
    if (length(c_out)) add(key, ".", k, " has covariate(s) not in step3_model_building.covariates: ", paste(c_out, collapse = ", "))
    all_pairs[[k]] <- paste(rep(names(x), lengths(x)), unlist(x, use.names = FALSE))
  }
  both <- intersect(all_pairs$include, all_pairs$exclude)
  if (length(both)) add(key, ": ", paste(sub(" ", "-", both), collapse = ", "), " is in both include and exclude")
  invisible()
}

is_blank <- function(x) is.null(x) || !nzchar(as.character(x)[1])
is_number <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
is_count <- function(x, zero_ok = FALSE) is_number(x) && x == round(x) && (x > 0 || (zero_ok && x == 0))
