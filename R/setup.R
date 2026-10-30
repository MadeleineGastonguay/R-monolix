# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Bug fixes since porting: fixed values matched by name; results read from the project's own
# folder (no globals); attempt cap and seed in initial_guess_sampling(); only missing betas
# set to 0 in set_base(); fixed Monolix seed; BICc as default model-building criterion;
# create_project() stops when newProject() fails and sets logit-normal limits; generate_report()
# takes groups, labels and rescaling from the config (template generalized).

#' Set scenario tasks and SAEM settings for the loaded project
#'
#' @param linearize Used only when `scenario` is `NULL`.
#' @param estimation The `estimation` section of the config (see
#'   [config_defaults()]); `NULL` uses the defaults.
#' @param scenario The `scenario` section of the config; `NULL` uses the
#'   default tasks with `linearization = linearize`.
#' @export
set_project_settings <- function(linearize = FALSE, estimation = NULL, scenario = NULL){
  require_monolix()
  if(is.null(estimation)) estimation <- config_defaults()$estimation
  if(is.null(scenario)){
    scenario <- config_defaults()$scenario
    scenario$linearization <- linearize
  }
  saem <- estimation$saem

  # Set scenario for Monolix
  apply_scenario(scenario)
  
  # One fixed Monolix seed so SAEM runs are reproducible
  lixoftConnectors::setProjectSettings(seed = MONOLIX_SEED)
  
  # Set Parameter Estimation Settings
  lixoftConnectors::setPopulationParameterEstimationSettings(
    nbBurningIterations = saem$burn_in,
    nbExploratoryIterations = saem$exploratory,
    nbSmoothingIterations = saem$smoothing,
    smoothingAlpha = saem$smoothing_alpha,
    exploratoryAutoStop = isTRUE(saem$exploratory_autostop),
    smoothingAutoStop = isTRUE(saem$smoothing_autostop)
  )
  
}

#' Create a new Monolix project from a model file, data file and statistical model
#'
#' With a config, use [setup_project()] instead, which takes everything from a [model_spec()].
#'
#' @param distributions Named character vector or list of parameter
#'   distributions, e.g. `c(param1 = "normal", param2 = "logitnormal")`
#'   (`lognormal`, `normal`, `logitnormal`, or Monolix's spelling, e.g.
#'   `logNormal`). Stops if Monolix doesn't apply one of them.
#' @param limits Named list of `c(min, max)` bounds for logit-normal
#'   parameters, e.g. `list(param2 = c(0, 100))`; Monolix's default is (0, 1).
#' @param estimation,scenario Passed to [set_project_settings()].
#' @param model_file Path to the Mlxtran structural model.
#' @param data_file Path to the data file.
#' @param headerTypes Monolix column type of each data column, in data order
#'   (e.g. `"id"`, `"time"`, `"observation"`, `"catcov"`).
#' @param error_model Error model of each observation, e.g.
#'   `list(y = "combined1")`.
#' @param exclude_random_effects Parameters without a random effect.
#' @param treatment_effects Named list of covariates per parameter, e.g.
#'   `list(param1 = "cov1")`; `NULL` for none.
#' @param normal_dist Parameters with a normal distribution (older form of
#'   `distributions`).
#' @param fixed Named values of parameters to fix, e.g. `c(param1_pop = 10)`.
#' @param linearize Use linearization for standard errors and the likelihood.
#' @return Called for its effect: the new project is loaded in Monolix.
#' @export
create_project <- function(model_file, data_file, 
                           headerTypes, error_model, 
                           exclude_random_effects, treatment_effects = NULL, 
                           normal_dist = NULL, fixed = NULL,
                           linearize = FALSE, distributions = NULL, limits = NULL,
                           estimation = NULL, scenario = NULL){
  
  # Create a new project. newProject() reports failures as a message and returns FALSE
  # instead of raising an error, so check its value.
  created <- lixoftConnectors::newProject(modelFile = model_file,
             data = list(dataFile = data_file,
                         headerTypes = headerTypes
             )
  )
  if (isFALSE(created)) {
    stop("Monolix could not create the project from ", model_file,
         "; see the [ERROR] message above.", call. = FALSE)
  }
  
  # Set Error model
  lixoftConnectors::setErrorModel(error_model)
  
  # Remove random effects if needed
  if(!is.null(exclude_random_effects)){
    var_model <- lixoftConnectors::getIndividualParameterModel()$variability$id
    var_model[exclude_random_effects] <- FALSE
    lixoftConnectors::setIndividualParameterVariability(var_model)
  }
  
  if(!is.null(treatment_effects)){
    # treatment_effects should be a list of vectors of categorical variables to add to parameters
    covariateModel <- lixoftConnectors::getIndividualParameterModel()$covariateModel
    for(parameter in names(treatment_effects)){
      covariateModel[[parameter]][treatment_effects[[parameter]]] <- TRUE
    }
    lixoftConnectors::setCovariateModel(covariateModel)
  }
  
  # Monolix sets logit limits only when the starting value is already inside them, so move it
  # there first (before the distribution too, to avoid Monolix adjusting the default limits)
  if(length(limits)) start_inside_limits(limits)

  # `distributions` (named, e.g. c(ta = "normal")) generalizes `normal_dist`
  if(!is.null(normal_dist) || length(distributions) > 0){
    current <- lixoftConnectors::getIndividualParameterModel()$distribution
    current[normal_dist] <- "normal"
    current[names(distributions)] <- monolix_distribution(unlist(distributions))
    lixoftConnectors::setIndividualParameterDistribution(current)
    # Monolix can ignore a distribution without an R error (probitNormal in 2024R1): check
    got <- lixoftConnectors::getIndividualParameterModel()$distribution
    got <- as.character(got[names(current)])
    bad <- names(current)[is.na(got) | tolower(got) != tolower(current)]
    if (length(bad)) {
      stop("Monolix did not apply the distribution of ", paste0(bad, " (", current[bad], ")", collapse = ", "),
           "; it kept ", paste(got[match(bad, names(current))], collapse = ", "), ". See any [ERROR] message above.", call. = FALSE)
    }
  }
  if(length(limits)) set_logit_limits(limits)
  
  if(!is.null(fixed)){
    popparams <- lixoftConnectors::getPopulationParameterInformation()
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- unlist(fixed)[popparams$name[popparams$name %in% names(fixed)]]
    lixoftConnectors::setPopulationParameterInformation(popparams)
  }
  
  
  print(lixoftConnectors::getIndividualParameterModel())
  
  set_project_settings(linearize, estimation, scenario)
  
}

# Monolix's spelling of a distribution name (lognormal -> logNormal); other names unchanged
monolix_distribution <- function(x) {
  out <- unname(DISTRIBUTIONS[tolower(x)])
  ifelse(is.na(out), x, out)
}

# Move the starting value of each limited parameter (`<p>_pop`) to the middle of its limits
# when it's outside them (Monolix's default is 1). The IGS sets its own starting values later.
start_inside_limits <- function(limits) {
  pop <- lixoftConnectors::getPopulationParameterInformation()
  moved <- FALSE
  for (p in names(limits)) {
    l <- as.numeric(limits[[p]])
    i <- pop$name == paste0(p, "_pop")
    if (any(i) && outside_bounds(pop$initialValue[i], l[1], l[2])) {
      pop$initialValue[i] <- mean(l)
      moved <- TRUE
    }
  }
  if (moved) lixoftConnectors::setPopulationParameterInformation(pop)
  invisible(moved)
}

# Set the bounds of logit-normal parameters on the loaded project (their starting values must
# already be inside them: start_inside_limits()), and check that Monolix took them.
set_logit_limits <- function(limits) {
  limits <- lapply(limits, as.numeric)
  do.call(lixoftConnectors::setIndividualLogitLimits, limits)
  got <- lixoftConnectors::getIndividualParameterModel()$limits
  ok <- vapply(names(limits), function(p) isTRUE(all.equal(as.numeric(unlist(got[[p]])), limits[[p]])), logical(1))
  if (!all(ok)) {
    stop("Monolix did not set the limits of ", paste(names(limits)[!ok], collapse = ", "),
         "; see the [ERROR] message above. Limits apply only to logitnormal parameters.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Set initial values (and fixed parameters) on the loaded project
#' @param initial_values Data frame with `parameter` and `estimate`: the
#'   starting values. Parameters missing from it keep the project's values;
#'   treatment effects (`beta_`) missing from it start at 0.
#' @param omega_under_1 Allow random-effect SDs (`omega_`) to start below 1?
#'   By default they start at 1 or more, which helps SAEM converge.
#' @param fixed Named values of parameters to fix, e.g. `c(param1_pop = 10)`.
#' @return Called for its effect on the loaded project.
#' @export
set_base <- function(initial_values, omega_under_1 = F, fixed = NULL){
  
  # Usually, we don't want an initial guess for omega parameters less than 1 to help with convergence of SAEM
  if(!omega_under_1){
    initial_values <- initial_values %>% 
      rowwise %>% 
      # don't use an initial condition less than 1 for random effects to help with convergence
      mutate(estimate = ifelse(grepl("omega", parameter), max(1, estimate), estimate)) 
  }
  
  initial_values <- initial_values %>% 
    column_to_rownames("parameter") 
  
  popparams <- lixoftConnectors::getPopulationParameterInformation()
  # Parameters without a value (e.g. the fixed exponent c of a proportional error model, which
  # has no estimate) keep the project's own value: an NA here becomes `value=NaN` in the
  # project, and Monolix then can't parse the model-building files it writes from it.
  given <- initial_values[popparams$name, "estimate"]
  keep <- is.na(given) & !grepl("beta_", popparams$name)
  popparams$initialValue[!keep] <- given[!keep]
  
  # Set values for fixed parameters
  if(!is.null(fixed)){
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- unlist(fixed)[popparams$name[popparams$name %in% names(fixed)]]
  }
  
  # If an initial condition for a treatment effect is not provided, set it to 0
  # (only the missing ones; supplied treatment effects are kept)
  missing_beta <- grepl("beta_", popparams$name) & is.na(popparams$initialValue)
  popparams[missing_beta, "initialValue"] <- 0
  
  lixoftConnectors::setPopulationParameterInformation(popparams)
}

#' Render the HTML report for a Monolix project
#'
#' The groups, labels and rescaling come from the config. Needs Monolix ([init_monolix()]) and a project that has been run.
#'
#' @param project_file Path to the .mlxtran project.
#' @param output_file File name of the report, written to the project's
#'   results folder.
#' @param config Optional analysis config ([read_analysis_config()]); supplies
#'   the defaults for the arguments below.
#' @param stratify_by Covariate columns that define the report groups
#'   (`report.stratify_by`). `NULL`: no groups.
#' @param group_labels Group key -> label (`report.group_labels`), see
#'   [report_groups()].
#' @param report_transforms Rescaling of population estimates
#'   (`model.report_transforms`), see [transform_estimates()].
#' @param title Report title.
#' @param simulate_variables If `TRUE`, simulate the model variables on a
#'   fine grid with Simulx first ([simulate_model_variables()]). If that
#'   fails, the report plots them at the observation times (tables.txt).
#' @param variables Passed to [simulate_model_variables()].
#' @return Path to the report, invisibly.
#' @export
generate_report <- function(project_file, output_file = "report.html", config = NULL,
                            stratify_by = config$report$stratify_by,
                            group_labels = config$report$group_labels,
                            report_transforms = config$model$report_transforms,
                            title = "Monolix report", simulate_variables = TRUE, variables = NULL){
  require_monolix()
  project_file <- normalizePath(project_file, mustWork = TRUE)
  # The template sets ggplot's global theme for Monolix's plots; put the caller's back after
  old_theme <- ggplot2::theme_get()
  on.exit(ggplot2::theme_set(old_theme), add = TRUE)
  # Fine-grained model-variable curves; stale curves are removed so a failed simulation
  # doesn't leave curves from an earlier fit
  sim_file <- simulated_variables_file(str_remove(project_file, "\\.mlxtran$"))
  if (file.exists(sim_file)) file.remove(sim_file)
  if (simulate_variables) {
    tryCatch(simulate_model_variables(project_file, variables),
             error = function(e) warning("Simulx simulation of the model variables failed (the report ",
                                         "uses tables.txt instead): ", conditionMessage(e), call. = FALSE))
  }
  output_file <- file.path(str_remove(project_file, "\\.mlxtran$"), output_file)
  template <- system.file("rmarkdown", "templates", "monolix_report", "skeleton", "skeleton.Rmd",
                          package = "monolixPipelineR")
  out <- rmarkdown::render(
    template,
    params = list(project = project_file, title = title, stratify_by = unlist(stratify_by),
                  group_labels = group_labels, report_transforms = report_transforms),
    output_file = output_file,
    # knit in a temporary folder, so intermediate files don't land in the package
    intermediates_dir = tempdir(),
    # the template uses package helpers and imports (dplyr, ggplot2, ...)
    envir = new.env(parent = asNamespace("monolixPipelineR")),
    quiet = TRUE
  )
  invisible(out)
}

