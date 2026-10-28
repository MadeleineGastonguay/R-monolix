# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Bug fixes since porting: fixed values matched by name; results read from the project's own
# folder (no globals); attempt cap and seed in initial_guess_sampling(); only missing betas
# set to 0 in set_base(); fixed Monolix seed; BICc as default model-building criterion;
# create_project() stops when newProject() fails; generate_report() takes groups, labels and
# rescaling from the config (template generalized).

#' Set scenario tasks and SAEM settings for the loaded project
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
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
#' Ported from the SIV_HIV_rebound_kinetics functions.R. With a config, use
#' [setup_project()] instead, which takes everything from a [model_spec()].
#'
#' @param distributions Named character vector or list of parameter
#'   distributions, e.g. `c(ta = "normal")`.
#' @param estimation,scenario Passed to [set_project_settings()].
#' @export
create_project <- function(model_file, data_file, 
                           headerTypes, error_model, 
                           exclude_random_effects, treatment_effects = NULL, 
                           normal_dist = NULL, fixed = NULL,
                           linearize = FALSE, distributions = NULL, estimation = NULL, scenario = NULL){
  
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
  
  # `distributions` (named, e.g. c(ta = "normal")) generalizes `normal_dist`
  if(!is.null(normal_dist) || length(distributions) > 0){
    current <- lixoftConnectors::getIndividualParameterModel()$distribution
    current[normal_dist] <- "normal"
    current[names(distributions)] <- unlist(distributions)
    lixoftConnectors::setIndividualParameterDistribution(current)
  }
  
  if(!is.null(fixed)){
    popparams <- lixoftConnectors::getPopulationParameterInformation()
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- unlist(fixed)[popparams$name[popparams$name %in% names(fixed)]]
    lixoftConnectors::setPopulationParameterInformation(popparams)
  }
  
  
  print(lixoftConnectors::getIndividualParameterModel())
  
  set_project_settings(linearize, estimation, scenario)
  
}

#' Set initial values (and fixed parameters) on the loaded project
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
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
  popparams$initialValue <- initial_values[popparams$name, "estimate"]
  
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
#' Ported from the SIV_HIV_rebound_kinetics functions.R and generalized: the
#' groups, labels and rescaling come from the config instead of the SIV study
#' columns. Needs Monolix ([init_monolix()]) and a project that has been run.
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

