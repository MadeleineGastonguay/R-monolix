# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Bug fixes since porting: fixed values matched by name; results read from the project's own
# folder (no globals); attempt cap and seed in initial_guess_sampling(); only missing betas
# set to 0 in set_base(); fixed Monolix seed; BICc as default model-building criterion.

#' Set scenario tasks and SAEM settings for the loaded project
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
set_project_settings <- function(linearize = FALSE){
  # Set scenario for Monolix
  scenario <- lixoftConnectors::getScenario()
  scenario$tasks["standardErrorEstimation"] <- TRUE
  scenario$tasks["logLikelihoodEstimation"] <- TRUE
  scenario$linearization <- linearize
  lixoftConnectors::setScenario(scenario)
  
  # Set Parameter Estimation Settings
  # One fixed Monolix seed so SAEM runs are reproducible
  lixoftConnectors::setProjectSettings(seed = MONOLIX_SEED)
  
  lixoftConnectors::setPopulationParameterEstimationSettings(
    nbBurningIterations=50,
    nbExploratoryIterations=1000,
    nbSmoothingIterations=500,
    smoothingAlpha=0.7,
    exploratoryAutoStop = TRUE,
    smoothingAutoStop = TRUE
  )
  
}

#' Create a new Monolix project from a model file, data file and statistical model
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
create_project <- function(model_file, data_file, 
                           headerTypes, error_model, 
                           exclude_random_effects, treatment_effects = NULL, 
                           normal_dist = NULL, fixed = NULL,
                           linearize = FALSE){
  
  # Create a new project
  lixoftConnectors::newProject(modelFile = model_file,
             data = list(dataFile = data_file,
                         headerTypes = headerTypes
             )
  )
  
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
  
  if(!is.null(normal_dist)){
    distributions <- lixoftConnectors::getIndividualParameterModel()$distribution
    distributions[normal_dist] <- "normal"
    lixoftConnectors::setIndividualParameterDistribution(distributions)
  }
  
  if(!is.null(fixed)){
    popparams <- lixoftConnectors::getPopulationParameterInformation()
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- unlist(fixed)[popparams$name[popparams$name %in% names(fixed)]]
    lixoftConnectors::setPopulationParameterInformation(popparams)
  }
  
  
  print(lixoftConnectors::getIndividualParameterModel())
  
  set_project_settings(linearize)
  
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
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
generate_report <- function(project_file, output_file = "report.html"){
  output_file <- file.path(str_remove(project_file, "\\.mlxtran$"), output_file)
  rmarkdown::render(
    system.file("rmarkdown", "templates", "monolix_report", "skeleton", "skeleton.Rmd", package = "monolixPipelineR"),
    params = list(project = project_file),
    output_file = output_file
  )
}

