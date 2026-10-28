# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Bug fixes since porting: fixed values matched by name; results read from the project's own
# folder (no globals); seed and predefined-run handling in initial_guess_sampling();
# results read with the scenario's method (linearization or not); final_model.mlxtran always saved (before any rerun, so the ModelBuilding results are kept); only missing betas
# set to 0 in set_base(); fixed Monolix seed; SCM as default model-building method;
# condition number is NA instead of an error when the correlation matrix has nan;
# each optional record step has its own error handler; sampling moved to draw_initial_values().

#' Fit the loaded project from many sampled initial values
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#'
#' @param predefined_initialvalues Optional table of initial values (columns
#'   `parameter`, `run`, `initial_value`), e.g. shared across candidate models in
#'   step 4. When given, each predefined set is fitted exactly once, under its own
#'   run number, and `nits` is ignored. Failed runs are not replaced, so paired
#'   comparisons only use the runs that succeeded in every model.
#' @param set_error Named starting values for error-model parameters (see
#'   [error_starts()]). They are still estimated; only the start is set. If
#'   `NULL`, error parameters follow `set_variance` or their sampling range, as
#'   in the SIV code.
#' @param seed Seed for sampling initial values (`NULL` leaves the RNG alone).
#' @export
initial_guess_sampling <- function(initialValue_ranges, log10 = NULL,
                                   nits = 100, set_variance = NULL, set_treatment = NULL, set_error = NULL,
                                   record_shrinkage = FALSE,
                                   fixed = NULL, predefined_initialvalues = NULL,
                                   record_pvals = FALSE, record_indiv_params = FALSE,
                                   record_condition_number = FALSE, record_table = FALSE,
                                   seed = NULL){

  require_monolix()
  if(!is.null(seed)) set.seed(seed)

  popparams <- lixoftConnectors::getPopulationParameterInformation()
  methods <- result_methods()  # linearization or not, from the project's scenario

  estimates <- tabestimates <- tabiters <-  phase <- initialvalues <- shrinkage_df <- treatment_pvals <- indivParams <- predictions <- condition_numbers <- table_predictions <- NULL;

  # With predefined initial values, fit each predefined set once and only those
  if(!is.null(predefined_initialvalues)){
    run_ids <- sort(unique(predefined_initialvalues$run))
    if(length(run_ids) != nits){
      warning("nits = ", nits, " but ", length(run_ids), " predefined initial value sets were given; ",
              "using only the ", length(run_ids), " predefined sets.", call. = FALSE)
    }
  }
  
  run_condition <- TRUE
  k <- 1
  i <- if(is.null(predefined_initialvalues)) 1 else run_ids[1]
  while(run_condition){
    print(i)

    if(is.null(predefined_initialvalues)){
      
      # sample new initial estimates (see draw_initial_values() in sampling.R)
      popini <- draw_initial_values(popparams$name, initialValue_ranges, log10 = log10, fixed = fixed,
                                    set_error = set_error, set_variance = set_variance,
                                    set_treatment = set_treatment)
      
      # set sampled values as new initial estimates
      newpopparams <- popparams
      newpopparams$initialValue <- popini
      initialvalues <- rbind(initialvalues, newpopparams %>% mutate(run = i) %>% rename(parameter = name, initial_value = initialValue))
      
      if(!is.null(fixed)){
        newpopparams[newpopparams$name %in% names(fixed),"method"] <- "FIXED"
        newpopparams[newpopparams$name %in% names(fixed),"initialValue"] <- unlist(fixed)[newpopparams$name[newpopparams$name %in% names(fixed)]]
      }
      
    }else{
      
      igs <- predefined_initialvalues %>% filter(run == i) %>% as.data.frame()
      newpopparams <- popparams %>% rowwise %>%  
        # Only set initial conditions for parameters that have saved initial conditions. Otherwise, use default values
        mutate(initialValue = ifelse(name %in% igs$parameter, igs[igs$parameter == name, "initial_value"], initialValue)) %>% 
        as.data.frame() 
      initialvalues <- rbind(initialvalues, newpopparams %>% mutate(run = i) %>% rename(parameter = name, initial_value = initialValue))
    }
    
    lixoftConnectors::setPopulationParameterInformation(newpopparams)
    
    # run the estimation
    tryCatch({
      lixoftConnectors::runScenario()
      
      # store the estimates and s.e. in a table
      estimates <- as.data.frame(lixoftConnectors::getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      # standard errors and likelihood from the method the scenario uses
      rses <- data.frame(rse = read_rses(methods, rownames(estimates)), row.names = rownames(estimates))
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$run <- i
      ll <- read_likelihood(methods)
      estimates$ofv <- ll[["ofv"]]
      estimates$aic <- ll[["aic"]]
      estimates$bic <- ll[["bic"]]
      estimates$bicc <- ll[["bicc"]]
      names(estimates)[names(estimates) == "Row.names"] <- "parameter"
      tabestimates <- rbind(tabestimates, estimates)
      
      # store the iterations
      iters <- lixoftConnectors::getSAEMiterations()$estimates %>% mutate(iteration = 1:nrow(.))
      iters$run <- i
      tabiters <- rbind(tabiters, iters)
      
      # Optional records. Each has its own error handler, so a failure in one (e.g. the
      # condition number) doesn't lose the others for this run; the run itself is kept.
      if(record_pvals) record_step("Wald tests", i, {
        p <- lixoftConnectors::getTests()$wald
        if(!is.null(p)) treatment_pvals <- rbind(treatment_pvals, p %>% mutate(run = i))
      })

      # Store Phase Iteration
      phase <- rbind(phase, lixoftConnectors::getSAEMiterations()$iterationNumbers)

      if(record_indiv_params) record_step("individual parameters", i, {
        indivParams <- rbind(
          indivParams,
          lixoftConnectors::getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>%
            mutate(run = i)
        )
        predictions <- rbind(
          predictions,
          read_delim(file.path(project_results_dir(), "predictions.txt"), show_col_types = FALSE) %>%
            mutate(run = i)
        )
      })

      if(record_shrinkage) record_step("shrinkage", i, {
        shrinkage_df <- rbind(shrinkage_df, read_shrinkage(i))
      })

      if(record_condition_number) record_step("condition number", i, {
        condition_number <- condition_number_of(lixoftConnectors::getCorrelationOfEstimates()[[methods$se]], run = i)
        condition_numbers <- rbind(condition_numbers, data.frame(run = i, condition_number = condition_number))
      })

      if(record_table) record_step("tables.txt", i, {
        table_predictions <- rbind(
          table_predictions,
          read_delim(file.path(project_results_dir(), "tables.txt"), show_col_types = FALSE) %>%
            mutate(run = i)
        )
      })

    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
    if(is.null(predefined_initialvalues)){
      # Random sampling: keep drawing until nits runs succeed
      i <- i + 1
      if(length(unique(tabestimates$run)) >= nits) run_condition <- FALSE
    }else{
      # Predefined values: stop after the last predefined set
      k <- k + 1
      if(k > length(run_ids)) run_condition <- FALSE else i <- run_ids[k]
    }
    
  }
  
  if(is.null(tabestimates)) stop("All runs failed; see the errors printed above.", call. = FALSE)
  
  if(!is.null(predefined_initialvalues)){
    failed <- setdiff(run_ids, unique(tabestimates$run))
    if(length(failed) > 0){
      warning(length(failed), " of ", length(run_ids), " predefined runs failed (runs ",
              paste(failed, collapse = ", "), "). They are not replaced.", call. = FALSE)
    }
  }
  
  out <- list(tabestimates = tabestimates,
              tabiters = tabiters,
              phase = phase,
              initialvalues = initialvalues)
  
  if(record_shrinkage){
    out$shrinkage <- shrinkage_df
  }
  
  if(record_pvals){
    out$pvals <- treatment_pvals
  }
  
  if(record_indiv_params){
    out$indivParams <- indivParams
    out$indivFit <- predictions
  }
  
  if(record_condition_number){
    out$condition_number <- condition_numbers
  }
  
  if(record_table){
    out$table_predictions <- table_predictions
  }
  
  # Rerun the model with the initial values that gave the best fit:
  best_run <- tabestimates %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique
  igs <- initialvalues %>% filter(run == best_run) %>% rename(name = parameter, initialValue = initial_value) %>% select(-run)
  if(!is.null(fixed)){
    igs[igs$name %in% names(fixed),"method"] <- "FIXED"
    igs[igs$name %in% names(fixed),"initialValue"] <- unlist(fixed)[igs$name[igs$name %in% names(fixed)]]
  }
  lixoftConnectors::setPopulationParameterInformation(igs)
  lixoftConnectors::runScenario()
  
  return(out)
}

# Run one optional record step of initial_guess_sampling(); on error print it and carry on.
# Assignments inside `expr` (with <-) are made in the caller's frame.
record_step <- function(what, run, expr) {
  env <- parent.frame()
  tryCatch(eval(substitute(expr), env),
           error = function(e) cat("ERROR recording ", what, " for run ", run, ": ", conditionMessage(e), "\n", sep = ""))
  invisible(NULL)
}

#' Profile likelihood for one population parameter
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#'
#' @param linearize `NULL` (default) uses the project's scenario; `TRUE` or
#'   `FALSE` overrides it.
#' @export
profile_likelihood <- function(parameter, range, igs, nits, use_estimates_as_ig = F, linearize = NULL, record_predictions = F){
  
  values <- seq(range[1], range[2], length.out = nits)
  
  initialValues <- igs
  if("estimate" %in% names(igs)){
    initialValues <- igs %>% select(name = parameter, initialValue = estimate) %>%
      mutate(method = "MLE", initialValue = ifelse(grepl("pop", name), initialValue, pmax(initialValue, 1))) 
  }
  if("parameter" %in% names(igs)){
    names(initialValues)[names(initialValues) == "parameter"] <- "name"
  }
  if("initial_value" %in% names(initialValues)){
    names(initialValues)[names(initialValues) == "initial_value"] <- "initialValue"
  }
  
  
  initialValues[initialValues$name == parameter, "method"] <- "FIXED"
  
  # linearize = NULL keeps the project's scenario; TRUE/FALSE overrides it
  if(!is.null(linearize)){
    scenario <- lixoftConnectors::getScenario()
    scenario$linearization <- linearize
    lixoftConnectors::setScenario(scenario)
  }
  methods <- result_methods()
  
  tabestimates <- NULL
  predictions <- NULL
  for(i in 1:nits){
    print(i)
    
    # Set initial conditions
    value <- values[i]
    
    initialValues[initialValues$name == parameter, "initialValue"] <- value
    lixoftConnectors::setPopulationParameterInformation(initialValues)
    
    # print(lixoftConnectors::getPopulationParameterInformation())
    
    tryCatch({
      # Model fitting
      lixoftConnectors::runScenario()
      
      # Get standard errors and likelihoods:
      likelihoods <- read_likelihood(methods)

      # store the estimates and s.e. in a table
      estimates <- as.data.frame(lixoftConnectors::getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      rses <- data.frame(rse = read_rses(methods, rownames(estimates)), row.names = rownames(estimates))
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$profiled_parameter <- parameter
      estimates$profiled_value <- value
      estimates$ofv <- likelihoods[1]
      estimates$aic <- likelihoods[2]
      estimates$bic <- likelihoods[3]
      estimates$bicc <- likelihoods[4]
      names(estimates)[names(estimates) == "Row.names"] <- "parameter"
      
      print(paste0(parameter, ": ", value, " bicc: ", round(unique(estimates$bicc))))
      
      # Combine results
      tabestimates <- rbind(tabestimates, estimates)
      
      if(record_predictions){
        vl_pred <- lixoftConnectors::getEstimatedPredictions()$tot  %>% mutate(id = as.character(id)) 
        EBE <- lixoftConnectors::getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode  %>% mutate(id = as.character(id))
        table_pred <- read_delim(file.path(project_results_dir(), "tables.txt")) %>% 
          select(id, time, contains("mode"))  %>% mutate(id = as.character(id))
        predictions_temp <- right_join(
          EBE,
          left_join(
            vl_pred, table_pred
          )
        ) %>% mutate(profiled_parameter = parameter, profiled_value = value) %>% 
          select(profiled_parameter, profiled_value, everything())
        
        # Combine resuts
        predictions <- rbind(predictions, predictions_temp)
      }
      
      
      if(use_estimates_as_ig){
        initialValues <- estimates %>% select(name = parameter, initialValue = estimate) %>% 
          mutate(method = "MLE", initialValue = ifelse(grepl("pop", name), initialValue, pmax(initialValue, 1))) 
        initialValues <- rbind(initialValues, data.frame(name = parameter, initialValue = value, method = "FIXED"))
      }
      
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
  }
  
  if(record_predictions){
    return(list(tabestimates = tabestimates, predictions = predictions))
  }else{
    return(tabestimates)
  }
  
}

#' Run Monolix covariate model building (SCM or COSSAC)
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R. Results go to
#' `<project>/ModelBuilding` and `<project>/final_model.mlxtran`, so to run
#' both SCM and COSSAC, give each its own `projectFile` (e.g. `scm.mlxtran`
#' and `cossac.mlxtran`).
#'
#' @param scenario Optional `scenario` section of the config, used when the
#'   best model's tasks have to be rerun. `NULL` runs every task without
#'   linearization, as the SIV code did.
#' @export
covariate_model_building <- function(projectFile, method = "SCM", criteria = "BIC", covariates = "all", 
                                     linearization = FALSE, relationships = NULL, parameters = NULL, lrt_thresold = c(0.05,0.01),
                                     scenario = NULL){
  
  if(dir.exists(file.path(gsub(".mlxtran","",projectFile), "ModelBuilding"))){ 
    stop("Model building directory already exists. Aborting to avoid overwriting previous results.")
  }else if (!file.exists(projectFile)){
    lixoftConnectors::saveProject(projectFile)
  }
  
  # get the default covariate search settings ####
  set <- lixoftConnectors::getModelBuildingSettings()
  
  # modify the settings
  set$strategy <- tolower(method) # Monolix's names: "scm" or "cossac"
  set$useSambaBeforeCossac <- FALSE # if true, corresponds to "covSAMBA-COSSAC" in the GUI
  set$criterion <- criteria # "BIC" (Monolix uses this for BICc) or "LRT"
  set$useLin <- linearization
  set$threshold$lrt <- lrt_thresold # first element is for forward selection, second for backward elimination
  
  if(!is.null(parameters)){
    set$parameters <- parameters
  }
  
  if(!is.null(relationships)){
    set$relationships <- relationships
  }
  
  if(!any(covariates == "all")){
    set$covariates <- covariates
  }
  
  #####
  # Run search
  lixoftConnectors::runModelBuilding(settings=set)
  model_building_results <- lixoftConnectors::getModelBuildingResults()
  
  #####
  # Generate report for best model
  # Get best model
  best_model <- which(unlist(model_building_results %>% lapply("[[", "bestModel")))
  best_model_file <- sprintf("ModelBuilding_%03d", best_model)
  best_cov <- model_building_results[[best_model]]$individualModels
  
  final_projectFile <- file.path(gsub(".mlxtran", "", projectFile), "ModelBuilding", sprintf("ModelBuilding_%03d", best_model), paste0(sprintf("ModelBuilding_%03d", best_model), ".mlxtran") )
  lixoftConnectors::loadProject(final_projectFile)
  # Save as final_model before any rerun: a run writes results to the loaded project's folder,
  # which would overwrite the ModelBuilding_00x results with a different scenario.
  final_model_file <- file.path(gsub(".mlxtran", "", projectFile), "final_model.mlxtran")
  lixoftConnectors::saveProject(final_model_file)

  if(any(as.character(unlist(lixoftConnectors::getLaunchedTasks())) == "FALSE")){
    # rerun with the analysis scenario (or, without one, every task as in the SIV code)
    if(is.null(scenario)){
      sc <- lixoftConnectors::getScenario()
      sc$tasks[] <- TRUE
      sc$linearization <- FALSE
      lixoftConnectors::setScenario(sc)
    }else{
      apply_scenario(scenario)
    }

    lixoftConnectors::runScenario()
    lixoftConnectors::computeChartsData()
    lixoftConnectors::saveProject(final_model_file)
  }

  # Generate report with results
  generate_report(project_file = final_model_file, output_file = paste0(toupper(method), "_final_model_report.html"))
  
}

#' Create, fit and report an alternative treatment-effect model
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
alt_mod <- function(out_folder, project_name, 
                    model_file, data_file, 
                    initialvalues, treatment_effects, 
                    headerTypes, error_model, 
                    exclude_random_effects = NULL, normal_dist = NULL, 
                    fixed = NULL, linearize = FALSE){
  
  
  create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, treatment_effects, normal_dist, fixed, linearize)
  
  set_base(initialvalues, fixed = fixed)
  
  lixoftConnectors::saveProject(file.path(out_folder, paste0(project_name, ".mlxtran")))
  
  lixoftConnectors::runScenario()
  
  print(lixoftConnectors::getEstimatedLogLikelihood())
  
  # Generate report with results
  generate_report(project_file = file.path(out_folder, paste0(project_name, ".mlxtran")))
  
}

