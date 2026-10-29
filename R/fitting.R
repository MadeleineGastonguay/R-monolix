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
#' @param predefined_initialvalues Optional table of initial values (columns
#'   `parameter`, `run`, `initial_value`), e.g. shared across candidate models in
#'   step 4. When given, each predefined set is fitted exactly once, under its own
#'   run number, and `nits` is ignored. Failed runs are not replaced, so paired
#'   comparisons only use the runs that succeeded in every model.
#' @param set_error Named starting values for error-model parameters (see
#'   [error_starts()]). They are still estimated; only the start is set. If
#'   `NULL`, error parameters follow `set_variance` or their sampling range.
#' @param set_treatment Starting value for every `beta_*` treatment effect.
#'   Without it, treatment effects start at 0 (they have no sampling range).
#' @param seed Seed for sampling initial values (`NULL` leaves the RNG alone).
#' @param initialValue_ranges Sampling ranges: a data frame with one row per
#'   parameter (row names are the parameter names) and columns `min` and `max`.
#' @param log10 Names of the parameters sampled uniformly on a log10 scale;
#'   the others are sampled uniformly.
#' @param nits Number of runs (failed runs are replaced, up to a limit).
#' @param set_variance Starting value of every random-effect SD (`omega_`),
#'   the same in every run; `NULL` samples them from their range.
#' @param record_shrinkage Record eta shrinkage per run (`shrinkage`).
#' @param fixed Named values of parameters to fix, e.g.
#'   `c(param1_pop = 10)`.
#' @param record_pvals Record the Wald-test p-values of the treatment effects
#'   per run (`pvals`).
#' @param record_indiv_params Record the individual parameters
#'   (`indivParams`) and individual fits (`indivFit`) per run.
#' @param record_condition_number Record the condition number of the
#'   correlation matrix of the estimates per run (`condition_number`).
#' @param record_table Record Monolix's `tables.txt` per run
#'   (`table_predictions`), for models with a `table = {...}` output.
#' @return List with `tabestimates` (estimates, RSEs and likelihood of every
#'   run), `initialvalues`, `tabiters` (SAEM iterations), `phase` and the
#'   records asked for. After the runs, the project is left at the best run's
#'   starting values, fitted again.
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
  # Which observations are censored, joined onto each run's predictions (predictions.txt
  # doesn't say); NULL when the data has no cens column
  censoring <- if (record_indiv_params) tryCatch({
    d <- lixoftConnectors::getData()
    censoring_table(read_monolix_data(d$dataFile), d$header, d$headerTypes)
  }, error = function(e) NULL)

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
        pred_i <- read_delim(file.path(project_results_dir(), "predictions.txt"), show_col_types = FALSE)
        if (!is.null(censoring) && all(names(censoring)[1:2] %in% names(pred_i))) {
          pred_i <- left_join(pred_i, censoring, by = names(censoring)[1:2])
        }
        predictions <- rbind(predictions, pred_i %>% mutate(run = i))
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

# Round values for log10 grids: within each factor of 10, the first of these lists that gives
# enough values. 1-2-5, then denser.
LOG_GRID_STEPS <- list(c(1, 2, 5), c(1, 2, 3, 5, 7), c(1, 1.5, 2, 3, 4, 5, 7),
                       c(1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8),
                       c(1, 1.1, 1.2, 1.3, 1.5, 1.7, 2, 2.2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9))

#' Grid of values for a profile likelihood
#'
#' Round numbers between `range[1]` and `range[2]`, plus both ends, with at
#' least `n` values:
#' * linear: a round step (1, 2, 2.5 or 5 times a power of 10), the largest
#'   that gives `n` values, e.g. 0.5, 1, 2, ..., 20;
#' * log10: within each factor of 10, values from 1-2-5 or a denser list of
#'   round numbers (1, 1.5, 2, 3, 4, 5, 7, ...), the sparsest that gives `n`
#'   values, e.g. 0.5, 0.7, 1, 1.5, 2, 3, 4, 5, 7, 10, 15, 20. A range too
#'   narrow for that (e.g. 90 to 110) uses a round linear step instead.
#'
#' @param range `c(min, max)`.
#' @param n Smallest number of values; the grid can have a few more.
#' @param log10 Round numbers on a log10 scale (`TRUE`) or a linear one.
#' @return Sorted numeric vector.
#' @export
profile_grid <- function(range, n, log10 = FALSE) {
  lo <- range[1]; hi <- range[2]
  with_ends <- function(x) sort(unique(signif(c(lo, x[x > lo & x < hi], hi), 12)))
  if (isTRUE(log10)) {
    decades <- 10^(floor(base::log10(lo)):ceiling(base::log10(hi)))
    for (steps in LOG_GRID_STEPS) {
      g <- with_ends(as.vector(outer(steps, decades)))
      if (length(g) >= n) return(g)
    }
  }
  # linear, or a log range too narrow for the round-number lists
  k <- floor(base::log10(hi - lo))
  steps <- sort(as.vector(outer(c(1, 2, 2.5, 5), 10^((k - 6):(k + 1)))), decreasing = TRUE)
  for (step in steps) {
    g <- with_ends(seq(ceiling(lo / step), floor(hi / step)) * step)
    if (length(g) >= n) return(g)
  }
  g
}

#' Profile likelihood for one population parameter
#'
#' Fits the loaded project once per value of `parameter` on a grid of round
#' numbers over `range` ([profile_grid()]), or on `values`, with `parameter` fixed at that value and everything else estimated.
#'
#' @param parameter Population parameter to profile, e.g. `"param1_pop"`.
#' @param range `c(min, max)` of the grid.
#' @param igs Initial values: a table with `name`/`initialValue` (or
#'   `parameter`/`initial_value`), or a table of estimates (`parameter`,
#'   `estimate`), whose non-`_pop` values are raised to at least 1.
#' @param nits Smallest number of grid values (the grid can have a few more).
#' @param linearize `NULL` (default) uses the project's scenario; `TRUE` or
#'   `FALSE` overrides it.
#' @param log10 With `values = NULL`, `TRUE` makes the grid round numbers on a
#'   log10 scale, `FALSE` (default) on a linear one ([profile_grid()]).
#' @param values The grid itself, instead of `range`, `nits` and `log10`.
#' @return Table with one row per grid value and estimated parameter
#'   (`parameter`, `estimate`, `rse`, `profiled_parameter`, `profiled_value`,
#'   `ofv`, `aic`, `bic`, `bicc`).
#' @param use_estimates_as_ig Start each grid value's fit from the previous
#'   fit's estimates instead of from `igs`.
#' @export
profile_likelihood <- function(parameter, range, igs, nits, use_estimates_as_ig = F, linearize = NULL,
                               log10 = FALSE, values = NULL){
  
  if (is.null(values)) values <- profile_grid(range, nits, log10)
  nits <- length(values)
  
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
      
      if(use_estimates_as_ig){
        initialValues <- estimates %>% select(name = parameter, initialValue = estimate) %>% 
          mutate(method = "MLE", initialValue = ifelse(grepl("pop", name), initialValue, pmax(initialValue, 1))) 
        initialValues <- rbind(initialValues, data.frame(name = parameter, initialValue = value, method = "FIXED"))
      }
      
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
  }
  
  return(tabestimates)
  
}

#' Run Monolix covariate model building (SCM or COSSAC)
#'
#' Results go to `<project>/ModelBuilding` and `<project>/final_model.mlxtran`, so to run
#' both SCM and COSSAC, give each its own `projectFile` (e.g. `scm.mlxtran`
#' and `cossac.mlxtran`).
#'
#' @param report_config Optional config, for the report's grouping and labels.
#' @return The model-building results (`getModelBuildingResults()`), invisibly.
#' @param scenario Optional `scenario` section of the config, used when the
#'   best model's tasks have to be rerun. `NULL` runs every task without
#'   linearization.
#' @param projectFile Project to search from (saved first if it doesn't
#'   exist). Stops if its `ModelBuilding` folder already exists.
#' @param method `"SCM"` or `"COSSAC"`.
#' @param criteria `"BIC"` (Monolix's name for BICc) or `"LRT"`.
#' @param covariates Covariates to search, or `"all"`.
#' @param linearization Use linearization during the search.
#' @param relationships Parameter-covariate pairs to keep in every model or in
#'   none, from [model_building_relationships()]; `NULL` for none.
#' @param parameters Parameters to search; `NULL` for all.
#' @param lrt_thresold With `criteria = "LRT"`, the p-value thresholds for
#'   adding and for removing a covariate.
#' @export
covariate_model_building <- function(projectFile, method = "SCM", criteria = "BIC", covariates = "all", 
                                     linearization = FALSE, relationships = NULL, parameters = NULL, lrt_thresold = c(0.05,0.01),
                                     scenario = NULL, report_config = NULL){
  
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
  generate_report(project_file = final_model_file, output_file = paste0(toupper(method), "_final_model_report.html"),
                  config = report_config)

  invisible(model_building_results)
}

#' Create, fit and report an alternative treatment-effect model
#' @param report_config Optional config, for the report's grouping and labels
#'   (see [generate_report()]).
#' @inheritParams create_project
#' @param out_folder Folder to save the project in.
#' @param project_name Project name; it's saved as
#'   `<out_folder>/<project_name>.mlxtran`.
#' @param initialvalues Starting values, as for [set_base()]
#'   (`parameter`, `estimate`).
#' @return Called for its effects: the fitted project and its report.
#' @export
alt_mod <- function(out_folder, project_name, 
                    model_file, data_file, 
                    initialvalues, treatment_effects, 
                    headerTypes, error_model, 
                    exclude_random_effects = NULL, normal_dist = NULL, 
                    fixed = NULL, linearize = FALSE, report_config = NULL){
  
  
  create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, treatment_effects, normal_dist, fixed, linearize)
  
  set_base(initialvalues, fixed = fixed)
  
  lixoftConnectors::saveProject(file.path(out_folder, paste0(project_name, ".mlxtran")))
  
  lixoftConnectors::runScenario()
  
  print(lixoftConnectors::getEstimatedLogLikelihood())
  
  # Generate report with results
  generate_report(project_file = file.path(out_folder, paste0(project_name, ".mlxtran")), config = report_config)
  
}

