#####
# Functions for parameterizing nonlinear mixed effects models using monolix
#####

require(tidyverse)
require(here)
require(Rsmlx)
require(ggbeeswarm)
require(patchwork)
require(tune)
require(gridExtra)
require(GGally)
require(ggpubr)
require(ggpmisc)


#####
## Convergence assessment function:
# Requires a project to be loaded before running
# INPUTS:
# initialValue_ranges: a tibble with ranges of values
# log10: a vector of parameters that should be sampled on a log-scale (NULL if none)
# nits: number of iterations
# fixed: a vector of parameters that should be fixed
# OUTPUTS:
# 

convergence_assessment <- function(initialValue_ranges, log10 = NULL, nits = 100, fixed = NULL){
  
  initialValue_ranges <- initialValue_ranges %>%
    column_to_rownames("parameter")
  
  popparams <- getPopulationParameterInformation()
  
  estimates <- NULL;
  tabestimates <- NULL; tabiters <- NULL; pvals <- NULL; phase <- NULL; initialvalues <- NULL;
  
  run_condition <- TRUE
  i <- 1
  while(run_condition){
    print(i)
    
    # sample new initial estimates
    popini <- sapply(popparams$name, function(j){
      if(j %in% log10){
        10^(runif(1, log10(initialValue_ranges[j,1]), log10(initialValue_ranges[j,2])))
      }else{
        runif(1, initialValue_ranges[j,1], initialValue_ranges[j,2])
      }
    })

    # set sampled values as new initial estimates
    newpopparams <- popparams
    newpopparams$initialValue <- popini
    initialvalues <- rbind(initialvalues, newpopparams %>% mutate(iteration = i))
    
    if(!is.null(fixed)){
      # parameters to estimate as fixed
      newpopparams[newpopparams$name %in% names(fixed),"method"] <- "FIXED"
      newpopparams[newpopparams$name %in% names(fixed),"initialValue"] <- fixed
    }
    
    setPopulationParameterInformation(newpopparams)
    
    # run parameter estimation
    tryCatch({
      runScenario()
      
      # store the estimates and s.e. in a table
      estimates <- as.data.frame(getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      rses <- getEstimatedStandardErrors()$stochasticApproximation$rse
      names(rses) <- getEstimatedStandardErrors()$stochasticApproximation$parameter
      rses <- as.data.frame(rses)
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$run <- i
      estimates$ofv <- getEstimatedLogLikelihood()$importanceSampling[1]
      estimates$aic <- getEstimatedLogLikelihood()$importanceSampling[2]
      estimates$bic <- getEstimatedLogLikelihood()$importanceSampling[3]
      estimates$bicc <- getEstimatedLogLikelihood()$importanceSampling[4]
      names(estimates)[names(estimates) == "Row.names"] <- "param"
      tabestimates <- rbind(tabestimates, estimates)
      
      # store the covariate parameter p values
      if(!is.null(getTests()$wald$p.value)){
        pvals_temp <- as.data.frame(getTests()$wald$p.value)
        colnames(pvals_temp) <- 'pvalue'
        pvals_temp$parameters <- getTests()$wald$parameter
        pvals_temp$run <- i
        pvals <- rbind(pvals, pvals_temp)
      }
      
      # store the iterations
      iters <- getSAEMiterations()$estimates %>% mutate(iteration = 1:nrow(.))
      iters$run <- i
      tabiters <- rbind(tabiters, iters) 
      
      # Store Phase Iteration
      phase <- rbind(phase, getSAEMiterations()$iterationNumbers)
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
    i <- i + 1
    
    # Stop when we reach the number of successful iterations
    if(length(unique(tabestimates$run)) >= nits) run_condition <- FALSE
    
  }
  
  return(
    list(
      tabestimates = tabestimates, 
              tabiters = tabiters, 
              pvals = pvals, 
              phase = phase, 
              initialvalues = initialvalues
      )
    )
}

convergence_assessment2 <- function(initialValue_ranges, log10 = NULL, nits = 100, fixed = NULL){
  
  initialValue_ranges <- initialValue_ranges %>%
    column_to_rownames("parameter")
  
  popparams <- getPopulationParameterInformation()
  
  estimates <- NULL;
  tabestimates <- NULL; tabiters <- NULL; pvals <- NULL; phase <- NULL; initialvalues <- NULL;
  
  run_condition <- TRUE
  i <- 1
  while(run_condition){
    print(i)
    
    # sample new initial estimates
    popini <- sapply(popparams$name, function(j){
      if(j %in% log10){
        10^(runif(1, log10(initialValue_ranges[j,1]), log10(initialValue_ranges[j,2])))
      }else{
        runif(1, initialValue_ranges[j,1], initialValue_ranges[j,2])
      }
    })
    
    # set sampled values as new initial estimates
    newpopparams <- popparams
    newpopparams$initialValue <- popini
    initialvalues <- rbind(initialvalues, newpopparams %>% mutate(iteration = i))
    
    if(!is.null(fixed)){
      # parameters to estimate as fixed
      newpopparams[newpopparams$name %in% names(fixed),"method"] <- "FIXED"
      newpopparams[newpopparams$name %in% names(fixed),"initialValue"] <- fixed
    }
    
    setPopulationParameterInformation(newpopparams)
    
    # run parameter estimation
    tryCatch({
      runScenario()
      
      # store the estimates and s.e. in a table
      estimates <- as.data.frame(getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      rses <- getEstimatedStandardErrors()$stochasticApproximation$rse
      names(rses) <- getEstimatedStandardErrors()$stochasticApproximation$parameter
      rses <- as.data.frame(rses)
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$run <- i
      estimates$ofv <- getEstimatedLogLikelihood()$importanceSampling[1]
      estimates$aic <- getEstimatedLogLikelihood()$importanceSampling[2]
      estimates$bic <- getEstimatedLogLikelihood()$importanceSampling[3]
      estimates$bicc <- getEstimatedLogLikelihood()$importanceSampling[4]
      names(estimates)[names(estimates) == "Row.names"] <- "param"
      tabestimates <- rbind(tabestimates, estimates)
      
      # store the covariate parameter p values
      if(!is.null(getTests()$wald$p.value)){
        pvals_temp <- as.data.frame(getTests()$wald$p.value)
        colnames(pvals_temp) <- 'pvalue'
        pvals_temp$parameters <- getTests()$wald$parameter
        pvals_temp$run <- i
        pvals <- rbind(pvals, pvals_temp)
      }
      
      # store the iterations
      iters <- getSAEMiterations()$estimates %>% mutate(iteration = 1:nrow(.))
      iters$run <- i
      tabiters <- rbind(tabiters, iters) 
      
      # Store Phase Iteration
      phase <- rbind(phase, getSAEMiterations()$iterationNumbers)
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
    i <- i + 1
    
    # Stop when we reach the number of successful iterations
    if(length(unique(tabestimates$run)) >= nits) run_condition <- FALSE
    
  }
  
  return(
    list(
      tabestimates = tabestimates, 
      tabiters = tabiters, 
      pvals = pvals, 
      phase = phase, 
      initialvalues = initialvalues
    )
  )
}

#####
# Set base model with best initial values
# Requires a project to be loaded before running

set_base <- function(initial_values, omega_under_1 = FALSE, fixed = NULL){
  # Unless specified, don't use an initial condition less than 1 for random effects to help with convergence
  if(!omega_under_1){ 
    initial_values <- initial_values %>% 
      rowwise %>% 
      mutate(estimate = ifelse(grepl("omega", param), max(1, estimate), estimate)) 
  }
  
  initial_values <- initial_values %>% 
    column_to_rownames("param") 
  
  popparams <- getPopulationParameterInformation()
  popparams$initialValue <- initial_values[popparams$name, "estimate"]
  
  if(!is.null(fixed)){
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- fixed
  }
  
  setPopulationParameterInformation(popparams)
}


#####
# Model building function
## Can run SCM or COSSAC
## By defualt, BICc is used for the slection criteria (criteria = "BIC")

covariate_model_building <- function(projectFile, method, criteria = "BIC", covariates = "all"){
  
  
  ###### get the default covariate search settings
  set <- getModelBuildingSettings()
  
  # modify the settings
  set$strategy <- method # "cossac" or "scm"
  set$useSambaBeforeCossac <- FALSE # if true, corresponds to "covSAMBA-COSSAC" in the GUI
  set$criterion <- criteria # "BIC" or "LRT"
  set$useLin <- FALSE
  
  if(!any(covariates == "all")){
    set$covariates <- covariates
  }
  
  #####
  # Run search
  runModelBuilding(settings=set)
  model_building_results <- getModelBuildingResults()
  
  #####
  # Generate report for best model
  # Get best model
  best_model <- which(unlist(model_building_results %>% lapply("[[", "bestModel")))
  best_model_file <- sprintf("ModelBuilding_%03d", best_model)
  best_cov <- model_building_results[[best_model]]$individualModels
  
  final_projectFile <- here(gsub(".mlxtran", "", projectFile), "ModelBuilding", sprintf("ModelBuilding_%03d", best_model), paste0(sprintf("ModelBuilding_%03d", best_model), ".mlxtran") )
  loadProject(final_projectFile)
  
  if(any(as.character(unlist(getLaunchedTasks())) == "FALSE")){
    scenario <- getScenario()
    scenario$tasks[] <- TRUE
    scenario$linearization <- FALSE
    setScenario(scenario)
    
    runScenario()
    computeChartsData()
    
    saveProject( here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran"))  
  }
  
  # Generate report with results
  rmarkdown::render(
    here('scripts', 'Monolix_template.Rmd'), 
    params = list(project = here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran")), 
    output_file = here(gsub(".mlxtran", "", projectFile), paste0(method, "_final_report.html"))
  )
  
}

#####
# Define project settings and scenario

project_settings <- function(){
  # Set scenario for Monolix
  scenario <- getScenario()
  scenario$tasks["standardErrorEstimation"] <- TRUE
  scenario$tasks["logLikelihoodEstimation"] <- TRUE
  scenario$linearization <- FALSE
  setScenario(scenario)
  
  # Set Parameter Estimation Settings
  setPopulationParameterEstimationSettings(
    nbBurningIterations=50,
    nbExploratoryIterations=1000,
    nbSmoothingIterations=500,
    smoothingAlpha=0.7,
    exploratoryAutoStop = TRUE,
    smoothingAutoStop = TRUE
  )
  
}

#####
# Set up a new project

create_project <- function(model_file, data_file, headerTypes, error_model, exclude_random_effects, treatment_effects = NULL, normal_dist = NULL){
  
  # Create a new project
  newProject(modelFile = model_file,
             data = list(dataFile = here("data", data_file),
                         headerTypes = headerTypes
             )
  )
  
  # Set Error model
  # Assumes that observation column is named "vl"
  setErrorModel(vl = error_model)
  
  # Remove random effects if needed
  if(!is.null(exclude_random_effects)){
    var_model <- getIndividualParameterModel()$variability$id
    var_model[exclude_random_effects] <- FALSE
    setIndividualParameterVariability(var_model)
  }
  
  if(!is.null(treatment_effects)){
    # treatment_effects should be a list of vectors of categorical variables to add to parameters
    covariateModel <- getIndividualParameterModel()$covariateModel
    for(parameter in names(treatment_effects)){
      covariateModel[[parameter]][treatment_effects[[parameter]]] <- TRUE
    }
    setCovariateModel(covariateModel)
  }
  
  if(!is.null(normal_dist)){
    distributions <- getIndividualParameterModel()$distribution
    distributions[normal_dist] <- "normal"
    setIndividualParameterDistribution(distributions)
  }
  
  print(getIndividualParameterModel())
  
  project_settings()
  
}

#####
# Find narrower parameter ranges
get_narrow_ranges <- function(tabestimates){
  ## Find 25th percentile 
  q25 <- quantile(distinct(tabestimates, run, bicc)$bicc, 0.25)
  
  initialValue_ranges_narrow <- tabestimates %>% 
    group_by(param) %>% 
    filter(!outlier_detector(estimate)) %>% # Remove estimate outliers
    filter(bicc <= q25)  %>%  # subset to best BICcs
    summarise(min = min(estimate, na.rm = T), max = max(estimate, na.rm = T)) %>% 
    rename(parameter = param)
  
  initialValue_ranges_narrow
  
}

#####
# Function to detect outliers
outlier_detector <- function(value, upper_only = F){
  if(upper_only){
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) 
  }else{
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) |
      (value < quantile(value, probs = 0.25, na.rm = T) - 1.5*IQR(value, na.rm = T) )  
  }
  
}

#####
# Data formatting for visualization
filter_results <- function(tabestimates, remove_estimate_outliers, remove_bicc_outliers){
  
  ## Calculate 25th percentile of biCcs:
  q25 <- quantile(distinct(tabestimates, run, bicc)$bicc, 0.25)
  
  ## Filter data
  caption <- ""
  # Runs with bicc outliers
  if(remove_bicc_outliers){
    remove_runs <- tabestimates %>% distinct(run, bicc) %>% filter(outlier_detector(bicc, upper_only = T)) %>% pull(run)
    if(length(remove_runs) >= 1){
      caption <- paste0(caption, "Runs with outlier BICcs removed\n")  
    }else{
      cat("No runs with outlier BICc to removed\n")
      caption <- caption
    }
    
  }else{
    remove_runs <- c()  
  }
  
  # Remove Estimate outliers 
  if(remove_estimate_outliers){
    tabestimates <- tabestimates %>% group_by(param) %>% filter(!outlier_detector(estimate)) %>% ungroup()
    caption <- paste0(caption, "Estimate outliers removed")
  }
  
  # Remove runs with bicc outliers
  tabestimates <- tabestimates %>% filter(!run %in% remove_runs)
  
  return(list(tabestimates = tabestimates, caption = caption, q25 = q25))
}


#####
# Plot BIC Versus estimates

plot_bicc_v_estimate <- function(tabestimates, out_folder = "", project_name = "", narrow = TRUE,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE, log10_x = TRUE, 
                                 treatment_effects = FALSE, ylim = NULL,
                                 file_addendum = "", save = TRUE, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  q25 <- filtered_results$q25
  
  if(!is.null(parameter)){
    tabestimates <- tabestimates %>% filter(param == parameter)
  }else if(!treatment_effects){
    tabestimates <- tabestimates %>% filter(!grepl("beta_", param))
  }else{
    tabestimates <- tabestimates %>% filter(grepl("beta_", param))
  }
  
  if(!is.null(ylim)){
    tabestimates <- tabestimates %>% filter(bicc <= ylim[2])
  }
  
  ## Generate plot
  popparams <- getPopulationParameterInformation()
  
  
  if(narrow){
    estimates_vs_bicc <- tabestimates %>%
      mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot(aes(estimate, bicc, shape = bicc == min(bicc), fill = bicc)) +
      geom_point(size = 2, alpha = 0.75) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) 
  }else{
    estimates_vs_bicc <- tabestimates %>%
      mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot(aes(estimate, bicc,  fill = bicc)) +
      geom_point(size = 2, alpha = 0.75, shape = 21) +
      geom_text(data = . %>% filter(estimate == min(estimate, na.rm = T)), 
                aes(y = q25, x = estimate, label = paste0("25th\npercentile:\n", round(q25))),
                hjust = 0, vjust = -0.1, size = 3) + 
      geom_hline(yintercept = q25, lty = "dashed", alpha = 0.5) + 
      labs(caption = paste0("Dashed line indicates 25th percentile of BICc\n", caption))
  }
  estimates_vs_bicc <- estimates_vs_bicc + 
    facet_wrap(~param, scales = "free_x") + 
    labs(x = "Final estimate", y = "BICc", shape = "Best fit", fill = "BICc",
         title = "Comparing BIC and parameter estimates",
         caption = caption) + 
    scale_shape_manual(values = c(21, 24)) + 
    scale_fill_viridis_c(direction = -1)
  
  if(any(grepl("beta_", tabestimates$param))) estimates_vs_bicc <- estimates_vs_bicc + labs(subtitle = ifelse(treatment_effects, "treatment effects only", "treatment effects excluded"))
  if(!is.null(ylim)) estimates_vs_bicc <- estimates_vs_bicc + ylim(ylim) + labs(caption = paste0(caption, "\nBICc greater than ", max(ylim), " removed"))
  if(!is.null(parameter)){
    if(!grepl("beta_", parameter) & log10_x) estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  }else if(log10_x & !treatment_effects){ 
    estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  }
  
  
  if(save){
    file_name <- ifelse(narrow, 
                        ifelse(treatment_effects, "BIC_vs_estimates_narrow_range_treatment_effects", "BIC_vs_estimates_narrow_range"),
                        ifelse(treatment_effects, "BIC_vs_estimates_treatment_effects", "BIC_vs_estimates"))
    file_name <- paste0(file_name, file_addendum, ".pdf")
    
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), file_name), estimates_vs_bicc, width = 10, height = 7)
    
  }
  
  return(estimates_vs_bicc)
  
}

#####
# Plot initial value and final estimate distributions

plot_distributions <- function(initialvalues, tabestimates, out_folder = "", project_name = "", parameter_set = "pop",
                               remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                               save = T, parameter = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  subtitles <- c("pop" = "Population parameters",
                 "omega" = "Standard deviations of the random effect",
                 "beta" = "Treatment effects",
                 "error" = "Error model parameters")
  subtitle <- paste(subtitles[parameter_set], collapse = " and ")
  
  if(is.null(parameters)){
    params <- unlist(lapply(parameter_set, function(z){ 
      p <- c()
      if(z == "error"){
        p <- c(p, unique(tabestimates$param)[unique(tabestimates$param) %in% c("a", "b", "c")])
      }else{
        if(z == "beta") z <- paste0(z, "_")
        p <- c(p, unique(tabestimates$param)[grepl(z, unique(tabestimates$param))])
      }
      return(p)
    }))
    
    df <- merge(initialvalues %>% select(param = name, run = iteration, initialValue), tabestimates) %>% 
      filter(param %in% params)
  }else{
    df <-  merge(initialvalues %>% select(param = name, run = iteration, initialValue), tabestimates) %>% 
      filter(param == parameter) 
  }
  
  popparams <- getPopulationParameterInformation()
  
  plot <- df %>% 
    pivot_longer(c(initialValue, estimate)) %>% 
    mutate(name = factor(name, levels = c("initialValue", "estimate"), labels = c("Initial Value", "Final Estimate")),
           param = factor(param, levels = popparams$name)) %>% 
    ggplot(aes(name, value)) + geom_violin(scale = "width") + 
    facet_wrap(~param, scales = "free_y") +
    geom_boxplot(position = position_dodge(width=0.9), width = 0.1, fill = "white", outlier.shape = NA) +
    geom_beeswarm(aes( fill = bicc, color = bicc, group = name), size = 2, alpha = 0.75, shape = 21) + 
    scale_fill_viridis_c(direction = -1) +
    scale_color_viridis_c(direction = -1) +
    geom_point(data = . %>% filter(bicc == min(bicc)), aes(fill = bicc), shape = 24, size = 3, alpha = 0.75) +
    labs(x = "", y = "Parameter Value", fill = "BICc",
         title = "Comparing initial value and final estimate distributions",
         caption = caption,
         subtitle = subtitle) + 
    guides(color = "none")
  
  if(!is.null(parameter)){
    if(grepl("beta_", parameter) | grepl("omega_", parameter)){
      ylim <- layer_scales(plot)$y$range$range
      if(min(ylim) > 0) ylim[1] <- 0
      if(min(ylim) < 0 & max(ylim) < 0) ylim[2] <- 0
      plot <- plot + coord_cartesian(ylim = ylim)
    }
  }
  
  if(save){
    file_name <- paste0("iv_estimate_dist_", paste(parameter_set, collapse = "_"), ".pdf")
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), file_name), plot, width = 10, height = 7)  
  }
  
  return(plot)
}

#####
# Plot final estimates versus initial values

plot_iv_vs_estimates <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                                 remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                                 save = T, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  df <- merge(initialvalues %>% select(param = name, run = iteration, initialValue), tabestimates) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    arrange(param)
  
  params <- as.list(unique(df$param))
  
  if(is.null(parameter)){
    plots <- lapply(params, function(p){
      df %>% filter(param == p) %>% 
        ggplot(aes(initialValue, estimate, fill = bicc, shape = bicc == min(bicc))) + 
        geom_point(size = 2, alpha = 0.75) + 
        geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) + 
        coord_obs_pred() +
        labs(title = p) + 
        theme(plot.title = element_text(hjust = 0.5))
    })
    
    color_limit <- range(df$bicc)
    
    combined_plot <- wrap_plots(plots)
    
    plot <- combined_plot + 
      plot_layout(guides = "collect", axis_titles = "collect") &
      labs( color = "BICc", shape = "Best fit") &
      plot_annotation(title = "Comparing initial value and final estimates", caption = caption) &
      theme(axis.title = element_blank()) &
      scale_shape_manual(values = c(21, 24)) &
      scale_fill_viridis_c(direction = -1, limits = color_limit)
    
    plot <- grid.arrange(patchworkGrob(plot), left = "Final Estimate", bottom = "Initial Value")
  }else{
    color_limit <- range(df$bicc)
    
    plot <- df %>% 
      filter(param == parameter) %>% 
      ggplot(aes(initialValue, estimate, fill = bicc, shape = bicc == min(bicc))) + 
      geom_point(size = 2, alpha = 0.75) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) + 
      coord_obs_pred() +
      labs(color = "BICc", y = "Final Estimate", x = "Initial Value", shape = "Best fit") + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_shape_manual(values = c(21, 24)) + 
      scale_fill_viridis_c(direction = -1, , limits = color_limit) + 
      ggpmisc::stat_correlation(method = "pearson",  label.x = 0, label.y = 1,
                                use_label("R", "P", sep = "\n"),
                                output.type = "text")
  }
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "estimates_vs_ivs.pdf"), plot, width = 12, height = 10)  
  }
  
  return(plot)
}

#####
# Correlate final and initial values

cor_iv_estimates <- function(tabestimates, initialvalues, out_folder, project_name,
                             remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, save = T){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  df <- merge(initialvalues %>% select(param = name, run = iteration, initialValue), tabestimates) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    group_by(param) %>% 
    summarise(cor = cor.test(initialValue, estimate)$estimate, pval = cor.test(initialValue, estimate)$p.value)
  
  breaks <- round(range(log10(df$pval)))
  breaks <- seq(breaks[1], breaks[2], by = 2)
  
  plot <- df %>% 
    ggplot(aes(cor, fct_rev(param), fill = log10(pval), shape = pval <= 0.05)) + 
    geom_vline(xintercept = 0, lty = "dashed") + 
    geom_point(size = 3) + 
    xlim(c(-1,1)) + 
    labs(x = "Correlation", y = "Parameter", caption = caption, 
         title = "Correlation between initial values and final estimates") +
    scale_fill_viridis_c(option = "plasma", direction = -1, breaks = breaks) + 
    scale_shape_manual(values = c(21, 24))
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "cor_iv_estimates.pdf"), plot, width = 7, height = 8)  
  }
  
  return(plot)
}


#####
# Plot parameters against eachother

correlate_estimates <- function(tabestimates, out_folder = "", project_name = "",
                                remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, save = T, file_addendum = ""){
  
  # Filter tabestimates to remove estimate outliers and/or BICc outliers
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  correlation_df <- tabestimates %>% 
    select(param, estimate, run)  %>% 
    pivot_wider(names_from = param, values_from = estimate) %>% 
    select(-run) %>% 
    cor(use = "complete.obs") %>% # only keep runs that have no estimate outliers (if estimate outliers were removed)
    as.data.frame %>% 
    rownames_to_column("parameter2") %>% 
    pivot_longer(!parameter2, names_to = "parameter1") 
  
  correlation_heatmap <- correlation_df %>% 
    mutate(parameter2 = factor(parameter2, levels = popparams$name),
           parameter1 = factor(parameter1, levels = popparams$name)) %>%
    ggplot(aes(parameter1, parameter2, fill = value)) + 
    geom_tile() + 
    geom_text(data = . %>% filter(value > 0.5 | value < -0.5),
              aes(label = round(value, 2))) +
    scale_fill_gradient2(limits = c(-1,1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) + 
    labs(fill = "correlation", title = "Correlation of SAEM estimates", caption = caption)
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), paste0("SAEM_correlations", file_addendum, ".pdf")), correlation_heatmap, width =10, height = 10)  
  }
  
  return(correlation_heatmap)
  
}

### Plot population estiamates against eachother
compare_estimates <- function(tabestimates, out_folder = "", project_name = "",
                              remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, save = T, file_addendum = ""){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  plot <- tabestimates %>% 
    select(param, estimate, run, bicc)  %>% 
    filter(grepl("pop", param)) %>% 
    pivot_wider(names_from = param, values_from = estimate) %>% 
    GGally::ggpairs(columns = 3:ncol(.), 
                    lower = list(continuous = function(data, mapping){
                      ggplot(data, mapping) + 
                        geom_point(aes(fill = bicc, shape = bicc == min(bicc)), alpha = 0.75) + 
                        geom_point(data = . %>% filter(bicc == min(bicc)), aes(fill = bicc, shape = bicc == min(bicc)), size = 2) + 
                        scale_fill_viridis_c(direction = -1) + 
                        scale_shape_manual(values = c(21, 24)) + 
                        labs(fill = "BICc", shape = "Best fit") 
                    } 
                    )
    ) + 
    labs(caption = caption) 
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), paste0("paired_SAEM_population_estimates", file_addendum, ".pdf")), plot, width = 7, height = 7)  
  }
  
  return(plot)
  
}

#####
# Calculate coefficient of variance

calc_CV <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                    remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                    save = T, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  cv <- merge(
    tabestimates, 
    initialvalues %>% rename(param = name, run = iteration)
  )  %>% 
    group_by(param) %>% 
    summarise(`Final Estimate` = sd(estimate)/mean(estimate),
              `Initial Value` = sd(initialValue)/mean(initialValue)) 
  
  if(!is.null(parameter)) cv <- cv %>% filter(param == parameter)
  
  cv <- cv %>% pivot_longer(!param) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    ggplot(aes(factor(param), abs(value), fill = name)) + 
    geom_col(position = position_dodge(),alpha = 0.7) + 
    geom_text(aes(label = abs(round(value, 2))), position = position_dodge(width = 0.91), hjust = -0.1) + 
    labs(x = "Parameter", y = "Absolute Value of Coefficient of Variation", caption = caption, title = "Comparing CV of initial values and final estimates") + 
    coord_flip() + 
    theme(legend.position.inside = c(0.9,0.9), legend.title = element_blank(), 
          legend.position = "inside", legend.background = element_rect(color = "white")) + 
    scale_fill_brewer(palette = "Set1") 
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "coefficient_of_variation.pdf"), cv, width = 8.5, height = 10)  
  }
  
  return(cv)
}



#####
# Generate report for an alternate model

alt_mod <- function(treatment_effects, project_name, initialvalues,
                    model_file, data_file, headerTypes, error_model, exclude_random_effects, out_folder){
  
  create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, treatment_effects)
  
  set_base(initialvalues, treatment = 0)
  
  saveProject(here(out_folder, paste0(project_name, ".mlxtran")))
  
  runScenario()
  
  print(getEstimatedLogLikelihood())
  
  # Generate report with results
  rmarkdown::render(
    here('scripts', 'Monolix_template.Rmd'), 
    params = list(project = here(out_folder, paste0(project_name, ".mlxtran"))), 
    output_file = here(out_folder, project_name, "report.html")
  )
  
}

###
# Plot SAEM iterations

plot_SAEM <- function(tabiters, phase, tabestimates, out_folder = "", project_name = "",
                      remove_bicc_outliers = FALSE, save = T, filter = NULL, parameter = NULL){
  
  # browser()
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  # caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  
  if(is.null(parameter)){
    plotList <- list()
    i <- 1
    for (param in c(popparams$name, "convergenceIndicator")) {
      if(param != "convergenceIndicator"){
        if (popparams[popparams$name == param, ]$method == "FIXED") next
      }
      p <- param
      
      temp <- tabiters %>% 
        left_join(data.frame(run = unique(tabiters$run), phase = phase[,1])) %>% 
        filter(run %in% unique(tabestimates$run)) %>% 
        left_join(tabestimates %>% distinct(run, bicc)) 
      
      if(!is.null(filter)) temp <- temp %>% filter(iteration-phase >= filter)
      
      plotList[[i]] <- temp %>% 
        ggplot( aes(x = iteration-phase, y = !!sym(param))) +
        geom_line(aes(group = run, color = bicc)) +
        theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
        # geom_vline(xintercept = changePhase, color = 1:length(changePhase)) +
        labs(title = param, x = NULL, y = NULL) + 
        # scale_y_log10() + 
        # coord_cartesian(xlim = c(-200, 314)) + 
        geom_vline(xintercept = 0, lty = "dashed", color = "black") + 
        geom_hline(yintercept = tabestimates %>% filter(bicc == min(bicc), param == p) %>% pull(estimate),
                   alpha = 0.5) + 
        scale_color_viridis_c(direction = -1) 
      i <- i + 1
    }
    
    plot <- grid.arrange(grobs = plotList, ncol = 5)
  }else{
    
    
    temp <- tabiters %>% 
      left_join(data.frame(run = unique(tabiters$run), phase = phase[,1])) %>% 
      filter(run %in% unique(tabestimates$run)) %>% 
      left_join(tabestimates %>% distinct(run, bicc)) 
    
    
    if(!is.null(filter)) temp <- temp %>% filter(iteration-phase >= filter)
    
    plot <- temp %>% 
      ggplot( aes(x = iteration-phase, y = !!sym(parameter))) +
      geom_line(aes(group = run, color = bicc)) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs( x = "Iteration - phase switch", y = "Estimate") + 
      geom_vline(xintercept = 0, lty = "dashed", color = "black") + 
      geom_hline(yintercept = tabestimates %>% filter(bicc == min(bicc), param == parameter) %>% pull(estimate),
                 alpha = 0.5) + 
      scale_color_viridis_c(direction = -1) 
    
    ## adjust y axis if SAEM portion is much more variable than estimates
    us <- range(temp[parameter])
    ue <- range(tabestimates %>% filter(param == parameter) %>% pull(estimate))
    lower <- min(c(us, ue))
    if(log10(abs(us[2]/ue[2])) > 1) plot <- plot + coord_cartesian(ylim = c(lower, ue[2]*1.2))
  }
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "SAEM_trajectories.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

plot_estimates <- function(tabiters, phase, tabestimates, out_folder = "", project_name = "",
                           remove_bicc_outliers = FALSE, save = T, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  
  popparams <- getPopulationParameterInformation()
  
  # Get plot y-limit, excluding unusally large rses:
  ylim_range <-  tabestimates %>% 
    group_by(param) %>%
    filter(!outlier_detector(rses)) %>% 
    mutate(upper = estimate * (1+as.numeric(rses)/100), lower = estimate * (1-as.numeric(rses)/100)) %>% 
    summarise(lower = min(lower), upper = max(upper)) %>% 
    column_to_rownames("param")
  
  if(is.null(parameter)){
    # plot population parameters
    plotList <- list()
    i <- 1
    for (param in popparams$name) {
      if (popparams[popparams$name == param, ]$method == "FIXED") next
      estimates <- tabestimates[tabestimates$param == param, ]
      plotList[[i]] <- ggplot(estimates %>% filter(run %in% unique(tabestimates$run)), 
                              aes(x = run, y = estimate, color = bicc)) +
        geom_point() +
        geom_errorbar(aes(ymax = estimate * (1+as.numeric(rses)/100), 
                          ymin = estimate * (1-as.numeric(rses)/100))) +
        theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
        labs(title = param, x = NULL, y = NULL) +
        geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate),
                   alpha = 0.5) + 
        scale_color_viridis_c(direction = -1) 
      i <- i + 1
    }
    
    plot <- grid.arrange(grobs = plotList, ncol = 5)
  }else{
    estimates <- tabestimates[tabestimates$param == parameter, ]
    plot <- ggplot(estimates %>% filter(run %in% unique(tabestimates$run)), 
                   aes(x = run, y = estimate, color = bicc)) +
      geom_point() +
      geom_errorbar(aes(ymax = estimate * (1+as.numeric(rses)/100), 
                        ymin = estimate * (1-as.numeric(rses)/100))) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs(x = NULL, y = "Final Estimate") +
      geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate), alpha = 0.5) + 
      scale_color_viridis_c(direction = -1) +
      coord_cartesian(ylim = unname(unlist(ylim_range[parameter,])))
  }
  
  if(grepl("beta_", parameter) | grepl("omega_", parameter)){
    ylim <- layer_scales(plot)$y$range$range
    if(min(ylim) > 0) ylim[1] <- 0
    if(min(ylim) < 0 & max(ylim) < 0) ylim[2] <- 0
    plot <- plot + coord_cartesian(ylim = ylim)
  }
  
  if(save){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "Convergence_estimates.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

###
# Function for combined plots

parameter_plots <- function(parameter, remove_bicc_outliers, remove_estimate_outliers,
                            log10_tabestimates, tabestimates, tabiters,
                            phase, initialvalues){
  
  tabestimates <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
  tabestimates <- tabestimates %>% filter(param == parameter)
  
  if(!is.null(log10_tabestimates)){
    log10_tabestimates <- filter_results(log10_tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
    log10_tabestimates <- log10_tabestimates %>% filter(param == parameter)
  }
  
  SAEM_traj <- plot_SAEM(
    tabiters, phase, tabestimates, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    save = F, parameter = parameter#, filter = -200
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  SAEM_convergence <- plot_estimates(
    tabiters, phase, tabestimates,
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    save = F, parameter = parameter
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    labs(x = "Run")
  
  narrow_range <- tabestimates %>% filter(param == parameter) %>% 
    pull(estimate) %>% range
  
  
  if(!is.null(log10_tabestimates)){
    wide_plot <- plot_bicc_v_estimate(
      log10_tabestimates, narrow = FALSE, 
      # remove_bicc_outliers =  remove_bicc_outliers, 
      # remove_estimate_outliers = remove_estimate_outliers,
      save = F, parameter = parameter
    ) + 
      scale_fill_viridis_c(direction = -1, option = "magma",
                           name = "BICc wide\nparameter\nsearch")  +
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) + 
      annotate("rect", xmin = narrow_range[1], xmax = narrow_range[2], 
               ymin = -Inf, ymax = Inf,
               alpha = 0.3, fill = "lightgray")
  }
  
  
  narrow_plot <- plot_bicc_v_estimate(
    tabestimates,  narrow = TRUE, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    log10_x = F, save = F, parameter = parameter
  ) + 
    labs(fill = "BICc narrow\nparameter\nsearch") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  iv_vs_final <- plot_iv_vs_estimates(
    tabestimates, initialvalues, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    save = F, parameter = parameter
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  # Get coefficients of variation
  cv <- merge(
    tabestimates, 
    initialvalues %>% rename(param = name, run = iteration)
  )  %>% 
    group_by(param) %>% summarise(`Final Estimate` = sd(estimate)/mean(estimate),
                                  `Initial Value` = sd(initialValue)/mean(initialValue)) %>% 
    pivot_longer(!param)
  
  distributions <- plot_distributions(
    initialvalues, tabestimates, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    parameter = parameter, save = F
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch", 
         x = "Set") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    geom_text(data = cv, aes(label = round(value, 2), y = Inf), vjust = 1) + 
    annotate(geom = "text", label = "CV:", y = Inf, x = -Inf, vjust = 1, hjust = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  cv <- calc_CV(
    tabestimates, initialvalues, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers,
    save = F, parameter = parameter
  )
  
  caption <- ""
  if(remove_bicc_outliers) caption <- paste0(caption, "Runs with BICc outliers removed\n")
  if(remove_estimate_outliers) caption <- paste0(caption, "Estimate outliers removed")
  
  if(!is.null(log10_tabestimates)){
    final_plot <- wide_plot +  iv_vs_final +  SAEM_traj + narrow_plot + distributions + SAEM_convergence +
      plot_layout(nrow = 3, byrow = F, guides = "collect", widths = c(1, 1.1)) &
      plot_annotation(title = parameter, theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = caption) &
      labs(caption = "") 
  }else{
    final_plot <- plot_spacer() + iv_vs_final +  SAEM_traj + narrow_plot + distributions + SAEM_convergence +
      plot_layout(nrow = 3, byrow = F, guides = "collect", widths = c(1, 1.1)) &
      plot_annotation(title = parameter, theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = caption) &
      labs(caption = "") 
  }
  
  return(final_plot)
  
}

