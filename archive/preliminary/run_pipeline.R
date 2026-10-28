#####
# This script runs model building and convergence assessment for the viral rebound kinetics model
#
# INPUTS:
# - dataset 
# - type for each column in the dataset
# - structural model file
# - random effects to exclude
# - type of error model
# - ranges (over orders of magnitude) of biologically plausible values for initial population parameter values
# - name of folder for results to be saved in (default is Monolix runs/data_file/)
# - criteria for model building (default: BICc)
#
# OUTPUTS:
# - best initial estimates for population parameters
# - best covariate model from model building with SCM and COSSAC
# - convergence assessment for best fit models
#
#####

run_pipeline <- function(data_file, headerTypes, model_file, exclude_random_effects,
                         error_model, initialValue_ranges, nits,
                         out_folder, criteria, tasks = "all", normal_dist = NULL, fixed = NULL){
  
  if(any(tasks == "all")) tasks <- c("initial values", "SCM building", "COSSAC building", "SCM convergence", "COSSAC convergence")
  
  
  if("initial values" %in% tasks){
    #####
    # Identify region of population parameter space where the global minimum lies
    # Define parent project
    project_name <- "establish_initial_values.mlxtran"
    
    
    if(file.exists(here(out_folder, project_name))){
      loadProject(here(out_folder, project_name))
    }else{
      create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, normal_dist = normal_dist)
      
      saveProject(here(out_folder, project_name))
    }
    
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "log10_initial_value_assessment.RData")
    
    if(file.exists(file)){
      load(file)
    }else{
      
      # Run convergence over orders of magnitude
      initialValue_ranges_log10 <- initialValue_ranges %>% 
        t() %>% as.data.frame %>% 
        setNames(c("min", "max")) %>% 
        rownames_to_column("parameter") %>% 
        mutate(parameter = paste0(parameter, "_pop")) 
      
      log10_parameter_results <- convergence_assessment(initialValue_ranges_log10, log10 = T, nits = nits, fixed = fixed)
      
      for(name in names(log10_parameter_results)){
        assign(name, log10_parameter_results[[name]])
      }
      
      save(tabestimates, tabiters, pvals, phase, initialvalues,
           file = file)
      
      
    }
    log10_tabestimates <- tabestimates
    log10_initialvalues <- initialvalues
    
    plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = FALSE, remove_bicc_outliers = T, remove_estimate_outliers = T)
    
    
    #####
    # Define the ranges for a more focused parameter search
    # Taking ranges from estimates with best 25% of BICc
    initialValue_ranges_narrow <- get_narrow_ranges(tabestimates)
    
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "narrow_parameter_ranges.csv")
    if(!file.exists(file)) write_csv(initialValue_ranges_narrow, file)
    
    #####
    ## Identify global minimum for parameter estimates
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "initial_value_assessment.RData")
    if(file.exists(file)){
      load(file)
    }else{
      
      # Run convergence for narrower range of values
      parameter_results <- convergence_assessment(initialValue_ranges_narrow, log10 = F, nits = nits, fixed = fixed)
      
      for(name in names(parameter_results)){
        assign(name, parameter_results[[name]])
      }
      
      save(tabestimates, tabiters, pvals, phase, initialvalues, 
           file = file)
      
      correlate_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = T, remove_bicc_outliers = T)
      
      
      plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE, remove_bicc_outliers = T, remove_estimate_outliers = T)
      plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE, remove_bicc_outliers = T, remove_estimate_outliers = T,
                           file_addendum = "_nolog", log10_x = F)
      
      plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE, remove_bicc_outliers = T, remove_estimate_outliers = T, 
                           ylim = c(min(tabestimates$bicc), 1050), file_addendum = "_filtered")
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = T, 
            remove_estimate_outliers = T,
            log10_tabestimates, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
      
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots_outliers.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = F, 
            remove_estimate_outliers = F,
            log10_tabestimates, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
      
    }
    
    
    #####
    # Select best estimates as initial values
    
    best_parameter_values <- tabestimates %>% 
      filter(bicc == min(bicc))
    
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "best_initial_values.csv")
    if(!file.exists(file)) write_csv(best_parameter_values, file)
    
    #####
    # Generate report for final model with best initial values
    project_name <- "best_initial_values.mlxtran"
    
    if(!file.exists(here(out_folder, project_name))){
    saveProject(here(out_folder, project_name))
    set_base(best_parameter_values, fixed = fixed, omega_under_1 = F)
    runScenario()
    saveProject(here(out_folder, project_name))
    
    # Generate report with results
    rmarkdown::render(
      here('scripts', 'Monolix_template.Rmd'), 
      params = list(project = here(out_folder, project_name)), 
      output_file = here(out_folder, gsub(".mlxtran", "", project_name), "best_initial_values_report.html")
    )
    }
  }
  
  if("SCM building" %in% tasks){
    #####
    # Run SCM
    
    # Define parent project
    project_name <- "model_building_SCM.mlxtran"
    projectFile <- here(out_folder, project_name)
    
    if(file.exists(projectFile)){
      
      loadProject(projectFile)
      
    }else{
      
      create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, normal_dist = normal_dist)
      
      # Set initial values to best 
      set_base(best_parameter_values, fixed = fixed)
      
      saveProject(projectFile)
      
      covariate_model_building(projectFile, "SCM", criteria)
      
    }
  }
  
  if("SCM convergence" %in% tasks){
    
    #####
    # Convergence assessment for best model from SCM:
    project_name <- "convergence_assessment_SCM_final_model.mlxtran"
    
    if(file.exists(here(out_folder, project_name))){
      loadProject(here(out_folder, project_name))
    }else{
      # Define parent project based on prior results
      loadProject(here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran"))
      
      # Save as new project for convergence assessment
      
      saveProject(here(out_folder, project_name))
      
      # Set project settings and scenario
      project_settings()
    }
    
  
    # Convergence assessment over smaller ranges 
    initialValue_ranges_narrow <- read_csv(here(out_folder, "establish_initial_values", "narrow_parameter_ranges.csv"))
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "initial_value_assessment.RData")
    
    if(file.exists(file)){
      load(file)
    }else{
      parameter_results <- convergence_assessment(initialValue_ranges_narrow, log10 = F, nits = nits, fixed = fixed)
      
      for(name in names(parameter_results)){
        assign(name, parameter_results[[name]])
      }
      
      save(tabestimates, tabiters, pvals, phase, initialvalues, 
           file = file)
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = T, 
            remove_estimate_outliers = T,
            NULL, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
      
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots_outliers.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = F, 
            remove_estimate_outliers = F,
            NULL, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
    }
    
    
    
    plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE)
    plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE, log10_x = F, file_addendum = "_nolog", remove_estimate_outliers = T)
    
    correlate_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = T, remove_bicc_outliers = T)
    correlate_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = F, remove_bicc_outliers = F, file_addendum = "_outliers")
    
    compare_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = T, remove_bicc_outliers = T)
    
    
    #####
    # Select best estimates as initial values
    
    best_parameter_values <- tabestimates %>% 
      filter(bicc == min(bicc))
    
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "best_initial_values.csv")
    if(!file.exists(file)) write_csv(best_parameter_values, file)
    
    #####
    # Generate report for final model with best initial values
    project_name <- "SCM_final_model_best_initial_values.mlxtran"
    saveProject(here(out_folder, project_name))
    set_base(best_parameter_values, omega_under_1 = T, fixed = fixed)
    runScenario()
    saveProject(here(out_folder, project_name))
    
    # Generate report with results
    rmarkdown::render(
      here('scripts', 'Monolix_template.Rmd'), 
      params = list(project = here(out_folder, project_name)), 
      output_file = here(out_folder, gsub(".mlxtran", "", project_name), "SCM_best_model_report.html")
    )
  }
  
  if("COSSAC building" %in% tasks){
    #####
    # Run COSSAC
    
    # Define parent project
    project_name <- "model_building_COSSAC.mlxtran"
    projectFile <- here(out_folder, project_name)
    
    best_parameter_values <- read_csv(here(out_folder, "establish_initial_values", "best_initial_values.csv"))
    
    if(file.exists(projectFile)){
      loadProject(projectFile)
    }else{
      create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, normal_dist = normal_dist)
      
      # Set initial values to best 
      set_base(best_parameter_values, fixed = fixed)
      
      saveProject(projectFile)
      
      covariate_model_building(projectFile, "COSSAC", criteria)
    }
  }
  
  if("COSSAC convergence" %in% tasks){
    #####
    # Convergence assessment for best model from COSSAC:
    
    project_name <- "convergence_assessment_COSSAC_final_model.mlxtran"
    
    if(file.exists(here(out_folder, project_name))){
      loadProject(here(out_folder, project_name))
    }else{
      # Define parent project based on prior results
      loadProject(here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran"))
      
      # Save as new project for convergence assessment
      saveProject(here(out_folder, project_name))
      
      # Set project settings and scenario
      project_settings()
      
    }
    
    
    # Convergence assessment over smaller ranges 
    initialValue_ranges_narrow <- read_csv(here(out_folder, "establish_initial_values", "narrow_parameter_ranges.csv"))
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "initial_value_assessment.RData")
    
    if(file.exists(file)){
      load(file)
    }else{
      parameter_results <- convergence_assessment(initialValue_ranges_narrow, log10 = F, nits = nits, fixed = fixed)
      
      for(name in names(parameter_results)){
        assign(name, parameter_results[[name]])
      }
      
      save(tabestimates, tabiters, pvals, phase, initialvalues, 
           file = file)
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = T, 
            remove_estimate_outliers = T,
            NULL, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
      
      
      pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots_outliers.pdf"), width = 7, height = 8)
      for(parameter in getPopulationParameterInformation()$name){
        if(parameter %in% names(fixed)) next
        print(parameter)
        print(
          parameter_plots(
            parameter,
            # "lambda_pop",
            remove_bicc_outliers = F, 
            remove_estimate_outliers = F,
            NULL, tabestimates, tabiters,
            phase, initialvalues
          )
        )
      }
      dev.off()
    }
    
    
    
    plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE)
    plot_bicc_v_estimate(tabestimates, out_folder, project_name, narrow = TRUE, log10_x = F, file_addendum = "_nolog", remove_estimate_outliers = T)
    
    correlate_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = T, remove_bicc_outliers = T)
    correlate_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = F, remove_bicc_outliers = F, file_addendum = "_outliers")
    
    compare_estimates(tabestimates, out_folder, project_name, remove_estimate_outliers = T, remove_bicc_outliers = T)
    
    
    #####
    # Select best estimates as initial values
    
    best_parameter_values <- tabestimates %>% 
      filter(bicc == min(bicc))
    
    file <- here(out_folder, gsub(".mlxtran", "", project_name), "best_initial_values.csv")
    if(!file.exists(file)) write_csv(best_parameter_values, file)
    
    
    #####
    # Generate report for final model with best initial values
    project_name <- "COSSAC_final_model_best_initial_values.mlxtran"
    saveProject(here(out_folder, project_name))
    set_base(best_parameter_values, fixed = fixed)
    runScenario()
    saveProject(here(out_folder, project_name))
    
    # Generate report with results
    rmarkdown::render(
      here('scripts', 'Monolix_template.Rmd'), 
      params = list(project = here(out_folder, project_name)), 
      output_file = here(out_folder, gsub(".mlxtran", "", project_name), "COSSAC_best_model_report.html")
    )
  }
}
