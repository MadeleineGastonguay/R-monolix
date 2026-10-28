# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Changes since porting: parameter order comes from parameter_levels() instead of the loaded
# Monolix project, so plots can be made from saved results; error parameters found with
# is_error_param() (a1, b2 ... as well as a, b, c).

#' Plot BICc against each parameter estimate
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_bicc_v_estimate <- function(tabestimates, out_folder = "", project_name = "", narrow = TRUE,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE, log10_x = TRUE, 
                                 exclude_treatment_effects = FALSE, ylim = NULL,
                                 file_addendum = "",  parameter = NULL, BICc_threshold = NULL,
                                 color = "BICc", color_label = NULL, point_size = 1, add_best = TRUE){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  if(!is.null(parameter)){
    tabestimates <- tabestimates %>% filter(param == parameter)
  }else if(exclude_treatment_effects){
    tabestimates <- tabestimates %>% filter(!grepl("beta_", param))
  }
  
  if(!is.null(ylim)){
    tabestimates <- tabestimates %>% filter(bicc <= ylim[2])
  }
  
  ## Generate plot
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  if(any(color == "BICc")){
    tabestimates <- tabestimates %>% mutate(color = bicc)
  }else if(length(color) > 1){
    tabestimates <- tabestimates %>% mutate(color = run %in% color)
  }
  
  if(!is.factor(tabestimates$param)) tabestimates <- tabestimates %>% mutate(param = factor(param, levels = popparams$name))
  
  if(narrow){
    estimates_vs_bicc <- tabestimates %>%
      ggplot(aes(estimate, bicc, fill = color)) +
      geom_point(alpha = 0.75, size = point_size, shape = 21) 
      
      if(add_best) estimates_vs_bicc <- estimates_vs_bicc + geom_point(data = . %>% filter(bicc == min(bicc)), size = point_size+1, shape = 24)
  }else{
    estimates_vs_bicc <- tabestimates %>%
      # mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot(aes(estimate, bicc,  fill = color)) +
      geom_point(alpha = 0.75, shape = 21, size =  point_size + 1) +
      geom_hline(yintercept = BICc_threshold, lty = "dashed", alpha = 0.5) + 
      labs(caption = paste0("Dashed line indicates BICc below which fall the most probable models\n", caption))
  }
  estimates_vs_bicc <- estimates_vs_bicc + 
    facet_wrap(~param, scales = "free_x") + 
    labs(x = "Final estimate", y = "BICc", shape = "Best fit", 
         title = "Comparing BIC and parameter estimates",
         caption = caption) 
    # scale_shape_manual(values = c(21, 24))
  
  if(any(color == "BICc")){
    estimates_vs_bicc <- estimates_vs_bicc + scale_fill_viridis_c(direction = -1) + labs(fill = "BICc")
  }else if(length(color) > 1){
    estimates_vs_bicc <- estimates_vs_bicc + 
      scale_fill_manual(values = c(`TRUE` = "darkred", `FALSE` = "gray")) + labs(fill = color_label) + 
      guides(fill = guide_legend(override.aes = list(color = c("gray", "darkred"))))
  }
  
  
  if(any(grepl("beta_", tabestimates$param))) estimates_vs_bicc <- estimates_vs_bicc + labs(subtitle = ifelse(!exclude_treatment_effects, "", "treatment effects excluded"))
  if(!is.null(ylim)) estimates_vs_bicc <- estimates_vs_bicc + ylim(ylim) + labs(caption = paste0(caption, "\nBICc greater than ", max(ylim), " removed"))
  if(!is.null(parameter)){
    if(!grepl("beta_", parameter) & log10_x) estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  # }else if(log10_x & exclude_treatment_effects){ 
  }else if(log10_x){ 
    estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  }
  
  
  if(out_folder != ""){
    file_name <- ifelse(narrow, 
                        ifelse(!exclude_treatment_effects, "BIC_vs_estimates_narrow_range_treatment_effects", "BIC_vs_estimates_narrow_range"),
                        ifelse(!exclude_treatment_effects, "BIC_vs_estimates_treatment_effects", "BIC_vs_estimates_log_range"))
    file_name <- paste0(file_name, file_addendum, ".pdf")
    
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), file_name), estimates_vs_bicc, width = 10, height = 7)
    
  }
  
  return(estimates_vs_bicc)
  
}

#' Compare initial value and final estimate distributions
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_distributions <- function(initialvalues, tabestimates, out_folder = "", project_name = "", parameter_set = "pop",
                               remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                               save = T, parameter = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  subtitles <- c("pop" = "Population parameters",
                 "omega" = "Standard deviations of the random effect",
                 "beta" = "Treatment effects",
                 "error" = "Error model parameters")
  subtitle <- paste(subtitles[parameter_set], collapse = " and ")
  
  if(is.null(parameter)){
    params <- unlist(lapply(parameter_set, function(z){ 
      p <- c()
      if(z == "error"){
        p <- c(p, unique(tabestimates$param)[is_error_param(unique(tabestimates$param))])
      }else{
        if(z == "beta") z <- paste0(z, "_")
        p <- c(p, unique(tabestimates$param)[grepl(z, unique(tabestimates$param))])
      }
      return(p)
    }))
    
    df <- merge(initialvalues %>% select(param = parameter, run, initial_value), tabestimates) %>% 
      filter(param %in% params)
  }else{
    df <-  merge(initialvalues %>% select(param = parameter, run, initial_value), tabestimates) %>% 
      filter(param == parameter) 
  }
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  plot <- df %>% 
    pivot_longer(c(initial_value, estimate)) %>% 
    mutate(name = factor(name, levels = c("initial_value", "estimate"), labels = c("Initial Value", "Final Estimate")),
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
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), file_name), plot, width = 10, height = 7)  
  }
  
  return(plot)
}

#' Plot final estimates against initial values
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_iv_vs_estimates <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                                 remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                                 parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  df <- merge(initialvalues %>% select(param = parameter, run, initial_value), tabestimates) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    arrange(param)
  
  params <- as.list(unique(df$param))
  
  if(is.null(parameter)){
    plots <- lapply(params, function(p){
      df %>% filter(param == p) %>% 
        ggplot(aes(initial_value, estimate, fill = bicc, shape = bicc == min(bicc))) + 
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
      ggplot(aes(initial_value, estimate, fill = bicc, shape = bicc == min(bicc))) + 
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
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), "estimates_vs_ivs.pdf"), plot, width = 12, height = 10)  
  }
  
  return(plot)
}

#' Correlation between initial values and final estimates
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
cor_iv_estimates <- function(tabestimates, initialvalues, out_folder, project_name,
                             remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  df <- merge(initialvalues %>% select(parameter, run, initial_value), tabestimates) %>% 
    mutate(parameter = factor(parameter, levels = popparams$name)) %>% 
    group_by(parameter) %>% 
    # parameters that start at the same value every run (error, fixed) have no correlation
    filter(sd(initial_value) > 0, sd(estimate) > 0) %>% 
    summarise(cor = cor.test(initial_value, estimate)$estimate, pval = cor.test(initial_value, estimate)$p.value)
  
  breaks <- round(range(log10(df$pval)))
  breaks <- seq(breaks[1], breaks[2], by = 2)
  
  plot <- df %>% 
    ggplot(aes(cor, fct_rev(parameter), fill = log10(pval), shape = pval <= 0.05)) + 
    geom_vline(xintercept = 0, lty = "dashed") + 
    geom_point(size = 3) + 
    xlim(c(-1,1)) + 
    labs(x = "Correlation", y = "Parameter", caption = caption, 
         title = "Correlation between initial values and final estimates") +
    scale_fill_viridis_c(option = "plasma", direction = -1, breaks = breaks) + 
    scale_shape_manual(values = c(21, 24))
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), "cor_iv_estimates.pdf"), plot, width = 7, height = 8)  
  }
  
  return(plot)
}

#' Heatmap of correlations between estimates across runs
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
correlate_estimates <- function(tabestimates, method = "spearman", out_folder = "", project_name = "",
                                remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE,  file_addendum = "",
                                label_threshold = 0.5, text_color = "black"){
  
  # Filter tabestimates to remove estimate outliers and/or BICc outliers
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  correlation_df <- tabestimates %>% 
    select(param, estimate, run)  %>% 
    pivot_wider(names_from = param, values_from = estimate) %>% 
    select(-run) %>% 
    cor(use = "complete.obs", method = method) %>% # only keep runs that have no estimate outliers (if estimate outliers were removed)
    as.data.frame %>% 
    rownames_to_column("parameter2") %>% 
    pivot_longer(!parameter2, names_to = "parameter1") 
  
  correlation_heatmap <- correlation_df %>% 
    mutate(parameter2 = factor(parameter2, levels = popparams$name),
           parameter1 = factor(parameter1, levels = popparams$name)) %>%
    ggplot(aes(parameter1, parameter2, fill = value)) + 
    geom_tile() + 
    geom_text(data = . %>% filter(abs(value) >= label_threshold, parameter1 != parameter2),
              aes(label = round(value, 2)), color = text_color) +
    scale_fill_gradient2(limits = c(-1,1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) + 
    labs(fill = "correlation", title = str_interp("${str_to_title(method)} Correlation of SAEM estimates"), caption = caption)
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), paste0("SAEM_correlations", file_addendum, ".pdf")), correlation_heatmap, width =10, height = 10)  
  }
  
  return(correlation_heatmap)
  
}

#' Pairs plot of population estimates across runs
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
compare_estimates <- function(tabestimates, out_folder = "", project_name = "", pop_only = TRUE,
                              remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, file_addendum = ""){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  if(pop_only){
    df <- tabestimates %>% 
      select(param, estimate, run, bicc)  %>% 
      filter(grepl("pop", param)) 
  }else{
    df <- tabestimates %>% 
      select(param, estimate, run, bicc)  
  }
  
  plot <- df %>% 
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
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), paste0("paired_SAEM_population_estimates", file_addendum, ".pdf")), plot, width = 7, height = 7)  
  }
  
  return(plot)
  
}

#' Coefficient of variation of initial values and final estimates
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
calc_CV <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                    remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                    parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  cv <- merge(
    tabestimates, 
    initialvalues %>% rename(param = parameter)
  )  %>% 
    group_by(param) %>% 
    summarise(`Final Estimate` = sd(estimate)/mean(estimate),
              `Initial Value` = sd(initial_value)/mean(initial_value)) 
  
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
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), "coefficient_of_variation.pdf"), cv, width = 8.5, height = 10)  
  }
  
  return(cv)
}

#' Plot SAEM trajectories for each run
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_SAEM <- function(tabiters, phase, tabestimates, out_folder = "", project_name = "",
                      remove_bicc_outliers = FALSE, filter = NULL, parameter = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  # caption <- filtered_results$caption
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  
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
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), "SAEM_trajectories.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

#' Plot estimates and RSEs for each run
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_estimates <- function(tabestimates, out_folder = "", project_name = "",
                           remove_bicc_outliers = FALSE, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  # Get plot y-limit, excluding unusually large rse:
  ylim_range <-  tabestimates %>% 
    group_by(param) %>%
    filter(!outlier_detector(rse)) %>% 
    mutate(upper = estimate * (1+as.numeric(rse)/100), lower = estimate * (1-as.numeric(rse)/100)) %>% 
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
        geom_errorbar(aes(ymax = estimate * (1+as.numeric(rse)/100), 
                          ymin = estimate * (1-as.numeric(rse)/100))) +
        theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
        labs(title = param, x = NULL, y = NULL) +
        geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate),
                   alpha = 0.5) + 
        geom_hline(yintercept = estimates %>% pull(estimate) %>% median(),
                   alpha = 0.5, lty = "dashed") + 
        scale_color_viridis_c(direction = -1) 
      i <- i + 1
    }
    
    plot <- grid.arrange(grobs = plotList, ncol = 5)
  }else{
    estimates <- tabestimates[tabestimates$param == parameter, ]
    plot <- ggplot(estimates %>% filter(run %in% unique(tabestimates$run)), 
                   aes(x = run, y = estimate, color = bicc)) +
      geom_point() +
      geom_errorbar(aes(ymax = estimate * (1+as.numeric(rse)/100), 
                        ymin = estimate * (1-as.numeric(rse)/100))) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs(x = NULL, y = "Final Estimate") +
      geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate), alpha = 0.5) + 
      geom_hline(yintercept = estimates %>% pull(estimate) %>% median(),
                 alpha = 0.5, lty = "dashed") +
      scale_color_viridis_c(direction = -1) +
      coord_cartesian(ylim = unname(unlist(ylim_range[parameter,])))
    
    if(grepl("beta_", parameter) | grepl("omega_", parameter)){
      ylim <- layer_scales(plot)$y$range$range
      if(min(ylim) > 0) ylim[1] <- 0
      if(min(ylim) < 0 & max(ylim) < 0) ylim[2] <- 0
      plot <- plot + coord_cartesian(ylim = ylim)
    }
  }
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), "Convergence_estimates.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

#' Boxplots of estimates normalized to their median
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
estimate_boxplots <- function(tabestimates, remove_bicc_outliers = F, remove_estimate_outliers = F, 
                              out_folder = "", project_name = "", file_addendum = NULL){
  
  results_filtered <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers )
  tabestimates <- results_filtered$tabestimates
  
  CV <- tabestimates %>% group_by(parameter) %>% 
    summarise(CV = signif(sd(estimate)/mean(estimate)*100, 2))
  
  df <- tabestimates %>% 
    group_by(parameter) %>% 
    mutate(estimate = abs(estimate/median(estimate)) )
  
  limits <- range(df$estimate) 
  breaks <- c(1,2,5,10)
  bound <- ceiling(max(1/limits[1], limits[2]))
  breaks <- breaks[breaks <= bound]
  breaks <- sort(round(unique(c(1/breaks, breaks, bound, 1/bound)), 2))
  
  parameters <- unique(df$parameter)
  
  plot <- df %>% 
    mutate(parameter = factor(parameter, 
                              levels = c(parameters[grepl("pop", parameters)],
                                         parameters[grepl("omega", parameters)],
                                         parameters[grepl("beta_", parameters)],
                                         parameters[ parameters == "a"]))) %>% 
    ggplot(aes(factor(parameter), estimate)) + 
    geom_boxplot() + 
    geom_text(data = CV, aes(y = Inf, label = CV), vjust = 1) +
    annotate("text", -Inf, Inf, label = "CV (%)", vjust = 1, hjust = 0) +
    scale_y_log10(breaks = breaks, limits = range(breaks)) + 
    labs(x = "Parameter", y = "Normalized Estimate", caption = results_filtered$caption) 
  
  
  if(out_folder != ""){
    file <- ifelse(is.null(file_addendum), "normalized_estimate_boxplots.pdf", str_interp("normalized_estimate_boxplots${file_addendum}.pdf"))
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), file), plot, width = 12, height = 7)
    
  }
  
  return(plot)
}

#' Combined diagnostic panel for one parameter
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
parameter_plots <- function(parameter, remove_bicc_outliers, remove_estimate_outliers,
                            log10_tabestimates, tabestimates, tabiters,
                            phase, initialvalues){
  
  
  tabestimates <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
  tabestimates <- tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% rename(parameter = param)
  
  if(!is.null(log10_tabestimates)){
    log10_tabestimates <- filter_results(log10_tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
    log10_tabestimates <- log10_tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% rename(parameter = param)
  }
  
  SAEM_traj <- plot_SAEM(
    tabiters, phase, tabestimates, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    parameter = parameter#, filter = -200
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  SAEM_convergence <- plot_estimates(
    tabestimates,
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    parameter = parameter
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    labs(x = "Run")
  
  narrow_range <- tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% 
    pull(estimate) %>% range
  
  
  if(!is.null(log10_tabestimates)){
    wide_plot <- plot_bicc_v_estimate(
      log10_tabestimates, narrow = FALSE, 
      # remove_bicc_outliers =  remove_bicc_outliers, 
      # remove_estimate_outliers = remove_estimate_outliers,
      parameter = parameter
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
    log10_x = F,  parameter = parameter
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
    parameter = parameter
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  # Get coefficients of variation
  cv <- merge(
    tabestimates, 
    initialvalues
  )  %>% 
    group_by(parameter) %>% summarise(`Final Estimate` = sd(estimate)/mean(estimate),
                                      `Initial Value` = sd(initial_value)/mean(initial_value)) %>% 
    pivot_longer(!parameter)
  
  distributions <- plot_distributions(
    initialvalues, tabestimates, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    parameter = parameter
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch", 
         x = "Set") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    geom_text(data = cv, aes(label = round(value, 2)*100, y = Inf), vjust = 1) + 
    annotate(geom = "text", label = "CV (%):", y = Inf, x = -Inf, vjust = 1, hjust = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  cv <- calc_CV(
    tabestimates, initialvalues, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers,
    parameter = parameter
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

#' Plot RSEs against BICc or another column
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_vs_rses <- function(tabestimates, out_folder = "", project_name = "",  y = "bicc", file_addendum = NULL){
  names(tabestimates)[names(tabestimates) == y] <- "y"
  scales <- ifelse(y == "bicc", "free_x", "free")
  
  plot <- tabestimates %>% mutate(rse = abs(rse)) %>% 
    arrange(rse) %>% 
    mutate(color = fct_inorder(case_when(rse < 50~"<50", 
                                         rse >= 50 & rse < 100~ "50-100",
                                         rse > 100~ ">100"))) %>% 
    ggplot(aes(abs(rse), y, color = color)) +  
    geom_point() + 
    scale_x_log10() + 
    facet_wrap(~parameter, scales = scales) + 
    labs(x = "RSE (%)", y = y, color = "RSE (%)") +
    scale_color_brewer(palette = "Oranges")
  
  if(!out_folder == ""){
    file <- ifelse(is.null(file_addendum), str_interp("rses_vs_${y}.pdf"), str_interp("rses_vs_${y}${file_addendum}.pdf"))
    ggsave(file.path(out_folder, str_remove(project_name, "\\.mlxtran$"), file), 
           plot, width = 10, height = 7)
  }
  
  plot
}

#' Viridis colour scale for BICc
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
color_bicc <- function(...) scale_color_viridis_c("BICc", direction = -1, ...)

#' Viridis fill scale for BICc
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
fill_bicc <- function(...) scale_fill_viridis_c("BICc", direction = -1, ...)


# Parameters in Monolix's order: population values, random-effect SDs, correlations, covariate
# effects, error model; first-appearance order within each group.
parameter_levels <- function(tab) {
  p <- unique(as.character(if ("param" %in% names(tab)) tab$param else tab$parameter))
  group <- ifelse(grepl("_pop$", p), 1, ifelse(grepl("^omega", p), 2, ifelse(grepl("^corr_", p), 3,
                  ifelse(grepl("^beta_", p), 4, 5))))
  p[order(group, seq_along(p))]
}

# Parameter names in Monolix's order with a `method` column: "FIXED" for parameters with the
# same estimate in every run and no standard error, as Monolix reports fixed parameters
parameter_info <- function(tab) {
  if (!"param" %in% names(tab)) tab$param <- tab$parameter
  nm <- parameter_levels(tab)
  fixed <- vapply(nm, function(p) {
    x <- tab[tab$param == p, ]
    length(unique(x$estimate)) == 1 && "rse" %in% names(x) && all(is.na(suppressWarnings(as.numeric(x$rse))))
  }, logical(1))
  data.frame(name = nm, method = ifelse(fixed, "FIXED", "MLE"), stringsAsFactors = FALSE)
}
