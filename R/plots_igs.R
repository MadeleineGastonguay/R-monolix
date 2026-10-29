# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Changes since porting: parameter order comes from parameter_levels() instead of the loaded
# Monolix project, so plots can be made from saved results; error parameters found with
# is_error_param() (a1, b2 ... as well as a, b, c); estimate_boxplots() names the CV row in
# the subtitle (the label overlapped the first value) and angles the parameter names;
# plot_vs_rses() puts RSE = 100 in the >=100 group (it fell in none), leaves out estimates
# without a standard error (with a caption), uses darker oranges and the usual parameter order.

#' Plot BICc against each parameter estimate
#' @param tabestimates Estimates of an initial guess sampling run
#'   (`tabestimates` from [initial_guess_sampling()]): one row per run and
#'   parameter, with `parameter`, `run`, `estimate`, `rse` and `bicc`.
#' @param out_folder,project_name Where to save the figure:
#'   `file.path(out_folder, project_name)` (a `.mlxtran` ending is dropped).
#'   With `out_folder = ""` (default) nothing is saved; the plot is only
#'   returned.
#' @param narrow `TRUE` marks the best run with a triangle; `FALSE` draws the
#'   wide-search version, with a dashed line at `BICc_threshold`.
#' @param remove_bicc_outliers Leave out runs whose BICc is an upper outlier
#'   ([outlier_detector()]); the caption says so.
#' @param remove_estimate_outliers Leave out estimates that are outliers for
#'   their parameter ([outlier_detector()]); the caption says so.
#' @param log10_x Log10 x axis (not for treatment effects).
#' @param exclude_treatment_effects Leave out the `beta_` parameters.
#' @param ylim `c(min, max)` of the BICc axis; runs with a higher BICc are left
#'   out, and the caption says so.
#' @param file_addendum Text added to the end of the file name.
#' @param parameter One parameter to plot, instead of all of them.
#' @param BICc_threshold With `narrow = FALSE`, the BICc of the dashed line
#'   (e.g. the highest BICc among the best-supported runs).
#' @param color `"BICc"` (default) colours points by BICc; a vector of run
#'   numbers instead highlights those runs.
#' @param color_label Legend title when `color` is a vector of runs.
#' @param point_size Size of the points.
#' @param add_best With `narrow = TRUE`, mark the best run with a triangle.
#' @param note Text added at the start of the caption, e.g. which runs are
#'   shown (see [igs_figures()]).
#' @return A ggplot.
#' @export
plot_bicc_v_estimate <- function(tabestimates, out_folder = "", project_name = "", narrow = TRUE,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE, log10_x = TRUE, 
                                 exclude_treatment_effects = FALSE, ylim = NULL,
                                 file_addendum = "",  parameter = NULL, BICc_threshold = NULL,
                                 color = "BICc", color_label = NULL, point_size = 1, add_best = TRUE, note = NULL){
  
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
      ggplot_report(aes(estimate, bicc, fill = color)) +
      geom_point(alpha = 0.75, size = point_size, shape = 21) 
      
      if(add_best) estimates_vs_bicc <- estimates_vs_bicc + geom_point(data = . %>% filter(bicc == min(bicc)), size = point_size+1, shape = 24)
  }else{
    estimates_vs_bicc <- tabestimates %>%
      # mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot_report(aes(estimate, bicc,  fill = color)) +
      geom_point(alpha = 0.75, shape = 21, size =  point_size + 1) +
      geom_hline(yintercept = BICc_threshold, lty = "dashed", alpha = 0.5)
  }
  caption <- join_caption(note,
                          if (!narrow && !is.null(BICc_threshold)) "Dashed line: BICc at or below which the top runs fall.",
                          if (narrow && add_best) "Triangle: best run (lowest BICc).",
                          if (any(color == "BICc")) "Colour: BICc.", caption)
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
  if(!is.null(ylim)) estimates_vs_bicc <- estimates_vs_bicc + ylim(ylim) + labs(caption = join_caption(caption, paste0("BICc greater than ", max(ylim), " removed")))
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
#' @inheritParams plot_bicc_v_estimate
#' @param initialvalues Starting values of the same runs (`initialvalues` from
#'   [initial_guess_sampling()]): `parameter`, `run`, `initial_value`.
#' @param parameter_set Which parameters to show: any of `"pop"`, `"omega"`,
#'   `"beta"` (treatment effects) and `"error"`.
#' @param save Save the figure to `out_folder`?
#' @return A ggplot.
#' @export
plot_distributions <- function(initialvalues, tabestimates, out_folder = "", project_name = "", parameter_set = "pop",
                               remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                               save = T, parameter = NULL, note = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, paste("Violins and boxes: distribution across runs; points: runs, coloured by BICc;",
                                      "triangle: best run (lowest BICc)."), filtered_results$caption)
  
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
    ggplot_report(aes(name, value)) + geom_violin(scale = "width") + 
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
#' @inheritParams plot_distributions
#' @return A grob with one panel per parameter (draw it with
#'   [grid::grid.draw()] or save it with [ggplot2::ggsave()]), or, with
#'   `parameter`, a ggplot with the Pearson correlation.
#' @export
plot_iv_vs_estimates <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                                 remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                                 parameter = NULL, note = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, "Each point is a run; triangle: best run (lowest BICc). Colour: BICc.", filtered_results$caption)
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  df <- merge(initialvalues %>% select(param = parameter, run, initial_value), tabestimates) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    arrange(param)
  
  params <- as.list(unique(df$param))
  
  if(is.null(parameter)){
    plots <- lapply(params, function(p){
      df %>% filter(param == p) %>% 
        ggplot_report(aes(initial_value, estimate, fill = bicc, shape = bicc == min(bicc))) + 
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
    
    plot <- arrangeGrob(patchworkGrob(plot), left = "Final Estimate", bottom = "Initial Value")  # built, not drawn
  }else{
    color_limit <- range(df$bicc)
    
    plot <- df %>% 
      filter(param == parameter) %>% 
      ggplot_report(aes(initial_value, estimate, fill = bicc, shape = bicc == min(bicc))) + 
      geom_point(size = 2, alpha = 0.75) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) + 
      coord_obs_pred() +
      labs(color = "BICc", y = "Final Estimate", x = "Initial Value", shape = "Best fit", caption = caption) + 
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
#' @inheritParams plot_distributions
#' @return A ggplot.
#' @export
cor_iv_estimates <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                             remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, note = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- join_caption(note, paste("Pearson correlation across runs between each parameter's initial value and its final",
                                      "estimate; triangles: p <= 0.05. Parameters that start at the same value in every run",
                                      "are left out."), filtered_results$caption)
  
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
    ggplot_report(aes(cor, fct_rev(parameter), fill = log10(pval), shape = pval <= 0.05)) + 
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
#' @inheritParams plot_bicc_v_estimate
#' @param method Correlation method, passed to [stats::cor()].
#' @param label_threshold Correlations at least this large (in absolute value)
#'   are written on the tiles.
#' @param text_color Colour of those labels.
#' @return A ggplot.
#' @export
correlate_estimates <- function(tabestimates, method = "spearman", out_folder = "", project_name = "",
                                remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE,  file_addendum = "",
                                label_threshold = 0.5, text_color = "black", note = NULL){
  
  # Filter tabestimates to remove estimate outliers and/or BICc outliers
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, paste0(str_to_title(method), " correlation of the final estimates across runs; labelled where |r| >= ",
                                       label_threshold, "."), filtered_results$caption)
  
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
    ggplot_report(aes(parameter1, parameter2, fill = value)) + 
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
#' @inheritParams plot_bicc_v_estimate
#' @param pop_only Only the population (`_pop`) parameters?
#' @return A GGally `ggmatrix`.
#' @export
compare_estimates <- function(tabestimates, out_folder = "", project_name = "", pop_only = TRUE,
                              remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, file_addendum = "",
                              note = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, paste("Each point is a run; triangle: best run (lowest BICc). Colour: BICc.",
                                      "Diagonal: distribution of each estimate; upper panels: correlations."),
                          filtered_results$caption)
  
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
    GGally::ggpairs(columns = 3:ncol(.), legend = c(2, 1), 
                    lower = list(continuous = function(data, mapping){
                      ggplot_report(data, mapping) + 
                        geom_point(aes(fill = bicc, shape = bicc == min(bicc)), alpha = 0.75) + 
                        geom_point(data = . %>% filter(bicc == min(bicc)), aes(fill = bicc, shape = bicc == min(bicc)), size = 2) + 
                        scale_fill_viridis_c(direction = -1) + 
                        scale_shape_manual(values = c(21, 24)) + 
                        labs(fill = "BICc", shape = "Best fit") 
                    } 
                    )
    ) + 
    theme_report() +
    labs(caption = caption) 
  
  if(out_folder != ""){
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), paste0("paired_SAEM_population_estimates", file_addendum, ".pdf")), plot, width = 7, height = 7)  
  }
  
  return(plot)
  
}

#' Coefficient of variation of initial values and final estimates
#' @inheritParams plot_distributions
#' @return A ggplot.
#' @export
calc_CV <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                    remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                    parameter = NULL, note = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, "CV = SD / mean across the runs shown, of the initial values and of the final estimates.",
                          filtered_results$caption)
  
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
    ggplot_report(aes(factor(param), abs(value), fill = name)) + 
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
#' @inheritParams plot_bicc_v_estimate
#' @param tabiters SAEM iterations of every run (`tabiters` from
#'   [initial_guess_sampling()]): one column per parameter, plus `iteration`
#'   and `run`.
#' @param phase Iteration at which each run switched from the exploratory to
#'   the smoothing phase (`phase` from [initial_guess_sampling()]).
#' @param filter Show only iterations from this one on, counted from the
#'   switch (e.g. `-200`).
#' @return A patchwork (all parameters) or a ggplot (one parameter).
#' @export
plot_SAEM <- function(tabiters, phase, tabestimates, out_folder = "", project_name = "",
                      remove_bicc_outliers = FALSE, filter = NULL, parameter = NULL, note = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- join_caption(note, paste("Each line is one run's SAEM iterations, shifted so that 0 (dashed line) is the",
                                      "switch from the exploratory to the smoothing phase. Solid line: final estimate",
                                      "of the best run (lowest BICc). Colour: BICc."), filtered_results$caption)
  xlab <- "Iteration (0 = switch to the smoothing phase)"
  
  popparams <- parameter_info(tabestimates)  # order and fixed status without a loaded project
  
  
  if(is.null(parameter)){
    lims <- range(tabestimates$bicc)
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
        ggplot_report( aes(x = iteration-phase, y = !!sym(param))) +
        geom_line(aes(group = run, color = bicc)) +
        theme(plot.title = element_text(hjust = .5)) +
        labs(title = param, x = xlab, y = NULL, color = "BICc") + 
        # scale_y_log10() + 
        # coord_cartesian(xlim = c(-200, 314)) + 
        geom_vline(xintercept = 0, lty = "dashed", color = "black") + 
        geom_hline(yintercept = tabestimates %>% filter(bicc == min(bicc), param == p) %>% pull(estimate),
                   alpha = 0.5) + 
        scale_color_viridis_c(direction = -1, limits = lims)
      i <- i + 1
    }
    
    plot <- wrap_plots(plotList, ncol = min(5, length(plotList))) +
      plot_layout(guides = "collect", axis_titles = "collect") +
      plot_annotation(caption = caption)
  }else{
    
    
    temp <- tabiters %>% 
      left_join(data.frame(run = unique(tabiters$run), phase = phase[,1])) %>% 
      filter(run %in% unique(tabestimates$run)) %>% 
      left_join(tabestimates %>% distinct(run, bicc)) 
    
    
    if(!is.null(filter)) temp <- temp %>% filter(iteration-phase >= filter)
    
    plot <- temp %>% 
      ggplot_report( aes(x = iteration-phase, y = !!sym(parameter))) +
      geom_line(aes(group = run, color = bicc)) +
      theme(plot.title = element_text(hjust = .5)) +
      labs(x = xlab, y = "Estimate", color = "BICc", caption = caption) + 
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
#' @inheritParams plot_bicc_v_estimate
#' @return A patchwork (all parameters) or a ggplot (one parameter).
#' @export
plot_estimates <- function(tabestimates, out_folder = "", project_name = "",
                           remove_bicc_outliers = FALSE, parameter = NULL, note = NULL){
  
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
  
  explanation <- paste("Points: each run's final estimate; error bars: +/- RSE. Solid line: estimate of the best run",
                       "(lowest BICc); dashed line: median across the runs shown. Colour: BICc.")
  if(is.null(parameter)){
    # plot population parameters, one panel each, with one BICc legend
    lims <- range(tabestimates$bicc)
    plotList <- list()
    for (param in popparams$name) {
      if (popparams[popparams$name == param, ]$method == "FIXED") next
      estimates <- tabestimates[tabestimates$param == param, ]
      plotList[[param]] <- ggplot_report(estimates, aes(x = run, y = estimate, color = bicc)) +
        geom_point() +
        geom_errorbar(aes(ymax = estimate * (1+as.numeric(rse)/100), 
                          ymin = estimate * (1-as.numeric(rse)/100))) +
        theme(plot.title = element_text(hjust = .5)) +
        labs(title = param, x = "Run", y = NULL, color = "BICc") +
        geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate),
                   alpha = 0.5) + 
        geom_hline(yintercept = estimates %>% pull(estimate) %>% median(),
                   alpha = 0.5, lty = "dashed") + 
        scale_color_viridis_c(direction = -1, limits = lims)
    }
    
    plot <- wrap_plots(plotList, ncol = min(5, length(plotList))) +
      plot_layout(guides = "collect", axis_titles = "collect") +
      plot_annotation(caption = join_caption(note, explanation, filtered_results$caption))
  }else{
    estimates <- tabestimates[tabestimates$param == parameter, ]
    plot <- ggplot_report(estimates %>% filter(run %in% unique(tabestimates$run)), 
                   aes(x = run, y = estimate, color = bicc)) +
      geom_point() +
      geom_errorbar(aes(ymax = estimate * (1+as.numeric(rse)/100), 
                        ymin = estimate * (1-as.numeric(rse)/100))) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs(x = "Run", y = "Final Estimate", caption = join_caption(note, explanation, filtered_results$caption)) +
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
#' @inheritParams plot_bicc_v_estimate
#' @return A ggplot.
#' @export
estimate_boxplots <- function(tabestimates, remove_bicc_outliers = F, remove_estimate_outliers = F, 
                              out_folder = "", project_name = "", file_addendum = NULL, note = NULL){
  
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
  
  # Monolix's order: _pop, omegas, correlations, covariate effects, then error-model parameters
  plot <- df %>% 
    mutate(parameter = factor(parameter, levels = parameter_levels(df))) %>% 
    ggplot_report(aes(factor(parameter), estimate)) + 
    geom_boxplot() + 
    geom_text(data = CV, aes(y = Inf, label = CV), vjust = 1) +
    scale_y_log10(breaks = breaks, limits = range(breaks)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Parameter", y = "Normalized Estimate",
         caption = join_caption(note, paste("Each estimate divided by its median across the runs shown;",
                                            "numbers at the top: CV (%) across those runs."), results_filtered$caption)) 
  
  
  if(out_folder != ""){
    file <- ifelse(is.null(file_addendum), "normalized_estimate_boxplots.pdf", str_interp("normalized_estimate_boxplots${file_addendum}.pdf"))
    ggsave(file.path(out_folder, gsub(".mlxtran", "", project_name), file), plot, width = 12, height = 7)
    
  }
  
  return(plot)
}

#' Combined diagnostic panel for one parameter
#' @inheritParams plot_SAEM
#' @inheritParams plot_distributions
#' @param parameter The parameter to show.
#' @param log10_tabestimates Estimates of the wide IGS, for its BICc panel;
#'   `NULL` to leave it out.
#' @param tabestimates Estimates of the narrow IGS.
#' @return A patchwork.
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
#' @inheritParams plot_bicc_v_estimate
#' @param y Column plotted against the RSE, `"bicc"` by default.
#' @return A ggplot.
#' @export
plot_vs_rses <- function(tabestimates, out_folder = "", project_name = "",  y = "bicc", file_addendum = NULL, note = NULL){
  names(tabestimates)[names(tabestimates) == y] <- "y"
  scales <- ifelse(y == "bicc", "free_x", "free")
  
  n_missing <- sum(!is.finite(tabestimates$rse))
  groups <- c("<50", "50-100", ">=100")
  plot <- tabestimates %>% filter(is.finite(rse)) %>% mutate(rse = abs(rse)) %>% 
    mutate(parameter = factor(parameter, levels = parameter_levels(.)),
           color = factor(case_when(rse < 50 ~ "<50", rse < 100 ~ "50-100", TRUE ~ ">=100"), levels = groups)) %>% 
    ggplot_report(aes(rse, y, color = color)) +  
    geom_point(size = 2) + 
    scale_x_log10() + 
    facet_wrap(~parameter, scales = scales) + 
    labs(x = "RSE (%)", y = if (y == "bicc") "BICc" else y, color = "RSE (%)",
         caption = join_caption(note, if (n_missing) paste0(n_missing, " estimate(s) without a standard error left out"))) +
    # the darker half of Oranges, so every group shows on white
    scale_color_manual(values = stats::setNames(c("#FDAE6B", "#E6550D", "#7F2704"), groups), drop = FALSE)
  
  if(!out_folder == ""){
    file <- ifelse(is.null(file_addendum), str_interp("rses_vs_${y}.pdf"), str_interp("rses_vs_${y}${file_addendum}.pdf"))
    ggsave(file.path(out_folder, str_remove(project_name, "\\.mlxtran$"), file), 
           plot, width = 10, height = 7)
  }
  
  plot
}

#' Viridis colour scale for BICc
#' @param ... Passed to [ggplot2::scale_color_viridis_c()].
#' @return A ggplot2 scale.
#' @export
color_bicc <- function(...) scale_color_viridis_c("BICc", direction = -1, ...)

#' Viridis fill scale for BICc
#' @param ... Passed to [ggplot2::scale_fill_viridis_c()].
#' @return A ggplot2 scale.
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
