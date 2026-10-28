# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.

#' Plot sampling and estimate ranges for each IGS step
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_parameter_ranges <- function(initialValue_ranges, log10_tabestimates, tabestimates = NULL, p_threshold = NULL, bicc_cutoff = NULL){
  iv_range <- initialValue_ranges %>% t() %>% as.data.frame %>% set_names(c("min", "max")) %>% rownames_to_column("parameter")
  log10_range <- log10_tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup
  narrow_iv_range <- get_narrow_ranges(log10_tabestimates, p_threshold, bicc_cutoff) %>% select(-best)
  df <- rbind(
    iv_range %>% mutate(step = "1", range = "Initial Guess", parameter = paste(parameter, "pop", sep = "_")),
    log10_range %>% mutate(step = "2", range = "Estimate Range"),
    narrow_iv_range %>% mutate(step = "2", range = "Initial Guess")
  ) %>% 
    filter(grepl("pop", parameter))
  if(!is.null(tabestimates)){
    narrow_range <- tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup
    best <- tabestimates %>% filter(bicc == min(bicc), grepl("pop", parameter)) %>% select(parameter, estimate)
    df <- rbind(df, 
                narrow_range %>% mutate(step = "3", range = "Estimate Range")
    ) %>% 
      filter(grepl("pop", parameter))
  }
  
  plot <- ggplot(df, aes(x= step, ymin = min, ymax = max, color = range)) + 
    geom_linerange(linewidth = 2) + 
    # geom_linerange(data = filter(df, space == "Narrow Range"), linewidth = 2) + 
    facet_wrap(~parameter, scales = "free_y") + 
    scale_y_log10() + 
    scale_color_manual(values = c("gray", "darkred"))
  
  if(!is.null(tabestimates)){
    plot <- plot + 
      geom_point(data = best, aes(x = "3", y = estimate), 
                 color = "darkgreen", inherit.aes = F, size = 2) +
      geom_text(data = best, aes(x = "3", y = estimate, label = signif(estimate, 2)), 
                color = "darkgreen", inherit.aes = F, fontface = 2, vjust = 0, hjust = 1)
  }
  
  plot
  
}

#' Plot sampling and estimate ranges for each IGS step (publication version)
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_parameter_ranges_v2 <- function(initialValue_ranges, log10_tabestimates, tabestimates = NULL, p_threshold = NULL, bicc_cutoff = NULL){
  iv_range <- initialValue_ranges %>% t() %>% as.data.frame %>% set_names(c("min", "max")) %>% rownames_to_column("parameter")
  iv_range <- iv_range %>% mutate(min = ifelse(parameter == "NPs", min*1e5, min), max = ifelse(parameter == "NPs", max*1e5, max))
  log10_tabestimates <- log10_tabestimates %>% mutate(estimate = ifelse(parameter == "NPs_pop", estimate*1e5, estimate)) 
  log10_range <- log10_tabestimates %>% group_by(parameter) %>%  summarise(min = min(estimate), max = max(estimate)) %>% ungroup
  narrow_iv_range <- get_narrow_ranges(log10_tabestimates, p_threshold, bicc_cutoff) %>% select(-best)
  tabestimates_cut <- calculate_BICw(tabestimates) %>% filter(cumprob <= p_threshold) %>% pull(run)
  df <- rbind(
    iv_range %>% mutate(step = "1", range = "Initial Guess", parameter = paste(parameter, "pop", sep = "_")),
    log10_range %>% mutate(step = "2", range = "Estimate Range"),
    narrow_iv_range %>% mutate(step = "2", range = "Global Minimum")
  ) %>% 
    filter(grepl("pop", parameter))
  if(!is.null(tabestimates)){
    tabestimates <- tabestimates  %>% mutate(estimate = ifelse(parameter == "NPs_pop", estimate*1e5, estimate))
    narrow_range <- tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup
    best_range <- tabestimates %>% filter(run %in% tabestimates_cut) %>%
      group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup 
      
    best <- tabestimates %>% filter(bicc == min(bicc), grepl("pop", parameter)) %>% select(parameter, estimate) %>% 
      mutate(parameter = factor(parameter, levels = c("lambda_pop", "betas_pop", "ta_pop", "m_pop", "p_pop", "NPs_pop"),
                                labels = c("lambda['pop']", "beta['pop']~' (x'~10^-7~')'", "t['a, pop']", "m['pop']", "p['pop']", "N['p, pop']")))
    df <- rbind(df, 
                narrow_range %>% mutate(step = "3", range = "Estimate Range"),
                best_range %>% mutate(step = "3", range = "Best Candidates")
    ) %>% 
      filter(grepl("pop", parameter))
  }
  
  plot <- df %>% mutate(range = factor(range, levels = c("Initial Guess", "Estimate Range", "Global Minimum", "Best Candidates")),
                        step = factor(step, levels = rev(c("1", "2", "3")), labels =rev( c("Wide range of\ninitial conditions", "IGS from wide range", "IGS from 'good' range"))),
                        parameter = factor(parameter, levels = c("lambda_pop", "betas_pop", "ta_pop", "m_pop", "p_pop", "NPs_pop"),
                                           labels = c("lambda['pop']", "beta['pop']~' (x'~10^-7~')'", "t['a, pop']", "m['pop']", "p['pop']", "N['p, pop']"))) %>% 
    # ggplot(aes(x= step, ymin = log10(min), ymax = log10(max), color = range)) + 
    ggplot(aes(y= step, xmin = log10(min), xmax = log10(max), color = range)) + 
    geom_linerange(linewidth = 2) + 
    # geom_linerange(data = filter(df, space == "Narrow Range"), linewidth = 2) + 
    facet_wrap(~parameter, scales = "free_x", nrow = 1, labeller = "label_parsed", strip.position = "bottom") + 
    # scale_y_log10() + 
    scale_color_manual(values = c("gray39", "gray", "darkred", "darkgreen")) + 
    scale_x_continuous(labels = scales::label_math(10^.x))
  
  # if(!is.null(tabestimates)){
  #   plot <- plot + 
  #     geom_point(data = best, aes(y = "IGS from 'good' range", x = log10(estimate)), 
  #                color = "black", inherit.aes = F, size = 2) +
  #     geom_text(data = best, aes(y = "IGS from 'good' range", x = log10(estimate), label = signif(estimate, 2)), 
  #               color = "black", inherit.aes = F, fontface = 2, vjust = 0, hjust = 1.1)
  # }
  # 
  plot  + theme(strip.placement = "outside", axis.title = element_blank())
  
}

