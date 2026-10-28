# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.

#' Plot and summarise the steps of a stepwise covariate search
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
plot_model_building <- function(mb,  ...){
  mb_summary <- mb %>% 
    map_df(function(x) data.frame(
      param = x$tested[1], 
      treatment = x$tested[2],
      BICc = x$BICc, 
      LL = x$LL, 
      n_covar = sum(x$individualModels))
    )
  
  
  max_covar <- max(mb_summary$n_covar)
  mb_summary_no_back_elim <- mb_summary %>% 
    mutate(id = 1:nrow(.), t = c(0,diff(n_covar))) %>% 
    # filter(t >= 0) 
    filter(id <= max(id[n_covar == max_covar]))
  
  
  # best_mod_covar <- filter(mb_summary_no_back_elim, BICc == min(BICc))$n_covar
  
  mb_summary_no_back_elim <- mb_summary_no_back_elim %>% 
    # filter(id <= max(id[n_covar == best_mod_covar])) %>% 
    select(-id, -t)
  
  thresholds <- mb_summary_no_back_elim %>% 
    group_by(n_covar) %>% 
    summarise(threshold = BICc[BICc == min(BICc)]) %>% 
    ungroup %>% 
    mutate(n_covar = n_covar + 1)
  
  df <- mb_summary_no_back_elim %>% 
    filter(n_covar != 0) %>% 
    left_join(thresholds)
  
  p <- df %>% 
    filter(param != "") %>% 
    group_by(n_covar) %>% mutate(best = BICc == min(BICc)) %>% 
    ungroup %>% group_by(param, treatment) %>% mutate(final_mod = any(best)) %>% 
    ggplot(aes(BICc, treatment, color = param)) + 
    geom_rect(data = . %>% filter(best), 
              aes(xmin = BICc, xmax = BICc + 2, ymin = -Inf, ymax = Inf), 
              fill = "grey", alpha = 0.25, color = NA) +
    # geom_rect(data = . %>% filter(best), 
    #           aes(xmin = BICc, xmax = BICc + 5, ymin = -Inf, ymax = Inf), 
    #           fill = "grey", alpha = 0.15, color = NA) +
    geom_vline(aes(xintercept = threshold), alpha = 0.5) +
    geom_point(size = 2,  aes(shape = interaction(best, final_mod)), 
               position  = position_jitter(height = 0.2, seed = 100, width = 0)) +
    # geom_point(data = . %>% filter(final_mod), size = 2,  aes(shape = best)) +
    ggrepel::geom_text_repel(
      data = . %>% group_by(n_covar) %>%  mutate(label = ifelse(BICc <= (min(BICc) + 2), param, "")), aes(label = label), 
      hjust = 0, min.segment.length = 0, direction = "y",
      show.legend = F, position = position_jitter(height = 0.2, seed = 100, width = 0)
    ) +
    facet_wrap(~n_covar, labeller = "label_both", dir = "v", ...) + 
    ggsci::scale_color_lancet() +
    scale_shape_manual(values = c(1, 19, 17), labels = c("Not in final model", "Included in final model", "Best covariate at this step")) +
    # guides(shape = "none") + 
    labs(color = "Parameter", shape = "Treatment effect category", caption = "Forward selection only\nClosed shapes are effects included in the final model") + 
    theme(panel.border = element_rect(fill = NA)) 
  
  # Add custom breaks while preserving the default range
  return(list(plot = p, summary = left_join(mb_summary, thresholds)))
  
}

