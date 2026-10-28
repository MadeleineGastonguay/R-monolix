# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.

#' BICc weights and cumulative probability for each run
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
calculate_BICw <- function(tabestimates){
  # Calculate BICc weights 
  tabestimates %>% 
    distinct(run, bicc) %>% 
    mutate(
      t = exp(-(bicc - min(bicc))/2),
      w = t/sum(t)
    ) %>% 
    select(-t) %>% 
    arrange(desc(w)) %>% 
    mutate(cumprob = cumsum(w))
}

#' Narrower initial-value ranges from an initial guess sampling run
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
get_narrow_ranges <- function(tabestimates, p_threshold = NULL, bicc_cutoff = NULL){
  
  if(!is.null(p_threshold)){
    # First, get range of parameters that give "best" model
    bicw <- calculate_BICw(tabestimates)
    
    runs <- bicw %>% 
      filter(cumprob <= p_threshold) %>% 
      pull(run)
    
    bic_ranges <- tabestimates %>% 
      filter(run %in% runs) %>% 
      group_by(parameter) %>% 
      summarise(min = min(estimate), max = max(estimate), best = estimate[bicc == min(bicc)])
  }
  
  if(!is.null(bicc_cutoff)){
    bic_ranges <- tabestimates %>% 
      filter(bicc <= bicc_cutoff) %>% 
      group_by(parameter) %>% 
      summarise(min = min(estimate), max = max(estimate), best = estimate[bicc == min(bicc)])
  }
  
  # Now, calculate range as 1/2*best estimate to 2*best estimate, unless it is narrower than the above range
  search_ranges <- bic_ranges %>% 
    mutate(
      min = ifelse(best/2 > min, min, best/2),
      max = ifelse(best*2 < max, max, best*2)
    ) 
  
  return(search_ranges)
  
}

#' Flag values outside 1.5 IQR
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
outlier_detector <- function(value, upper_only = F){
  if(upper_only){
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) 
  }else{
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) |
      (value < quantile(value, probs = 0.25, na.rm = T) - 1.5*IQR(value, na.rm = T) )  
  }
  
}

#' Optionally drop outlier runs or estimates, returning a caption
#'
#' Ported from the SIV_HIV_rebound_kinetics functions.R.
#' @export
filter_results <- function(tabestimates, remove_estimate_outliers, remove_bicc_outliers){
  
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
    tabestimates <- tabestimates %>% group_by(parameter) %>% filter(!outlier_detector(estimate)) %>% ungroup()
    caption <- paste0(caption, "Estimate outliers removed")
  }
  
  # Remove runs with bicc outliers
  tabestimates <- tabestimates %>% filter(!run %in% remove_runs)
  
  return(list(tabestimates = tabestimates, caption = caption))
}

