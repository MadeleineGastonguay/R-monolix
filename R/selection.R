# Ported from SIV_HIV_rebound_kinetics/scripts/functions.R (Dec 2024).
# Behaviour unchanged apart from: package template path, file.path() instead of here(),
# and lixoftConnectors:: prefixes.
# Changes since porting: get_narrow_ranges() widens only when min_fold_width is given (the SIV
# code always widened to [best/2, 2*best]) and always keeps the best run (decided 2026-10-28:
# the SIV rule returned no runs when the best run alone carried more than p_threshold).

#' BICc weights and cumulative probability for each run
#' @param tabestimates Estimates of an initial guess sampling run, with `run`
#'   and `bicc`.
#' @return One row per run, best first: `run`, `bicc`, its BICc weight `w`
#'   and the cumulative weight `cumprob`.
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
#' The range of each parameter is the min-max of its estimates over the best-supported runs:
#' runs whose cumulative BICc weight is within `p_threshold` (the best run is
#' always kept, even when its own weight is above the threshold), or runs with
#' BICc at or below `bicc_cutoff`.
#'
#' @param tabestimates `tabestimates` from [initial_guess_sampling()].
#' @param p_threshold Cumulative BICc weight (see [calculate_BICw()]).
#' @param bicc_cutoff Alternative to `p_threshold`: keep runs with BICc at or
#'   below this value.
#' @param min_fold_width `NULL` (default) uses the min-max range as is. A number
#'   `f > 1` widens each range to at least `[best / f, best * f]` around the best
#'   run's estimate (positive estimates only).
#' @return Data frame with `parameter`, `min`, `max`, `best`.
#' @export
get_narrow_ranges <- function(tabestimates, p_threshold = NULL, bicc_cutoff = NULL, min_fold_width = NULL){

  if(is.null(p_threshold) == is.null(bicc_cutoff)) stop("Give one of p_threshold or bicc_cutoff.", call. = FALSE)

  if(!is.null(p_threshold)){
    bicw <- calculate_BICw(tabestimates)
    runs <- bicw %>% filter(cumprob <= p_threshold | row_number() == 1) %>% pull(run)
  }else{
    runs <- tabestimates %>% filter(bicc <= bicc_cutoff) %>% distinct(run) %>% pull(run)
    if(!length(runs)) stop("No run has BICc at or below ", bicc_cutoff, call. = FALSE)
  }
  best_run <- tabestimates$run[which.min(tabestimates$bicc)]

  search_ranges <- tabestimates %>%
    filter(run %in% runs) %>%
    group_by(parameter) %>%
    summarise(min = min(estimate), max = max(estimate), best = estimate[run == best_run][1], .groups = "drop")

  if(!is.null(min_fold_width)){
    stopifnot(is.numeric(min_fold_width), min_fold_width > 1)
    widen <- search_ranges$best > 0
    search_ranges$min[widen] <- pmin(search_ranges$min[widen], search_ranges$best[widen] / min_fold_width)
    search_ranges$max[widen] <- pmax(search_ranges$max[widen], search_ranges$best[widen] * min_fold_width)
  }

  search_ranges
}

#' Flag values outside 1.5 IQR
#' @param value Numeric vector.
#' @param upper_only Flag only values above the upper limit.
#' @return Logical vector: `TRUE` for values more than 1.5 times the
#'   interquartile range above the 3rd quartile (or below the 1st).
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
#' @param tabestimates Estimates of an initial guess sampling run.
#' @param remove_estimate_outliers Leave out estimates that are outliers for
#'   their parameter ([outlier_detector()]).
#' @param remove_bicc_outliers Leave out runs whose BICc is an upper outlier.
#' @return List with the filtered `tabestimates` and a `caption` saying what
#'   was left out.
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


#' Best-supported runs of an initial guess sampling run
#'
#' The runs with the lowest BICc whose BICc weights add up to `p_threshold`
#' ([calculate_BICw()]); the best run is always included. These are the runs
#' [get_narrow_ranges()] uses, and the ones the consistency plots show.
#'
#' @param tabestimates `tabestimates` from [initial_guess_sampling()].
#' @param p_threshold Cumulative BICc weight.
#' @return Run numbers, best first.
#' @export
top_runs <- function(tabestimates, p_threshold = 0.95) {
  bicw <- calculate_BICw(tabestimates)
  bicw$run[bicw$cumprob <= p_threshold | seq_len(nrow(bicw)) == 1]
}

# Caption line saying which runs a plot shows
runs_note <- function(tabestimates, runs = NULL, p_threshold = 0.95, sets_ranges = FALSE) {
  n_all <- length(unique(tabestimates$run))
  if (is.null(runs)) return(paste0("All ", n_all, " runs."))
  paste0("Top ", length(runs), " of ", n_all, " runs: the lowest-BICc runs making up ", signif(100 * p_threshold, 3),
         "% of the BICc weight", if (sets_ranges) " (these set the narrow ranges)" else "", ".")
}

# Non-empty caption pieces, one per line, wrapped so they fit narrow figures
join_caption <- function(..., width = 100) {
  x <- as.character(unlist(list(...)))
  x <- unlist(strsplit(x[!is.na(x)], "\n", fixed = TRUE))
  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x)) paste(unlist(lapply(x, strwrap, width = width)), collapse = "\n") else NULL
}
