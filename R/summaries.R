# Summaries for the checkpoints: identifiability after step 1 (which parameters to fix) and
# the paired comparison of candidate models in step 4 (which model to keep).

#' Identifiability summary of an initial guess sampling run
#'
#' One row per estimated population parameter, with the diagnostics used at
#' checkpoint 1 to decide which parameters to fix:
#' * spread of the estimate among the best-supported runs (cumulative BICc
#'   weight within `p_threshold`, as in [get_narrow_ranges()]);
#' * correlation between initial value and estimate across runs (an estimate
#'   that follows its start is not determined by the data);
#' * the strongest correlation with another parameter's estimate across runs;
#' * RSE of the best run and the median over runs;
#' * eta shrinkage (conditional mode) of the best run, on the `omega_` rows.
#'
#' The `flags` column lists which of the thresholds below are crossed. The
#' thresholds are rules of thumb to point at parameters worth a closer look
#' (profile likelihood, plots), not a decision.
#'
#' @param igs Output of [initial_guess_sampling()] (needs `tabestimates` and
#'   `initialvalues`; uses `shrinkage` if recorded).
#' @param p_threshold Cumulative BICc weight defining the best-supported runs.
#' @param method Correlation method.
#' @param thresholds Named list: `fold` (max/min among best runs), `cor_initial`,
#'   `cor_parameters` (absolute correlations), `rse` (%), `shrinkage` (%).
#' @return Data frame, one row per parameter.
#' @export
identifiability_summary <- function(igs, p_threshold = 0.95, method = "spearman",
                                    thresholds = list(fold = 2, cor_initial = 0.5, cor_parameters = 0.8,
                                                      rse = 50, shrinkage = 30)) {
  te <- igs$tabestimates
  iv <- igs$initialvalues
  th <- utils::modifyList(list(fold = 2, cor_initial = 0.5, cor_parameters = 0.8, rse = 50, shrinkage = 30), thresholds)

  # parameters fixed in every run aren't estimated
  fixed <- unique(iv$parameter[iv$method %in% "FIXED"])
  params <- setdiff(unique(te$parameter), fixed)
  te <- te[te$parameter %in% params, ]

  bicw <- calculate_BICw(te)
  top_runs <- bicw$run[bicw$cumprob <= p_threshold | seq_len(nrow(bicw)) == 1]
  best_run <- te$run[which.min(te$bicc)]

  # estimates as a run x parameter matrix
  wide <- tidyr::pivot_wider(te[, c("run", "parameter", "estimate")], names_from = parameter, values_from = estimate)
  est <- as.matrix(wide[, params, drop = FALSE])
  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = method))
  }

  rows <- lapply(params, function(p) {
    tp <- te[te$parameter == p, ]
    top <- tp$estimate[tp$run %in% top_runs]
    ivp <- iv[iv$parameter == p, c("run", "initial_value")]
    m <- merge(tp[, c("run", "estimate")], ivp, by = "run")
    others <- setdiff(params, p)
    cors <- vapply(others, function(q) safe_cor(est[, p], est[, q]), numeric(1))
    best_other <- if (length(cors) && any(is.finite(cors))) others[which.max(abs(cors))] else NA_character_
    data.frame(
      parameter = p,
      best = tp$estimate[tp$run == best_run][1],
      rse_best = suppressWarnings(as.numeric(tp$rse[tp$run == best_run][1])),
      median_rse = stats::median(suppressWarnings(as.numeric(tp$rse)), na.rm = TRUE),
      n_top_runs = length(top),
      top_min = min(top), top_max = max(top),
      top_fold = if (all(top > 0)) max(top) / min(top) else NA_real_,
      top_cv = if (length(top) > 1 && mean(top) != 0) stats::sd(top) / abs(mean(top)) else NA_real_,
      cor_initial = safe_cor(m$initial_value, m$estimate),
      max_abs_cor = if (is.na(best_other)) NA_real_ else abs(cors[[best_other]]),
      most_correlated_with = best_other,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  out$shrinkage <- NA_real_
  sh <- igs$shrinkage
  if (!is.null(sh) && "conditionalMode" %in% names(sh)) {
    sh <- sh[sh$run == best_run, ]
    for (k in seq_len(nrow(sh))) {
      nm <- sub("^eta_", "", sh$parameter[k])
      target <- intersect(c(nm, paste0("omega_", nm)), out$parameter)
      if (length(target)) out$shrinkage[out$parameter == target[length(target)]] <- sh$conditionalMode[k]
    }
  }

  flag <- function(cond, text) ifelse(!is.na(cond) & cond, text, NA_character_)
  flags <- cbind(
    flag(out$top_fold > th$fold, paste0("estimate spread > ", th$fold, "-fold among best runs")),
    flag(abs(out$cor_initial) > th$cor_initial, "estimate follows initial value"),
    flag(out$max_abs_cor > th$cor_parameters, paste0("correlated with ", out$most_correlated_with)),
    flag(out$rse_best > th$rse, paste0("RSE > ", th$rse, "%")),
    flag(out$shrinkage > th$shrinkage, paste0("shrinkage > ", th$shrinkage, "%"))
  )
  out$flags <- apply(flags, 1, function(f) paste(f[!is.na(f)], collapse = "; "))
  attr(out, "best_run") <- best_run
  attr(out, "top_runs") <- top_runs
  out
}

#' Paired comparison of candidate models
#'
#' Every candidate was fitted from the same sets of initial values (same run
#' numbers), so runs can be compared model against model. Only runs that
#' succeeded in every model are used.
#'
#' @param igs_list Named list of [initial_guess_sampling()] results, one per
#'   model.
#' @param delta BICc difference counted as a real difference.
#' @param fdr_method Multiple-testing adjustment for the Wald p-values
#'   ([stats::p.adjust()] method), applied over all models and runs together.
#' @return List of tables:
#'   * `summary`: per model, the number of paired runs, times best, times within
#'     `delta` of the best model of that run, lowest and median BICc, number of
#'     treatment effects and median condition number;
#'   * `pairwise`: for each ordered pair, the number of runs where `model` beats
#'     `other` by more than `delta`;
#'   * `per_run`: BICc of every model in every paired run;
#'   * `omegas`: random-effect SDs in each model's best run;
#'   * `wald`, `wald_summary`: adjusted Wald p-values of treatment effects, and
#'     per model and effect the share of runs significant at 0.05 and 0.01.
#' @export
compare_models <- function(igs_list, delta = 2, fdr_method = "fdr") {
  if (is.null(names(igs_list)) || any(!nzchar(names(igs_list)))) stop("igs_list must be a named list", call. = FALSE)
  models <- names(igs_list)
  te <- do.call(rbind, lapply(models, function(m) {
    x <- igs_list[[m]]$tabestimates
    x$model <- m
    x[, c("model", "run", "parameter", "estimate", "rse", "bicc")]
  }))
  te$model <- factor(te$model, levels = models)

  runs_by_model <- lapply(split(te$run, te$model), unique)
  paired <- Reduce(intersect, runs_by_model)
  dropped <- setdiff(Reduce(union, runs_by_model), paired)
  if (!length(paired)) stop("No run succeeded in every model", call. = FALSE)
  if (length(dropped)) message(length(dropped), " run(s) failed in at least one model and are left out: ",
                               paste(sort(dropped), collapse = ", "))

  per_run <- te[te$run %in% paired, ] %>%
    distinct(model, run, bicc) %>%
    group_by(run) %>%
    mutate(best_bicc = min(bicc), best_model = model[which.min(bicc)],
           is_best = model == best_model, within_delta = bicc - best_bicc <= delta) %>%
    ungroup() %>%
    arrange(run, model)

  cn <- do.call(rbind, lapply(models, function(m) {
    x <- igs_list[[m]]$condition_number
    if (is.null(x)) return(NULL)
    data.frame(model = m, run = x$run, condition_number = x$condition_number)
  }))

  summary <- per_run %>%
    group_by(model) %>%
    summarise(n_runs = n(), times_best = sum(is_best), times_within_delta = sum(within_delta),
              min_bicc = min(bicc), median_bicc = stats::median(bicc), .groups = "drop")
  summary$n_treatment_effects <- vapply(as.character(summary$model), function(m) {
    length(unique(te$parameter[te$model == m & grepl("^beta_", te$parameter)]))
  }, numeric(1))
  summary$median_condition_number <- vapply(as.character(summary$model), function(m) {
    if (is.null(cn)) return(NA_real_)
    v <- cn$condition_number[cn$model == m & cn$run %in% paired]
    if (!length(v) || all(is.na(v))) NA_real_ else stats::median(v, na.rm = TRUE)
  }, numeric(1))

  wide <- tidyr::pivot_wider(per_run[, c("run", "model", "bicc")], names_from = model, values_from = bicc)
  pairwise <- expand.grid(model = models, other = models, stringsAsFactors = FALSE)
  pairwise <- pairwise[pairwise$model != pairwise$other, ]
  pairwise$runs_better_by_delta <- mapply(function(a, b) sum(wide[[a]] < wide[[b]] - delta), pairwise$model, pairwise$other)
  rownames(pairwise) <- NULL

  best_runs <- te %>% filter(run %in% paired) %>% group_by(model) %>% filter(bicc == min(bicc)) %>%
    filter(run == min(run)) %>% ungroup()
  omegas <- best_runs %>% filter(grepl("^omega", parameter)) %>% select(model, parameter, estimate)

  wald <- do.call(rbind, lapply(models, function(m) {
    p <- igs_list[[m]]$pvals
    if (is.null(p) || !nrow(p)) return(NULL)
    pcol <- intersect(c("p.value", "pvalue", "p_value"), names(p))[1]
    if (is.na(pcol)) stop("No p-value column in the Wald tests of ", m, call. = FALSE)
    data.frame(model = m, run = p$run, parameter = p$parameter, p_value = p[[pcol]])
  }))
  wald_summary <- NULL
  if (!is.null(wald)) {
    wald <- wald[wald$run %in% paired, ]
    wald$p_adjusted <- stats::p.adjust(wald$p_value, method = fdr_method)
    wald_summary <- wald %>% group_by(model, parameter) %>%
      summarise(n_runs = n(), share_p05 = mean(p_adjusted < 0.05), share_p01 = mean(p_adjusted < 0.01),
                median_p_adjusted = stats::median(p_adjusted), .groups = "drop")
  }

  list(summary = as.data.frame(summary), pairwise = pairwise, per_run = as.data.frame(per_run),
       omegas = as.data.frame(omegas), wald = wald,
       wald_summary = if (is.null(wald_summary)) NULL else as.data.frame(wald_summary),
       paired_runs = sort(paired), dropped_runs = sort(dropped), delta = delta)
}
