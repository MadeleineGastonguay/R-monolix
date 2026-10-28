# The four-step workflow, one function per step, driven by analysis.yaml. The driver script
# (inst/scripts/run_pipeline.R) calls these in order. Each step:
#   - checks the config up to that step (stopping cleanly at a checkpoint decision),
#   - reuses its saved results when its run log matches the current config, or stops if the
#     settings changed, unless `rerun = TRUE` (which deletes and remakes the step's folder),
#   - saves its results (results.rds), diagnostics and run log in its own folder.
#
# Output layout under analysis.out_folder:
#   step1_igs/           wide.mlxtran, narrow.mlxtran (+ their Monolix folders), wide/ and
#                        narrow/ diagnostics, narrow_ranges.csv, identifiability_summary.csv,
#                        parameter_ranges.pdf, profile_likelihood/
#   step2_igs/           fixed.mlxtran, diagnostics/, best_estimates.csv
#   step3_model_building/search/      <method>.mlxtran per covariate search, results tables
#   step3_model_building/candidates/  <candidate>.mlxtran, candidates.yaml
#   step4_model_comparison/           initial_values.csv, <candidate>.mlxtran, comparison tables
#   final_model/         final_model.mlxtran and its report
#
# R seeds: wide IGS r_seed, narrow r_seed + 1, step 2 r_seed + 2, step 4 r_seed + 4, so the
# draws of different steps aren't the same uniform numbers.

#' Folders of each step's results
#'
#' @param config Output of [read_analysis_config()].
#' @return Named list of paths.
#' @export
pipeline_paths <- function(config) {
  root <- config$analysis$out_folder
  if (is_blank(root)) stop("analysis.out_folder is not set", call. = FALSE)
  list(root = root,
       step1 = file.path(root, "step1_igs"),
       step2 = file.path(root, "step2_igs"),
       step3_search = file.path(root, "step3_model_building", "search"),
       step3_candidates = file.path(root, "step3_model_building", "candidates"),
       step4 = file.path(root, "step4_model_comparison"),
       final = file.path(root, "final_model"))
}

#' Run the whole pipeline, stopping at the next checkpoint
#'
#' Runs steps 1 to 4 and the final model in order. At a checkpoint (a decision
#' still `~` in analysis.yaml) it prints what to decide and returns; record the
#' decision and call it again. Steps whose saved results match the config are
#' loaded, not refitted.
#'
#' @param config_path Path to analysis.yaml.
#' @param rerun Steps to refit even if saved results exist: any of `1`-`4` and
#'   `"final"`. Their folders are deleted first.
#' @return Invisibly, a list of the results of the steps that ran, with
#'   attribute `checkpoint` (the message) if it stopped at one.
#' @export
run_pipeline <- function(config_path = "analysis.yaml", rerun = NULL) {
  config <- read_analysis_config(config_path)
  out <- list()
  res <- tryCatch({
    out$step1 <- run_step1(config, rerun = 1 %in% rerun)
    out$step2 <- run_step2(config, rerun = 2 %in% rerun)
    out$step3 <- run_step3(config, rerun = 3 %in% rerun)
    out$step4 <- run_step4(config, rerun = 4 %in% rerun)
    out$final <- run_final_model(config, rerun = "final" %in% rerun)
    NULL
  }, pipeline_checkpoint = function(e) e)
  if (!is.null(res)) {
    message(conditionMessage(res))
    attr(out, "checkpoint") <- conditionMessage(res)
  }
  invisible(out)
}

# Reuse a step's saved results, or prepare its folder for a fresh run.
# TRUE = results can be loaded; FALSE = run the step. Stops if saved results don't match.
use_cached <- function(config, step, dir, rerun, before_checkpoint = FALSE) {
  if (rerun) {
    if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  } else if (check_run_log(config, step, dir, before_checkpoint)) {
    if (file.exists(file.path(dir, "results.rds"))) {
      message("Step ", step, ": loading saved results from ", dir)
      return(TRUE)
    }
  } else if (dir.exists(dir) && length(list.files(dir))) {
    # files but no run log: an interrupted run or someone else's files; don't overwrite silently
    stop(dir, " has files but no run log (an interrupted run?). Delete it or rerun the step.", call. = FALSE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  FALSE
}

save_step <- function(results, config, step, dir, before_checkpoint = FALSE) {
  saveRDS(results, file.path(dir, "results.rds"))
  write_run_log(config, step, dir, before_checkpoint = before_checkpoint)
  invisible(results)
}

# Save the loaded project to `file` (this also moves its results folder there)
save_project_as <- function(file) {
  lixoftConnectors::saveProject(normalizePath(file, mustWork = FALSE))
  invisible(file)
}

# Sampling-range table (rows = parameters) from get_narrow_ranges() output, _pop only
ranges_from_narrow <- function(narrow_ranges, log10) {
  nr <- narrow_ranges[grepl("_pop$", narrow_ranges$parameter), ]
  out <- data.frame(min = nr$min, max = nr$max, row.names = nr$parameter)
  attr(out, "log10") <- intersect(log10, nr$parameter)
  out
}

# initial_guess_sampling() with every record the pipeline uses
run_igs <- function(ranges, nits, seed, set_error, set_variance = NULL, set_treatment = NULL, fixed = NULL,
                    predefined = NULL) {
  initial_guess_sampling(ranges, log10 = attr(ranges, "log10"), nits = nits, seed = seed,
                         set_error = set_error, set_variance = set_variance, set_treatment = set_treatment,
                         fixed = fixed, predefined_initialvalues = predefined,
                         record_shrinkage = TRUE, record_pvals = TRUE, record_indiv_params = TRUE,
                         record_condition_number = TRUE)
}

diagnostics <- function(config, igs, dir) {
  pl <- config$step1_igs$plots
  save_igs_diagnostics(igs, dir, p_threshold = config$step1_igs$narrow$p_threshold,
                       remove_bicc_outliers = isTRUE(pl$remove_bicc_outliers),
                       remove_estimate_outliers = isTRUE(pl$remove_estimate_outliers))
}

#' Step 1: two-phase initial guess sampling
#'
#' Base model (no covariates, nothing fixed). A wide IGS over
#' `step1_igs.initial_ranges` (omegas over `omega_range`), narrow ranges from
#' its best-supported runs, then a narrow IGS with omegas starting at
#' `narrow.omega_initial`. Writes diagnostics for both, the identifiability
#' summary and parameter-range plot, and the optional profile likelihoods.
#' Checkpoint 1 (`decisions.fixed`) follows.
#'
#' @param config Output of [read_analysis_config()].
#' @param rerun Refit even if saved results exist (deletes them first).
#' @return List with `wide`, `narrow` (IGS results), `narrow_ranges`,
#'   `identifiability` and `profiles`.
#' @export
run_step1 <- function(config, rerun = FALSE) {
  validate_config(config, 1)
  dir <- pipeline_paths(config)$step1
  if (use_cached(config, 1, dir, rerun)) return(readRDS(file.path(dir, "results.rds")))
  s1 <- config$step1_igs
  seed <- config$analysis$r_seed

  spec <- model_spec(config, use_fixed = FALSE)
  setup_project(spec)
  save_project_as(file.path(dir, "wide.mlxtran"))
  pp <- lixoftConnectors::getPopulationParameterInformation()$name
  ranges <- sampling_ranges(spec, pp)
  starts <- error_starts(spec, pp)

  message("Step 1: wide initial guess sampling (", s1$wide_nits, " runs)")
  wide <- run_igs(ranges, s1$wide_nits, seed, set_error = starts)
  save_project_as(file.path(dir, "wide.mlxtran"))
  diagnostics(config, wide, file.path(dir, "wide"))

  narrow_ranges <- get_narrow_ranges(wide$tabestimates, p_threshold = s1$narrow$p_threshold,
                                     min_fold_width = s1$narrow$min_fold_width)
  utils::write.csv(narrow_ranges, file.path(dir, "narrow_ranges.csv"), row.names = FALSE)

  message("Step 1: narrow initial guess sampling (", s1$narrow$nits, " runs)")
  save_project_as(file.path(dir, "narrow.mlxtran"))
  narrow <- run_igs(ranges_from_narrow(narrow_ranges, attr(ranges, "log10")), s1$narrow$nits, seed + 1,
                    set_error = starts, set_variance = s1$narrow$omega_initial)
  save_project_as(file.path(dir, "narrow.mlxtran"))
  diagnostics(config, narrow, file.path(dir, "narrow"))

  identifiability <- identifiability_summary(narrow, p_threshold = s1$narrow$p_threshold)
  utils::write.csv(identifiability, file.path(dir, "identifiability_summary.csv"), row.names = FALSE)
  tryCatch(ggsave(file.path(dir, "parameter_ranges.pdf"),
                  plot_parameter_ranges(ranges, wide$tabestimates, narrow_ranges, narrow$tabestimates,
                                        p_threshold = s1$narrow$p_threshold,
                                        report_transforms = config$model$report_transforms),
                  width = 12, height = 6),
           error = function(e) message("Could not write parameter_ranges.pdf: ", conditionMessage(e)))

  profiles <- list()
  if (length(s1$profile_likelihood)) {
    best <- best_run_estimates(narrow$tabestimates)[, c("parameter", "estimate")]
    pdir <- file.path(dir, "profile_likelihood")
    dir.create(pdir, showWarnings = FALSE)
    for (p in names(s1$profile_likelihood)) {
      pl <- s1$profile_likelihood[[p]]
      message("Step 1: profile likelihood for ", p, " (", pl$n, " values)")
      profiles[[p]] <- profile_likelihood(p, unlist(pl$range), best, pl$n)
      utils::write.csv(profiles[[p]], file.path(pdir, paste0(p, ".csv")), row.names = FALSE)
    }
  }

  save_step(list(wide = wide, narrow = narrow, narrow_ranges = narrow_ranges, identifiability = identifiability,
                 profiles = profiles, sampling_ranges = ranges), config, 1, dir)
}

#' Step 2: initial guess sampling with the fixed parameters
#'
#' Skipped when `decisions.fixed` is `{}`: the narrow IGS of step 1 is used
#' instead. Otherwise the base model with `decisions.fixed` applied, sampled
#' over step 1's narrow ranges.
#'
#' @inheritParams run_step1
#' @return The IGS result used by step 3 (step 2's, or step 1's narrow IGS),
#'   with attribute `skipped`.
#' @export
run_step2 <- function(config, rerun = FALSE) {
  validate_config(config, 2)
  step1 <- run_step1(config)
  if (!length(config$decisions$fixed)) {
    message("Step 2: skipped (decisions.fixed is {}); using step 1's narrow IGS")
    return(structure(step1$narrow, skipped = TRUE))
  }
  dir <- pipeline_paths(config)$step2
  if (use_cached(config, 2, dir, rerun)) return(readRDS(file.path(dir, "results.rds")))

  spec <- model_spec(config)
  setup_project(spec)
  save_project_as(file.path(dir, "fixed.mlxtran"))
  pp <- lixoftConnectors::getPopulationParameterInformation()$name
  ranges <- ranges_from_narrow(step1$narrow_ranges, attr(step1$sampling_ranges, "log10"))
  message("Step 2: initial guess sampling with ", paste(names(spec$fixed), collapse = ", "), " fixed (",
          config$step2_igs$nits, " runs)")
  igs <- run_igs(ranges, config$step2_igs$nits, config$analysis$r_seed + 2, set_error = error_starts(spec, pp),
                 set_variance = config$step1_igs$narrow$omega_initial, fixed = spec$fixed)
  save_project_as(file.path(dir, "fixed.mlxtran"))
  diagnostics(config, igs, file.path(dir, "diagnostics"))
  utils::write.csv(best_run_estimates(igs$tabestimates), file.path(dir, "best_estimates.csv"), row.names = FALSE)
  save_step(structure(igs, skipped = FALSE), config, 2, dir)
}

# Base model (fixed parameters applied) or a candidate, started from step 2's best run
setup_from_best <- function(config, treatment_effects, initial_values, file) {
  spec <- model_spec(config, treatment_effects = treatment_effects)
  setup_project(spec)
  set_base(data.frame(parameter = initial_values$parameter, estimate = initial_values$initial_value),
           omega_under_1 = TRUE, fixed = spec$fixed)
  save_project_as(file)
  spec
}

#' Step 3: covariate search and candidate models
#'
#' Starts from step 2's best run (`_pop` values rounded, omegas at
#' `omega_initial`, treatment effects at 0). Runs each covariate search in
#' `step3_model_building.method` (none, SCM, COSSAC or both) and saves every
#' tested model; then checkpoint 2 (`candidate_models`); then builds and fits
#' once each candidate from [resolve_candidates()].
#'
#' @inheritParams run_step1
#' @return List with `searches` (model-building results per method),
#'   `candidates` (covariate models) and `initial_values`.
#' @export
run_step3 <- function(config, rerun = FALSE) {
  validate_config(config, 3, decisions_needed = FALSE)
  base <- run_step2(config)
  paths <- pipeline_paths(config)
  s3 <- config$step3_model_building
  iv <- prepare_initial_values(best_run_estimates(base$tabestimates), sig_figs = s3$initial_value_rounding,
                               omega_value = s3$omega_initial, treatment_value = 0)

  # covariate search(es), before checkpoint 2
  sdir <- paths$step3_search
  if (use_cached(config, 3, sdir, rerun, before_checkpoint = TRUE)) {
    searches <- readRDS(file.path(sdir, "results.rds"))
  } else {
    searches <- list()
    for (m in model_building_methods(config)) {
      message("Step 3: ", m, " covariate search")
      file <- file.path(sdir, paste0(tolower(m), ".mlxtran"))
      setup_from_best(config, NULL, iv, file)
      prm <- unlist(s3$parameters)
      searches[[m]] <- covariate_model_building(
        normalizePath(file), method = m, criteria = s3$criterion, covariates = unlist(s3$covariates),
        linearization = isTRUE(s3$use_linearization), relationships = s3$relationships,
        parameters = if (identical(prm, "all")) NULL else prm,
        lrt_thresold = c(s3$lrt_thresholds$forward, s3$lrt_thresholds$backward),
        scenario = config$scenario, report_config = config)
      utils::write.csv(scm_results_table(searches[[m]]), file.path(sdir, paste0(tolower(m), "_tested_models.csv")),
                       row.names = FALSE)
      tryCatch(ggsave(file.path(sdir, paste0(tolower(m), "_model_building.pdf")),
                      plot_model_building(searches[[m]])$plot, width = 10, height = 7),
               error = function(e) message("Could not plot the ", m, " search: ", conditionMessage(e)))
    }
    save_step(searches, config, 3, sdir, before_checkpoint = TRUE)
  }

  # checkpoint 2
  validate_config(config, 3)
  cdir <- paths$step3_candidates
  if (use_cached(config, 3, cdir, rerun)) return(readRDS(file.path(cdir, "results.rds")))
  candidates <- resolve_candidates(config, searches)
  writeLines(as_yaml_text(list(candidates = unclass(candidates)[names(candidates)],
                               duplicates = attr(candidates, "duplicates"))),
             file.path(cdir, "candidates.yaml"))
  single <- list()
  for (nm in names(candidates)) {
    message("Step 3: building candidate ", nm, " (", format_covariate_model(candidates[[nm]]), ")")
    setup_from_best(config, candidates[[nm]], iv, file.path(cdir, paste0(nm, ".mlxtran")))
    single[[nm]] <- tryCatch({
      lixoftConnectors::runScenario()
      save_project_as(file.path(cdir, paste0(nm, ".mlxtran")))
      data.frame(model = nm, covariate_model = format_covariate_model(candidates[[nm]]),
                 bicc = read_likelihood(result_methods())[["bicc"]])
    }, error = function(e) {
      message("Candidate ", nm, " failed to fit from step 2's best values: ", conditionMessage(e))
      data.frame(model = nm, covariate_model = format_covariate_model(candidates[[nm]]), bicc = NA_real_)
    })
  }
  utils::write.csv(do.call(rbind, single), file.path(cdir, "single_fits.csv"), row.names = FALSE)
  save_step(list(searches = searches, candidates = candidates, initial_values = iv), config, 3, cdir)
}

#' Step 4: paired comparison of the candidate models
#'
#' Draws one set of `nits` initial values (`_pop` over step 1's narrow ranges,
#' omegas at `omega_initial`, treatment effects at `treatment_initial`, fixed
#' parameters fixed) and fits every candidate from the same values, then
#' compares them with [compare_models()]. Checkpoint 3 (`final_model`)
#' follows.
#'
#' @inheritParams run_step1
#' @return List with `igs` (per candidate), `comparison` and `initial_values`.
#' @export
run_step4 <- function(config, rerun = FALSE) {
  validate_config(config, 4, decisions_needed = FALSE)
  step1 <- run_step1(config)
  step3 <- run_step3(config)
  dir <- pipeline_paths(config)$step4
  cdir <- pipeline_paths(config)$step3_candidates
  if (use_cached(config, 4, dir, rerun)) return(readRDS(file.path(dir, "results.rds")))
  s4 <- config$step4_model_comparison
  models <- names(step3$candidates)

  # every candidate's parameters, so one draw covers them all
  pp <- list()
  for (nm in models) {
    lixoftConnectors::loadProject(normalizePath(file.path(cdir, paste0(nm, ".mlxtran"))))
    pp[[nm]] <- lixoftConnectors::getPopulationParameterInformation()$name
  }
  all_params <- unique(unlist(pp))
  spec <- model_spec(config)
  ranges <- ranges_from_narrow(step1$narrow_ranges, attr(step1$sampling_ranges, "log10"))
  iv <- sample_initial_values(ranges, all_params, n = s4$nits, fixed = spec$fixed,
                              set_error = error_starts(spec, all_params), set_variance = s4$omega_initial,
                              set_treatment = s4$treatment_initial, seed = config$analysis$r_seed + 4)
  utils::write.csv(iv, file.path(dir, "initial_values.csv"), row.names = FALSE)

  igs <- list()
  for (nm in models) {
    message("Step 4: fitting ", nm, " from the shared initial values (", s4$nits, " runs)")
    lixoftConnectors::loadProject(normalizePath(file.path(cdir, paste0(nm, ".mlxtran"))))
    file <- file.path(dir, paste0(nm, ".mlxtran"))
    save_project_as(file)
    igs[[nm]] <- run_igs(NULL, s4$nits, NULL, set_error = NULL, fixed = spec$fixed,
                         predefined = iv[iv$parameter %in% pp[[nm]], ])
    save_project_as(file)
    saveRDS(igs[[nm]], file.path(dir, paste0(nm, "_igs.rds")))
  }

  comparison <- compare_models(igs, fdr_method = s4$fdr_method)
  for (k in c("summary", "pairwise", "omegas", "wald_summary")) {
    if (!is.null(comparison[[k]])) utils::write.csv(comparison[[k]], file.path(dir, paste0(k, ".csv")), row.names = FALSE)
  }
  tryCatch({
    ggsave(file.path(dir, "model_comparison.pdf"), plot_model_comparison(comparison), width = 8, height = 5)
    pv <- plot_treatment_pvalues(comparison)
    if (!is.null(pv)) ggsave(file.path(dir, "treatment_pvalues.pdf"), pv, width = 10, height = 7)
  }, error = function(e) message("Could not plot the comparison: ", conditionMessage(e)))
  save_step(list(igs = igs, comparison = comparison, initial_values = iv), config, 4, dir)
}

#' The final model and its report
#'
#' After checkpoint 3: copies the chosen candidate's step 4 project (left at
#' its best run) to `final_model/final_model.mlxtran` and writes its report.
#'
#' @inheritParams run_step1
#' @return Path to the final project, invisibly.
#' @export
run_final_model <- function(config, rerun = FALSE) {
  step4 <- run_step4(config)
  validate_config(config, 4)
  final <- config$step4_model_comparison$final_model
  dir <- pipeline_paths(config)$final
  file <- file.path(dir, "final_model.mlxtran")
  if (!rerun && file.exists(file)) {
    message("Final model: ", file, " already exists (rerun to remake it)")
    return(invisible(file))
  }
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  message("Final model: ", final)
  lixoftConnectors::loadProject(normalizePath(file.path(pipeline_paths(config)$step4, paste0(final, ".mlxtran"))))
  save_project_as(file)
  lixoftConnectors::runScenario()
  lixoftConnectors::computeChartsData()
  save_project_as(file)
  writeLines(as_yaml_text(list(final_model = final, rationale = config$step4_model_comparison$final_rationale,
                               best_run = step4$comparison$summary[step4$comparison$summary$model == final, ])),
             file.path(dir, "decision.yaml"))
  generate_report(file, config = config, title = paste0(config$analysis$name, ": final model (", final, ")"))
  invisible(file)
}
