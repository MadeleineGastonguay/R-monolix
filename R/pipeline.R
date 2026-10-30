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
#                        parameter_ranges.pdf, profile_likelihood/ (run_profiles(): profile.mlxtran,
#                        <parameter>.csv, profile_summary.csv, profile_likelihood.pdf, run log),
#                        igs_figures.html (every figure on one page)
#   step2_igs/           fixed.mlxtran, diagnostics/, best_estimates.csv, identifiability_summary.csv,
#                        igs_figures.html
#   step3_model_building/             model_building.html (the search and the candidate fits)
#   step3_model_building/search/      <method>.mlxtran per covariate search, results tables
#   step3_model_building/candidates/  <candidate>.mlxtran (scm_best and cossac_best copied from
#                                     the search's final model), candidates.yaml, single_fits.csv
#   step4_model_comparison/           initial_values.csv, <candidate>.mlxtran, comparison tables,
#                                     comparison.html
#   final_model/         final_model.mlxtran and its report, decision.yaml, paired_SAEM_population_estimates.pdf,
#                        normalized_estimate_boxplots.pdf,
#                        rses_vs_bicc.pdf (from the final model's step 4 runs)
#   history.yaml         one entry per finished step (read_history())
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
#' Runs steps 1 to 4 and the final model in order, with the optional profile
#' likelihoods ([run_profiles()]) after step 1. At a checkpoint (a decision
#' still `~` in analysis.yaml) it prints what to decide and returns; record the
#' decision and call it again. Steps whose saved results match the config are
#' loaded, not refitted.
#'
#' @param config_path Path to analysis.yaml.
#' @param rerun Steps to refit even if saved results exist: any of `1`-`4`,
#'   `"profiles"` and `"final"`. Their folders are deleted first.
#' @return Invisibly, a list of the results of the steps that ran, with
#'   attribute `checkpoint` (the message) if it stopped at one.
#' @export
run_pipeline <- function(config_path = "analysis.yaml", rerun = NULL) {
  config <- read_analysis_config(config_path)
  out <- list()
  res <- tryCatch({
    out$step1 <- run_step1(config, rerun = 1 %in% rerun)
    out$profiles <- run_profiles(config, out$step1, rerun = "profiles" %in% rerun)
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
    changes <- run_log_changes(config, step, dir, before_checkpoint)
    if (!is.null(changes)) note_replaced(dir, changes)   # for history.yaml
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
  record_history(config, dir, step, results)
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

diagnostics <- function(config, igs, dir, sets_ranges = FALSE) {
  pl <- config$step1_igs$plots
  save_igs_diagnostics(igs, dir, p_threshold = config$step1_igs$narrow$p_threshold,
                       remove_bicc_outliers = isTRUE(pl$remove_bicc_outliers),
                       remove_estimate_outliers = isTRUE(pl$remove_estimate_outliers), sets_ranges = sets_ranges)
}

#' Step 1: two-phase initial guess sampling
#'
#' Base model (no covariates, nothing fixed). A wide IGS over
#' `step1_igs.initial_ranges` (omegas over `omega_range`), narrow ranges from
#' its best-supported runs, then a narrow IGS with omegas starting at
#' `narrow.omega_initial`. Writes diagnostics for both, the identifiability
#' summary and parameter-range plot. The optional profile likelihoods run
#' next, in [run_profiles()]; checkpoint 1 (`decisions.fixed`) follows.
#'
#' @param config Output of [read_analysis_config()].
#' @param rerun Refit even if saved results exist (deletes them first).
#' @return List with `wide`, `narrow` (IGS results), `narrow_ranges`,
#'   `identifiability` and `sampling_ranges`.
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
  diagnostics(config, wide, file.path(dir, "wide"), sets_ranges = TRUE)

  narrow_ranges <- get_narrow_ranges(wide$tabestimates, p_threshold = s1$narrow$p_threshold,
                                     min_fold_width = s1$narrow$min_fold_width, bounds = spec$bounds)
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

  out <- save_step(list(wide = wide, narrow = narrow, narrow_ranges = narrow_ranges, identifiability = identifiability,
                        sampling_ranges = ranges), config, 1, dir)
  figures_page(config, 1)
  out
}

#' Profile likelihoods after step 1
#'
#' Optional, between step 1 and checkpoint 1: for each parameter in
#' `step1_igs.profile_likelihood`, fits the base model once per value on a
#' grid with that parameter fixed. The grid is either `values` as listed, or
#' round numbers over `range`, at least `n` of them, on a log10 scale with
#' `log10: true` and a linear one otherwise ([profile_grid()]).
#'
#' Each fit starts from the best narrow IGS run of step 1: population
#' parameters at its estimates, omegas at 1, error-model parameters at their
#' usual starts ([error_starts()]). Standard errors aren't computed, since
#' only the OFV is used.
#'
#' Nothing in step 1 is changed or refitted. The fits run in a copy of the
#' narrow project, `step1_igs/profile_likelihood/profile.mlxtran`, which
#' each fit overwrites; results are the tables. Writes `<parameter>.csv` per
#' parameter, `profile_summary.csv` ([profile_summary()]),
#' `profile_likelihood.pdf` ([plot_profile_likelihood()]) and a run log.
#' A parameter whose grid values haven't changed is loaded, not refitted, so
#' adding a parameter or changing one's grid refits only that one.
#'
#' @inheritParams run_step1
#' @param step1 Output of [run_step1()]; `NULL` loads its saved results.
#' @param rerun Refit every profile (deletes the profile folder first).
#' @return Named list of [profile_likelihood()] tables (empty when none are
#'   configured), with attribute `summary`.
#' @export
run_profiles <- function(config, step1 = NULL, rerun = FALSE) {
  validate_config(config, 1)
  # the grid itself is what's compared with the saved profiles, so a profile is refitted only
  # when its values change
  settings <- lapply(config$step1_igs$profile_likelihood, function(pl) {
    log10 <- isTRUE(pl$log10)
    values <- if (!is.null(pl$values)) sort(unique(unlist(pl$values))) else profile_grid(unlist(pl$range), pl$n, log10)
    list(values = values, log10 = log10)
  })
  dir <- file.path(pipeline_paths(config)$step1, "profile_likelihood")
  if (rerun && dir.exists(dir)) unlink(dir, recursive = TRUE)
  if (!length(settings)) return(invisible(list()))
  if (is.null(step1)) step1 <- run_step1(config)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # settings of the profiles already on disk, from the last run log
  log_file <- file.path(dir, "run_log.yaml")
  done <- if (file.exists(log_file)) yaml::read_yaml(log_file, handlers = yaml_handlers())$extra$profile_likelihood else NULL
  if (is.null(done)) done <- list()
  done <- done[vapply(names(done), function(p) file.exists(file.path(dir, paste0(p, ".csv"))), logical(1))]
  same <- function(p) !is.null(done[[p]]) &&
    identical(flatten_settings(done[[p]]), flatten_settings(yaml_roundtrip(settings[[p]])))
  todo <- names(settings)[!vapply(names(settings), same, logical(1))]
  for (p in setdiff(names(settings), todo)) message("Profile likelihood for ", p, ": loading saved results")

  best <- best_run_estimates(step1$narrow$tabestimates)
  fitted <- character()
  if (length(todo)) {
    lixoftConnectors::loadProject(normalizePath(file.path(pipeline_paths(config)$step1, "narrow.mlxtran")))
    save_project_as(file.path(dir, "profile.mlxtran"))  # a copy, so the narrow project keeps its results
    apply_scenario(utils::modifyList(config$scenario, list(standard_errors = FALSE)))
    iv <- profile_initial_values(config, best, lixoftConnectors::getPopulationParameterInformation())
    done <- done[intersect(names(done), names(settings))]
    for (p in todo) {
      pl <- settings[[p]]
      message("Profile likelihood for ", p, " (", length(pl$values), " values: ", paste(pl$values, collapse = ", "), ")")
      pr <- profile_likelihood(p, range(pl$values), iv, length(pl$values), values = pl$values)
      if (is.null(pr) || !nrow(pr)) { warning("No profile likelihood fit for ", p, " succeeded", call. = FALSE); next }
      utils::write.csv(pr, file.path(dir, paste0(p, ".csv")), row.names = FALSE)
      done[[p]] <- pl
      fitted <- c(fitted, p)
      # after every parameter, so an interrupted run keeps the finished ones
      write_run_log(config, 1, dir, extra = list(profile_likelihood = done))
    }
    save_project_as(file.path(dir, "profile.mlxtran"))
  }

  have <- names(settings)[file.exists(file.path(dir, paste0(names(settings), ".csv")))]
  have <- intersect(have, names(done))
  profiles <- lapply(stats::setNames(have, have), function(p) utils::read.csv(file.path(dir, paste0(p, ".csv"))))
  if (!length(profiles)) return(invisible(profiles))
  summary <- profile_summary(profiles, n_values = vapply(settings[have], function(pl) as.numeric(length(pl$values)), 1))
  utils::write.csv(summary, file.path(dir, "profile_summary.csv"), row.names = FALSE)
  if (length(fitted)) {
    record_history(config, dir, "profiles",
                   decisions = list(profile_likelihood = lapply(settings[fitted], function(pl) list(values = pl$values))),
                   result = lapply(split(summary[summary$profiled_parameter %in% fitted, , drop = FALSE],
                                         summary$profiled_parameter[summary$profiled_parameter %in% fitted]),
                                   function(r) as.list(r[1, setdiff(names(r), "profiled_parameter"), drop = FALSE])))
    figures_page(config, 1)   # adds the profiles to step 1's page
  }
  log10 <- names(settings)[vapply(settings, function(pl) pl$log10, logical(1))]
  tryCatch(ggsave(file.path(dir, "profile_likelihood.pdf"),
                  plot_profile_likelihood(profiles, log10 = log10, best = stats::setNames(best$estimate, best$parameter)),
                  width = 5 * min(length(profiles), 3), height = 4 * ceiling(length(profiles) / 3)),
           error = function(e) message("Could not write profile_likelihood.pdf: ", conditionMessage(e)))
  attr(profiles, "summary") <- summary
  invisible(profiles)
}

# Starting values for the profile fits: population parameters at the best run's estimates,
# omegas at 1, error-model parameters at their usual starts. `info` is the project's
# getPopulationParameterInformation(); anything else keeps the project's value.
profile_initial_values <- function(config, best, info) {
  pp <- info$name
  v <- stats::setNames(best$estimate, best$parameter)[pp]
  v[is.na(v)] <- info$initialValue[is.na(v)]
  v[grepl("^omega", pp)] <- 1
  err <- error_starts(model_spec(config, use_fixed = FALSE), pp)
  v[names(err)] <- err
  data.frame(name = pp, initialValue = unname(v), method = "MLE", stringsAsFactors = FALSE)
}

#' Step 2: initial guess sampling with the checkpoint 1 decisions applied
#'
#' Skipped when `decisions.fixed` is `{}` and `decisions.random_effects_removed`
#' is empty ([step2_runs()]): the narrow IGS of step 1 is used instead.
#' Otherwise the base model with those parameters fixed and those random
#' effects removed, sampled over step 1's narrow ranges (the ranges of removed
#' random effects are unused). The model keeps the reduced variability from
#' here on: steps 3 and 4 and the final model use it too. Writes diagnostics,
#' the best run's estimates and an identifiability summary of the parameters
#' still estimated
#' ([identifiability_summary()]). Checkpoint 1b (`decisions.fixed_confirmed`)
#' follows: confirm that nothing else needs fixing or removing before step 3.
#'
#' @inheritParams run_step1
#' @return The IGS result used by step 3 (step 2's, or step 1's narrow IGS),
#'   with attribute `skipped`.
#' @export
run_step2 <- function(config, rerun = FALSE) {
  validate_config(config, 2)
  step1 <- run_step1(config)
  if (!step2_runs(config)) {
    message("Step 2: skipped (nothing fixed or removed at checkpoint 1); using step 1's narrow IGS")
    return(structure(step1$narrow, skipped = TRUE))
  }
  dir <- pipeline_paths(config)$step2
  if (use_cached(config, 2, dir, rerun)) return(readRDS(file.path(dir, "results.rds")))

  spec <- model_spec(config)
  setup_project(spec)
  save_project_as(file.path(dir, "fixed.mlxtran"))
  pp <- lixoftConnectors::getPopulationParameterInformation()$name
  ranges <- ranges_from_narrow(step1$narrow_ranges, attr(step1$sampling_ranges, "log10"))
  rm <- removed_random_effects(config)
  changes <- c(if (length(spec$fixed)) paste(paste(names(spec$fixed), collapse = ", "), "fixed"),
               if (length(rm)) paste("no random effect on", paste(rm, collapse = ", ")))
  message("Step 2: initial guess sampling with ", paste(changes, collapse = " and "), " (", config$step2_igs$nits, " runs)")
  igs <- run_igs(ranges, config$step2_igs$nits, config$analysis$r_seed + 2, set_error = error_starts(spec, pp),
                 set_variance = config$step1_igs$narrow$omega_initial, fixed = spec$fixed)
  save_project_as(file.path(dir, "fixed.mlxtran"))
  diagnostics(config, igs, file.path(dir, "diagnostics"))
  utils::write.csv(best_run_estimates(igs$tabestimates), file.path(dir, "best_estimates.csv"), row.names = FALSE)
  utils::write.csv(identifiability_summary(igs, p_threshold = config$step1_igs$narrow$p_threshold),
                   file.path(dir, "identifiability_summary.csv"), row.names = FALSE)
  out <- save_step(structure(igs, skipped = FALSE), config, 2, dir)
  figures_page(config, 2)
  out
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
#' once each candidate from [resolve_candidates()]. `scm_best` and
#' `cossac_best` aren't fitted again: they're copied from the search's final
#' model, which was already fitted with the analysis scenario.
#'
#' @inheritParams run_step1
#' @param rerun Refit the covariate search (and so the candidates too).
#' @param rerun_candidates Refit only the candidates, keeping the search, e.g.
#'   after changing `candidate_models`.
#' @return List with `searches` (model-building results per method),
#'   `candidates` (covariate models), `initial_values` and `single_fits`
#'   (each candidate's single fit: BICc, condition number of the correlation
#'   matrix of the estimates, a note naming treatment effects with Wald
#'   p >= 0.05, and whether it was fitted or copied).
#' @export
run_step3 <- function(config, rerun = FALSE, rerun_candidates = rerun) {
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
        linearization = isTRUE(s3$use_linearization), relationships = model_building_relationships(s3$relationships),
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
    figures_page(config, 3)   # the search, for checkpoint 2
  }

  # checkpoint 2; while it's undecided, flag scm_delta2 covariates that Monolix's backward step removed
  if (is.null(s3$candidate_models) && "scm_delta2" %in% auto_candidates(config) && !is.null(searches$SCM)) {
    note <- scm_delta2_backward_note(searches$SCM, s3$scm_stop_rule_delta_bicc %||% 2)
    if (!is.null(note)) warning(note, call. = FALSE)
  }
  validate_config(config, 3)
  cdir <- paths$step3_candidates
  if (use_cached(config, 3, cdir, isTRUE(rerun) || isTRUE(rerun_candidates))) return(readRDS(file.path(cdir, "results.rds")))
  candidates <- resolve_candidates(config, searches)
  writeLines(as_yaml_text(list(candidates = unclass(candidates)[names(candidates)],
                               duplicates = attr(candidates, "duplicates"))),
             file.path(cdir, "candidates.yaml"))
  single <- list()
  for (nm in names(candidates)) {
    cov_model <- format_covariate_model(candidates[[nm]])
    file <- file.path(cdir, paste0(nm, ".mlxtran"))
    # a search's best model was already fitted with the analysis scenario (covariate_model_building()):
    # copy it instead of fitting it again
    src <- search_final_model(paths$step3_search, nm)
    if (!is.null(src)) {
      message("Step 3: candidate ", nm, " (", cov_model, ") copied from ", src)
      single[[nm]] <- tryCatch({
        lixoftConnectors::loadProject(normalizePath(src))
        row <- single_fit_row(nm, cov_model, "copied from the search")
        save_project_as(file)
        row
      }, error = function(e) {
        message("Could not copy ", src, " (", conditionMessage(e), "); fitting ", nm, " instead")
        NULL
      })
      if (!is.null(single[[nm]])) next
    }
    message("Step 3: building candidate ", nm, " (", cov_model, ")")
    setup_from_best(config, candidates[[nm]], iv, file)
    single[[nm]] <- tryCatch({
      lixoftConnectors::runScenario()
      save_project_as(file)
      single_fit_row(nm, cov_model, "fitted from step 2's best values")
    }, error = function(e) {
      message("Candidate ", nm, " failed to fit from step 2's best values: ", conditionMessage(e))
      data.frame(model = nm, covariate_model = cov_model, bicc = NA_real_, condition_number = NA_real_,
                 wald_note = "", source = "fit failed")
    })
  }
  single_fits <- single_fit_table(do.call(rbind, single))
  utils::write.csv(single_fits, file.path(cdir, "single_fits.csv"), row.names = FALSE)
  out <- save_step(list(searches = searches, candidates = candidates, initial_values = iv, single_fits = single_fits),
                   config, 3, cdir)
  figures_page(config, 3)   # adds the candidate fits, for checkpoint 2b
  out
}

# One row of single_fits.csv, from the loaded project's results: BICc, the condition number of
# the correlation matrix of the estimates (FIM), and a note on treatment effects with Wald p >= 0.05
single_fit_row <- function(model, covariate_model, source) {
  methods <- result_methods()
  bicc <- read_likelihood(methods)[["bicc"]]
  cn <- tryCatch(condition_number_of(lixoftConnectors::getCorrelationOfEstimates()[[methods$se]]),
                 error = function(e) NA_real_)
  note <- tryCatch(wald_note(lixoftConnectors::getTests()$wald, se_method = methods$se),
                   error = function(e) paste("Wald tests not read:", conditionMessage(e)))
  data.frame(model = model, covariate_model = covariate_model, bicc = bicc, condition_number = cn,
             wald_note = note, source = source)
}

# Step 3's single fits, best BICc first, with the difference from the best
single_fit_table <- function(sf) {
  sf$delta_bicc <- sf$bicc - min(sf$bicc, na.rm = TRUE)
  for (k in c("condition_number", "wald_note", "source")) if (is.null(sf[[k]])) sf[[k]] <- NA
  sf$wald_note[is.na(sf$wald_note)] <- ""
  sf[order(sf$bicc, na.last = TRUE),
     c("model", "covariate_model", "bicc", "delta_bicc", "condition_number", "wald_note", "source")]
}

# The search's final model for an automatic best-model candidate (scm_best, cossac_best), or NULL.
# covariate_model_building() saves it as <method>/final_model.mlxtran next to <method>.mlxtran.
search_final_model <- function(search_dir, candidate) {
  method <- c(scm_best = "scm", cossac_best = "cossac")[candidate]
  if (is.na(method)) return(NULL)
  f <- file.path(search_dir, method, "final_model.mlxtran")
  if (file.exists(f)) f else NULL
}

#' Step 4: paired comparison of the candidate models
#'
#' Checkpoint 2b first: while `step4_model_comparison.compare` is `~`, prints
#' step 3's single fits (one fit per candidate) and stops, so poor candidates
#' can be left out (`compare: all`, or a list; see [step4_models()]). Then draws
#' one set of `nits` initial values (`_pop` over step 1's narrow ranges,
#' omegas at `omega_initial`, treatment effects at `treatment_initial`, fixed
#' parameters fixed), fits every compared candidate from the same values, and
#' compares them with [compare_models()]. Checkpoint 3 (`final_model`)
#' follows.
#'
#' @inheritParams run_step1
#' @return List with `igs` (per candidate), `comparison` and `initial_values`.
#' @export
run_step4 <- function(config, rerun = FALSE) {
  validate_config(config, 3)
  step1 <- run_step1(config)
  step3 <- run_step3(config)
  cdir <- pipeline_paths(config)$step3_candidates
  # checkpoint 2b: which candidates go on to the comparison, decided from step 3's single fits
  if (is.null(config$step4_model_comparison$compare)) {
    sf <- step3$single_fits %||% single_fit_table(utils::read.csv(file.path(cdir, "single_fits.csv")))
    message("Step 3 single fits (one fit each; ", file.path(cdir, "single_fits.csv"), "):\n",
            paste(utils::capture.output(print(sf, row.names = FALSE)), collapse = "\n"))
  }
  validate_config(config, 4, decisions_needed = FALSE)
  dir <- pipeline_paths(config)$step4
  if (use_cached(config, 4, dir, rerun)) return(readRDS(file.path(dir, "results.rds")))
  s4 <- config$step4_model_comparison
  models <- step4_models(config, step3$candidates)

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
  save_comparison_plots(comparison, dir)
  out <- save_step(list(igs = igs, comparison = comparison, initial_values = iv), config, 4, dir)
  figures_page(config, 4)
  out
}

#' The final model and its report
#'
#' After checkpoint 3: copies the chosen candidate's step 4 project (left at
#' its best run) to `final_model/final_model.mlxtran` and writes its report.
#' Also draws the paired-estimates, normalized boxplot and RSE plots of the
#' final model's step 4 runs ([save_final_model_plots()]), every time it's
#' called.
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
  if (rerun && dir.exists(dir)) {
    old <- tryCatch(yaml::read_yaml(file.path(dir, "decision.yaml")), error = function(e) NULL)
    if (!is.null(old)) note_replaced(dir, if (!identical(old$final_model, final)) paste0("final_model: ", old$final_model, " -> ", final))
    unlink(dir, recursive = TRUE)
  }
  # plots of the final model's step 4 runs: quick and without Monolix, so redrawn every time
  pl <- config$step1_igs$plots
  save_final_model_plots(step4$igs[[final]]$tabestimates, dir, remove_bicc_outliers = isTRUE(pl$remove_bicc_outliers),
                         remove_estimate_outliers = isTRUE(pl$remove_estimate_outliers))
  if (!rerun && file.exists(file)) {
    message("Final model: ", file, " already exists (rerun to remake it)")
    return(invisible(file))
  }
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
  record_history(config, dir, "final")
  invisible(file)
}
