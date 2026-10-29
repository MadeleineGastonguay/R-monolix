# A history of every step the pipeline has finished, in <out_folder>/history.yaml. Results are
# overwritten when a step is rerun; the history keeps a short record of each version: when it
# ran, what changed from the version it replaced, the decisions it used and its main result.
# Written by the pipeline only, like the run logs.

# Changes between the run log in `dir` and the current config, as "key: old -> new" lines
# (character(0) if nothing differs, NULL if there's no run log)
run_log_changes <- function(config, step, dir, before_checkpoint = FALSE) {
  path <- file.path(dir, "run_log.yaml")
  if (!file.exists(path)) return(NULL)
  log <- yaml::read_yaml(path, handlers = yaml_handlers())
  then <- c(flatten_settings(log$settings), flatten_settings(list(files = log$files)))
  # logs written before the profiles became their own step recorded them with step 1
  then <- then[!startsWith(names(then), "step1_igs.profile_likelihood")]
  now <- c(flatten_settings(yaml_roundtrip(step_settings(config, step, before_checkpoint))),
           flatten_settings(list(files = file_hashes(config))))
  keys <- union(names(then), names(now))
  differ <- keys[vapply(keys, function(k) !identical(then[[k]], now[[k]]), logical(1))]
  structure(vapply(differ, function(k) paste0(k, ": ", then[[k]] %||% "(not set)", " -> ", now[[k]] %||% "(not set)"), ""),
            names = NULL, written = log$written)
}

# What a rerun replaced, noted just before the old folder is deleted and used when the new
# version is recorded
history_state <- new.env(parent = emptyenv())

note_replaced <- function(dir, changes, previous_run = NULL) {
  history_state[[normalizePath(dir, mustWork = FALSE)]] <- list(
    previous_run = previous_run %||% attr(changes, "written"),
    changed = if (length(changes)) as.list(changes) else "nothing (rerun with the same settings)")
}

take_replaced <- function(dir) {
  key <- normalizePath(dir, mustWork = FALSE)
  out <- history_state[[key]]
  if (!is.null(out)) rm(list = key, envir = history_state)
  out
}

history_path <- function(config) file.path(pipeline_paths(config)$root, "history.yaml")

#' Read the pipeline's history
#'
#' Every time a step finishes, the pipeline appends an entry to
#' `<out_folder>/history.yaml`: the part of the pipeline (e.g. `step2_igs`),
#' when it finished, what changed from the version it replaced (`replaced`,
#' only for reruns), the decisions it used and its main result (best BICc,
#' best model, ...). Rerunning a step overwrites its results, so this is the
#' record of the versions you tried. Failed runs aren't recorded.
#'
#' @param config Output of [read_analysis_config()].
#' @return List of entries, oldest first (empty if there's no history yet).
#' @export
read_history <- function(config) {
  path <- history_path(config)
  if (!file.exists(path)) return(list())
  yaml::read_yaml(path, handlers = yaml_handlers())$runs %||% list()
}

# Append one entry; never stops the pipeline
record_history <- function(config, dir, step, results = NULL, decisions = NULL, result = NULL) {
  tryCatch({
    root <- normalizePath(pipeline_paths(config)$root, mustWork = FALSE)
    full <- normalizePath(dir, mustWork = FALSE)
    part <- if (startsWith(full, paste0(root, "/"))) substring(full, nchar(root) + 2) else basename(full)
    entry <- list(part = part, step = if (is.numeric(step)) as.integer(step) else step, finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  replaced = take_replaced(dir),
                  decisions = decisions %||% history_decisions(config, part),
                  result = result %||% history_result(part, results))
    entry <- entry[!vapply(entry, is.null, logical(1))]
    runs <- c(read_history(config), list(entry))
    dir.create(dirname(history_path(config)), recursive = TRUE, showWarnings = FALSE)
    writeLines(c("# Written by the pipeline: one entry per finished step, oldest first. Don't edit.",
                 as_yaml_text(list(runs = runs))), history_path(config))
  }, error = function(e) message("Could not update history.yaml: ", conditionMessage(e)))
  invisible(NULL)
}

# The decisions each part of the pipeline depends on
history_decisions <- function(config, part) {
  d <- config$decisions
  s3 <- config$step3_model_building
  s4 <- config$step4_model_comparison
  keep <- function(x) { x <- x[!vapply(x, is.null, logical(1))]; if (length(x)) x else NULL }
  switch(sub("^.*/", "", part),
    step2_igs = keep(list(fixed = d$fixed, random_effects_removed = d$random_effects_removed,
                          fixed_rationale = d$fixed_rationale)),
    search = keep(list(method = s3$method, criterion = s3$criterion, covariates = s3$covariates,
                       parameters = s3$parameters, relationships = s3$relationships)),
    candidates = keep(list(candidate_models = s3$candidate_models)),
    step4_model_comparison = keep(list(compare = s4$compare)),
    final_model = keep(list(final_model = s4$final_model, final_rationale = s4$final_rationale)),
    NULL)
}

# The main result of each part, from what it saved
history_result <- function(part, results) {
  if (is.null(results)) return(NULL)
  best_igs <- function(igs) {
    te <- igs$tabestimates
    list(runs = length(unique(te$run)), best_bicc = signif(min(te$bicc, na.rm = TRUE), 6),
         best_run = as.integer(te$run[which.min(te$bicc)]))
  }
  tryCatch(switch(sub("^.*/", "", part),
    step1_igs = list(wide = best_igs(results$wide), narrow = best_igs(results$narrow)),
    step2_igs = best_igs(results),
    search = lapply(results, function(mb) format_covariate_model(scm_best_model(mb))),
    candidates = {
      sf <- results$single_fits
      stats::setNames(as.list(signif(sf$bicc, 6)), sf$model)
    },
    step4_model_comparison = {
      s <- results$comparison$summary
      list(best_model = as.character(s$model[which.max(s$times_best)]),
           times_best = stats::setNames(as.list(s$times_best), as.character(s$model)),
           min_bicc = stats::setNames(as.list(signif(s$min_bicc, 6)), as.character(s$model)))
    },
    NULL), error = function(e) NULL)
}
