# Run logs: a machine-written record of what each step actually ran, and a check
# that cached results still match the current config.

#' Settings that affect the results of a step
#'
#' Everything a step's fits depend on, including earlier steps. Labels,
#' rationales, plot options, report settings and file locations are left out,
#' so changing them doesn't invalidate cached fits (file contents are checked
#' by hash instead). `step1_igs.profile_likelihood` is left out too: the
#' profiles run after step 1 ([run_profiles()]) and are checked on their own.
#'
#' @param config Output of [read_analysis_config()].
#' @param step 1 to 4.
#' @param before_checkpoint Leave out the decision made partway through the
#'   step (step 3's `candidate_models`), for the part of the step that runs
#'   before it (the covariate search).
#' @return Nested list.
#' @export
step_settings <- function(config, step, before_checkpoint = FALSE) {
  stopifnot(step %in% 1:4)
  model <- config$model[c("error_model", "error_initial", "random_effects_excluded", "distributions")]
  out <- list(
    analysis = list(r_seed = config$analysis$r_seed),
    data = list(header_types = config$data$header_types),
    model = model,
    estimation = config$estimation,
    scenario = config$scenario,
    step1_igs = config$step1_igs[setdiff(names(config$step1_igs), c("plots", "profile_likelihood"))]
  )
  if (step >= 2) {
    out$decisions <- list(fixed = config$decisions$fixed)
    out$step2_igs <- config$step2_igs
  }
  if (step >= 3) {
    s3 <- config$step3_model_building
    out$step3_model_building <- if (step == 3 && before_checkpoint) s3[setdiff(names(s3), "candidate_models")] else s3
  }
  if (step >= 4) {
    s4 <- config$step4_model_comparison
    out$step4_model_comparison <- s4[setdiff(names(s4), c("final_model", "final_rationale"))]
  }
  out
}

#' Write the run log for a step
#'
#' Records the settings the step used ([step_settings()]), the seeds, MD5
#' hashes of the data and model files, software versions and `sessionInfo()`.
#' Written by the pipeline only; don't edit it.
#'
#' @param config Output of [read_analysis_config()].
#' @param step 1 to 4.
#' @param dir Folder of the step's results.
#' @param extra Optional named list of anything else to record.
#' @inheritParams step_settings
#' @return Path to the log, invisibly.
#' @export
write_run_log <- function(config, step, dir, extra = list(), before_checkpoint = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  log <- list(
    step = step,
    written = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    config_file = attr(config, "path"),
    seeds = list(r = config$analysis$r_seed, monolix = MONOLIX_SEED),
    files = file_hashes(config),
    versions = list(
      R = R.version.string,
      monolixPipelineR = as.character(utils::packageVersion("monolixPipelineR")),
      monolixPipelineR_commit = package_commit(),
      lixoftConnectors = if (monolix_available()) as.character(utils::packageVersion("lixoftConnectors")) else NULL,
      monolix = monolix_version() %||% config$monolix$version
    ),
    settings = step_settings(config, step, before_checkpoint),
    extra = extra,
    session_info = utils::capture.output(utils::sessionInfo())
  )
  path <- file.path(dir, "run_log.yaml")
  writeLines(as_yaml_text(log), path)
  invisible(path)
}

#' Check that a step's cached results match the current config
#'
#' Compares the settings and file hashes in the step's `run_log.yaml` with the
#' current config. Stops and lists every difference if they don't match, so
#' results from an older config are never reused by accident.
#'
#' @inheritParams write_run_log
#' @return `TRUE` if a matching log exists, `FALSE` if there is no log yet.
#' @export
check_run_log <- function(config, step, dir, before_checkpoint = FALSE) {
  path <- file.path(dir, "run_log.yaml")
  if (!file.exists(path)) return(FALSE)
  log <- yaml::read_yaml(path, handlers = yaml_handlers())
  then <- c(flatten_settings(log$settings), flatten_settings(list(files = log$files)))
  # logs written before the profiles became their own step recorded them with step 1
  then <- then[!startsWith(names(then), "step1_igs.profile_likelihood")]
  now <- c(flatten_settings(yaml_roundtrip(step_settings(config, step, before_checkpoint))),
           flatten_settings(list(files = file_hashes(config))))
  keys <- union(names(then), names(now))
  differ <- keys[vapply(keys, function(k) !identical(then[[k]], now[[k]]), logical(1))]
  if (length(differ)) {
    lines <- vapply(differ, function(k) paste0("  - ", k, ": ", then[[k]] %||% "(not set)", " -> ", now[[k]] %||% "(not set)"), "")
    stop("The cached results in ", dir, " were made with different settings:\n", paste(lines, collapse = "\n"),
         "\nSet RERUN for this step to TRUE (which overwrites them), or restore the old settings.", call. = FALSE)
  }
  TRUE
}

file_hashes <- function(config) {
  h <- function(p) if (!is.null(p) && file.exists(p)) unname(tools::md5sum(p)) else NULL
  list(data = h(config$data$file), model = h(config$model$structural_model))
}

# Write and read back through YAML so values compare like-for-like with a stored log
yaml_roundtrip <- function(x) yaml::yaml.load(as_yaml_text(x), handlers = yaml_handlers())

# YAML text with true/false for logicals (the yaml package writes yes/no by default)
as_yaml_text <- function(x) {
  yaml::as.yaml(x, handlers = list(logical = function(v) {
    out <- ifelse(is.na(v), "~", ifelse(v, "true", "false"))
    class(out) <- "verbatim"
    out
  }))
}

# Nested list -> named character vector, e.g. c("step1_igs.narrow.nits" = "100")
flatten_settings <- function(x, prefix = "") {
  if (is.null(x)) return(stats::setNames(list("~"), sub("\\.$", "", prefix)))
  if (is.list(x) && !is.null(names(x)) && length(x)) {
    return(do.call(c, lapply(names(x), function(k) flatten_settings(x[[k]], paste0(prefix, k, ".")))))
  }
  if (is.list(x) && !length(x)) return(stats::setNames(list("{}"), sub("\\.$", "", prefix)))
  val <- paste(vapply(unlist(x), function(v) format(v, digits = 15), ""), collapse = ", ")
  stats::setNames(list(val), sub("\\.$", "", prefix))
}

package_commit <- function() {
  sha <- utils::packageDescription("monolixPipelineR")$RemoteSha
  if (!is.null(sha)) return(sha)
  dir <- tryCatch(dirname(system.file("DESCRIPTION", package = "monolixPipelineR")), error = function(e) "")
  if (!nzchar(dir) || !dir.exists(file.path(dir, ".git"))) return(NULL)
  out <- tryCatch(suppressWarnings(system2("git", c("-C", shQuote(dir), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
                  error = function(e) character())
  if (length(out) == 1) out else NULL
}

monolix_version <- function() {
  if (!monolix_available()) return(NULL)
  # quietly: lixoftConnectors prints a warning when it hasn't been initialized
  state <- NULL
  utils::capture.output(state <- tryCatch(suppressWarnings(suppressMessages(lixoftConnectors::getLixoftConnectorsState(quietly = TRUE))),
                                          error = function(e) NULL))
  if (is.list(state) && !is.null(state$version)) as.character(state$version) else NULL
}
