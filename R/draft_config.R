# Write a pre-filled analysis.yaml from a model file and a data file.

#' Draft an analysis.yaml for a new analysis
#'
#' Reads the data file and parses the structural model (no Monolix needed),
#' then writes a config with every key filled in: guessed values are marked
#' `# guessed: <reason>, check`, values only you can supply are `~` and marked
#' `# REQUIRED`, and everything else uses [config_defaults()]. The file is
#' meant to be edited; [validate_config()] refuses to run until every
#' REQUIRED value is filled.
#'
#' @param model_file Path to the Mlxtran structural model.
#' @param data_file Path to the data file.
#' @param name Analysis name (used for the default output folder).
#' @param path Where to write the YAML. Paths inside it are written relative to
#'   this file's folder where possible.
#' @param overwrite Overwrite an existing file?
#' @return `path`, invisibly.
#' @export
draft_config <- function(model_file, data_file, name = NULL, path = "analysis.yaml", overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop(path, " already exists. Use overwrite = TRUE to replace it.", call. = FALSE)
  }
  parsed <- parse_mlxtran(model_file)
  data <- read_monolix_data(data_file)
  ht <- guess_header_types(data)
  transforms <- guess_report_transforms(parsed)
  d <- config_defaults()

  obs <- ht$column[ht$type == "observation"]
  covs <- ht$column[ht$type %in% c("catcov", "contcov")]
  catcovs <- ht$column[ht$type == "catcov"]
  # the report colours groups with Dark2 (8 colours): only covariates with at most 8 levels
  n_levels <- vapply(catcovs, function(cc) length(unique(stats::na.omit(data[[cc]]))), numeric(1))
  strat <- catcovs[!grepl("^int", catcovs) & n_levels <= 8]

  # examples in the comments use this model's own parameters and covariates
  p1 <- parsed$parameters[1]
  p2 <- if (length(parsed$parameters) > 1) parsed$parameters[2] else p1
  c1 <- if (length(covs)) covs[1] else "cov1"
  c2 <- if (length(covs) > 1) covs[2] else c1

  yaml_dir <- dirname(normalizePath(path, mustWork = FALSE))
  L <- character()
  ln <- function(...) L <<- c(L, paste0(...))
  cm <- function(text) if (nzchar(text)) paste0("  # ", text) else ""

  ln("# analysis.yaml drafted by monolixPipelineR::draft_config() on ", format(Sys.Date()))
  ln("# Fill in every REQUIRED value and check every guessed one, then run validate_config().")
  ln("# Keys you leave out use the defaults from config_defaults().")
  ln("# `~` means not set or not decided yet; `{}` means none.")
  ln("")
  ln("analysis:")
  ln("  name: ", yv(name), cm(if (is.null(name)) "REQUIRED" else ""))
  ln("  out_folder: ", yv(if (is.null(name)) NULL else file.path("results", name)),
     cm("relative to this file; default results/<name>"))
  ln("  r_seed: ", yv(d$analysis$r_seed), cm("seed for sampling initial values (the Monolix seed is fixed)"))
  ln("")
  ln("monolix:")
  ln("  version: ~", cm("optional, e.g. \"2024R1\"; the install path comes from the MONOLIX_PATH environment variable"))
  ln("")
  ln("data:")
  ln("  file: ", yv(rel_to(data_file, yaml_dir)))
  ln("  header_types:")
  for (i in seq_len(nrow(ht))) {
    ln("    ", yk(ht$column[i]), ": ", yv(ht$type[i]), cm(paste0("guessed: ", ht$reason[i], ", check")))
  }
  ln("")
  ln("model:")
  ln("  structural_model: ", yv(rel_to(model_file, yaml_dir)))
  if (length(obs)) {
    ln("  error_model:", cm("constant, proportional, combined1 or combined2"))
    for (o in obs) ln("    ", yk(o), ": constant", cm("guessed: default, check"))
    ln("  error_initial: {}", cm("starting values for error-model parameters, which are still estimated; default Monolix's: a = 1, b = 0.3, c = 1"))
  } else {
    ln("  error_model: ~", cm("REQUIRED: no observation column was found"))
  }
  ln("  random_effects_excluded: []", cm(paste0("parameters with no random effect; model parameters: ",
                                                 paste(parsed$parameters, collapse = ", "))))
  ln("  distributions: {}", cm(paste0("e.g. {", p1, ": normal}; default lognormal")))
  if (length(transforms)) {
    ln("  report_transforms:", cm("how estimates are shown in reports"))
    for (k in names(transforms)) {
      src <- parsed$rescalings[[sub("_pop$", "", k)]]
      ln("    ", yk(k), ": ", flow(transforms[[k]]),
         cm(paste0("guessed from `", src$derived, " = ", sub("^x", sub("_pop$", "", k), src$expr), "`, check")))
    }
  } else {
    ln("  report_transforms: {}", cm(paste0("e.g. {", p1, "_pop: {label: \"", p1, " (units)\", expr: \"x * 10^5\"}}")))
  }
  ln("")
  ln("estimation:", cm("Monolix's SAEM settings"))
  ln("  saem: ", flow(d$estimation$saem))
  ln("")
  sc <- d$scenario
  ln("scenario:", cm("Monolix tasks run for every fit"))
  ln("  population_parameters: ", yv(sc$population_parameters), cm("required"))
  ln("  conditional_mode: ", yv(sc$conditional_mode), cm("required: individual estimates used by the pipeline"))
  ln("  conditional_distribution: ", yv(sc$conditional_distribution), cm("conditional mean and its shrinkage"))
  ln("  standard_errors: ", yv(sc$standard_errors), cm("false: RSEs and condition numbers are NA"))
  ln("  log_likelihood: ", yv(sc$log_likelihood), cm("required: BICc ranks runs and models"))
  ln("  linearization: ", yv(sc$linearization), cm("standard errors and likelihood; false = stochastic approximation + importance sampling"))
  ln("  plots: ", yv(sc$plots))
  ln("")
  ln("step1_igs:", cm("convergence assessment without covariates"))
  ln("  initial_ranges:", cm("REQUIRED: plausible min/max for each parameter; range_around(x) gives [x/100, 100x]"))
  for (p in parsed$parameters) ln("    ", yk(p), ": {min: ~, max: ~, log10: true}", cm("REQUIRED"))
  ln("  omega_range: ", flow(d$step1_igs$omega_range), cm("wide search, random-effect SDs"))
  ln("  wide_nits: ", yv(d$step1_igs$wide_nits))
  ln("  narrow: ", flow(d$step1_igs$narrow), cm("min_fold_width: e.g. 2 widens ranges to at least [best/2, 2*best]"))
  ln("  profile_likelihood: {}", cm(paste0("run after step 1, e.g. {", p1, "_pop: {range: [0.5, 20], n: 10, log10: true}} or {", p1, "_pop: {values: [1, 2, 3, 5]}}")))
  ln("  plots: ", flow(d$step1_igs$plots))
  ln("")
  ln("decisions:", cm("filled in by you at the checkpoints"))
  ln("  fixed: ~", cm(paste0("CHECKPOINT 1: e.g. {", p1, "_pop: 50", if (p2 != p1) paste0(", ", p2, "_pop: 0.7"), "}, or {} to fix none")))
  ln("  fixed_rationale: ~")
  ln("  random_effects_removed: []", cm(paste0("CHECKPOINT 1 (optional): parameters to drop the random effect of from step 2 on, e.g. [", parsed$parameters[1], "]")))
  ln("  fixed_confirmed: ~", cm("CHECKPOINT 1b (after step 2): true if nothing else needs fixing or removing; set back to ~ when you change fixed or random_effects_removed"))
  ln("")
  ln("step2_igs: ", flow(list(nits = d$step2_igs$nits)), cm("skipped when decisions.fixed is {} and random_effects_removed is []"))
  ln("")
  s3 <- d$step3_model_building
  ln("step3_model_building:")
  ln("  method: ", yv(s3$method), cm("SCM, COSSAC, [SCM, COSSAC] to run and compare both, or none to compare only candidate_models"))
  ln("  criterion: ", yv(s3$criterion), cm("Monolix's name for BICc; or LRT"))
  ln("  lrt_thresholds: ", flow(s3$lrt_thresholds), cm("only used with LRT"))
  if (length(covs)) {
    ln("  covariates: ", flow(covs), cm("guessed: all covariate columns, check"))
  } else {
    ln("  covariates: all")
  }
  ln("  parameters: all", cm("or a list of parameters to search"))
  ln("  relationships: ~", cm(paste0("or {include: {", p1, ": [", c1, "]}, exclude: {", p2, ": [", c2, "]}}: pairs in every model / in none")))
  ln("  use_linearization: ", yv(s3$use_linearization))
  ln("  initial_value_rounding: ", yv(s3$initial_value_rounding), cm("significant figures for _pop initial values"))
  ln("  omega_initial: ", yv(s3$omega_initial))
  ln("  scm_stop_rule_delta_bicc: ", yv(s3$scm_stop_rule_delta_bicc), cm("for the automatic scm_delta2 candidate (criterion BIC only)"))
  ln("  candidate_models: ~", cm(paste0("CHECKPOINT 2: e.g. {model1: {", p1, ": [", c1, "], ", p2, ": [", c1, "]}}, or {}; ",
                                        "scm_best and scm_delta2 (SCM; scm_delta2 with BIC only) and cossac_best (COSSAC) are added automatically")))
  ln("")
  s4 <- d$step4_model_comparison
  ln("step4_model_comparison:")
  ln("  nits: ", yv(s4$nits))
  ln("  compare: ~", cm("CHECKPOINT 2b (after step 3's single fits): all, or a list of candidates to compare"))
  ln("  omega_initial: ", yv(s4$omega_initial))
  ln("  treatment_initial: ", yv(s4$treatment_initial))
  ln("  fdr_method: ", yv(s4$fdr_method))
  ln("  final_model: ~", cm("CHECKPOINT 3"))
  ln("  final_rationale: ~")
  ln("")
  ln("report:")
  if (length(strat)) {
    ln("  stratify_by: ", flow(strat), cm("guessed: categorical covariates with at most 8 levels, not starting with int, check"))
  } else {
    ln("  stratify_by: ~")
  }
  ln("  group_labels: ~", cm("names for the groups, e.g. {\"1.0\": treated, \"0.0\": control}"))

  writeLines(L, path)

  required <- c(if (is.null(name)) "analysis.name", "step1_igs.initial_ranges min/max",
                if (!length(obs)) "model.error_model")
  message("Wrote ", path, "\n  REQUIRED: ", paste(required, collapse = "; "),
          "\n  Check the guessed header types, error model",
          if (length(transforms)) ", report transforms", ", covariates and stratification.")
  invisible(path)
}

# --- YAML text helpers (the yaml package can't write comments) ---

# Scalar value as YAML text that reads back as the same R value
yv <- function(x) {
  if (is.null(x) || length(x) == 0) return("~")
  if (is.logical(x)) return(if (isTRUE(x)) "true" else "false")
  if (is.numeric(x)) return(fmt_num(x))
  sub("\n$", "", yaml::as.yaml(as.character(x)))
}

# Map key, quoted when YAML would otherwise read it as something else
yk <- function(x) sub("\n$", "", yaml::as.yaml(as.character(x)))

# Numbers that YAML reads back as numbers (it reads 1e-7 as text but 1.0e-07 as a number)
fmt_num <- function(x) {
  if (x == round(x) && abs(x) < 1e15) return(format(x, scientific = FALSE))
  s <- format(signif(x, 6), scientific = abs(x) < 1e-3 || abs(x) >= 1e6)
  if (grepl("e", s) && !grepl("\\.", s)) s <- sub("e", ".0e", s)
  s
}

# Flow-style YAML for a vector or (nested) list
flow <- function(x) {
  if (is.list(x) && !is.null(names(x))) {
    if (!length(x)) return("{}")
    return(paste0("{", paste0(vapply(names(x), yk, ""), ": ", vapply(x, flow, ""), collapse = ", "), "}"))
  }
  if (is.list(x) || length(x) > 1) {
    if (!length(x)) return("[]")
    return(paste0("[", paste(vapply(x, yv, ""), collapse = ", "), "]"))
  }
  yv(x)
}

# Path relative to `dir` when it is inside it, otherwise absolute
rel_to <- function(target, dir) {
  t <- normalizePath(target, mustWork = FALSE)
  d <- normalizePath(dir, mustWork = FALSE)
  if (startsWith(t, paste0(d, "/"))) substring(t, nchar(d) + 2) else t
}
