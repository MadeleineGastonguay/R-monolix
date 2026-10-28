# The Monolix scenario: which tasks run, and whether standard errors and the
# likelihood use linearization. Set from the `scenario` section of analysis.yaml.

# analysis.yaml key -> Monolix task keyword
SCENARIO_TASKS <- c(
  population_parameters = "populationParameterEstimation",
  conditional_mode = "conditionalModeEstimation",
  conditional_distribution = "conditionalDistributionSampling",
  standard_errors = "standardErrorEstimation",
  log_likelihood = "logLikelihoodEstimation",
  plots = "plots"
)

# Tasks the pipeline can't work without: estimates, individual (conditional mode)
# predictions, and BICc to rank runs and models
REQUIRED_TASKS <- c("population_parameters", "conditional_mode", "log_likelihood")

#' Set the Monolix scenario of the loaded project
#'
#' @param scenario The `scenario` section of the config (see
#'   [config_defaults()]): one true/false per task plus `linearization`.
#' @export
apply_scenario <- function(scenario) {
  require_monolix()
  sc <- lixoftConnectors::getScenario()
  for (k in names(SCENARIO_TASKS)) {
    if (!is.null(scenario[[k]])) sc$tasks[SCENARIO_TASKS[[k]]] <- isTRUE(scenario[[k]])
  }
  sc$linearization <- isTRUE(scenario$linearization)
  lixoftConnectors::setScenario(sc)
  invisible(sc)
}

#' Which method's results to read
#'
#' With linearization, standard errors and the likelihood are stored under
#' `linearization`; without it, under `stochasticApproximation` and
#' `importanceSampling`.
#'
#' @param linearization `TRUE`/`FALSE`, or `NULL` to read it from the loaded
#'   project's scenario.
#' @return List with `se` and `ll`, the method names to read results from.
#' @export
result_methods <- function(linearization = NULL) {
  if (is.null(linearization)) {
    require_monolix()
    linearization <- isTRUE(lixoftConnectors::getScenario()$linearization)
  }
  if (linearization) list(se = "linearization", ll = "linearization")
  else list(se = "stochasticApproximation", ll = "importanceSampling")
}

# Relative standard errors by parameter, or NA when standard errors weren't computed
read_rses <- function(methods, parameters) {
  se <- tryCatch(lixoftConnectors::getEstimatedStandardErrors()[[methods$se]], error = function(e) NULL)
  if (is.null(se)) return(stats::setNames(rep(NA_real_, length(parameters)), parameters))
  stats::setNames(se$rse, se$parameter)[parameters]
}

# OFV, AIC, BIC and BICc from the scenario's likelihood method
read_likelihood <- function(methods) {
  ll <- lixoftConnectors::getEstimatedLogLikelihood()[[methods$ll]]
  if (is.null(ll)) stop("No ", methods$ll, " log-likelihood; is the log_likelihood task on?", call. = FALSE)
  ll <- unlist(ll)
  c(ofv = ll[[1]], aic = ll[[2]], bic = ll[[3]], bicc = ll[[4]])
}

# Eta shrinkage for whichever individual estimates were computed (others are NA)
read_shrinkage <- function(run) {
  sh <- lixoftConnectors::getEtaShrinkage()
  kinds <- c("conditionalMean", "conditionalMode", "conditionalDist")
  first <- sh[[intersect(kinds, names(sh))[1]]]
  params <- first$parameters
  col <- function(k) {
    x <- sh[[k]]
    if (is.null(x)) return(rep(NA_real_, length(params)))
    # with some parameters lacking random effects, Monolix returns one wide row
    if ("shrinkage" %in% names(x)) x$shrinkage else unlist(unname(x[1, -1]))
  }
  data.frame(run = run, parameter = params, conditionalMean = col("conditionalMean"),
             conditionalMode = col("conditionalMode"), conditionalDist = col("conditionalDist"))
}
