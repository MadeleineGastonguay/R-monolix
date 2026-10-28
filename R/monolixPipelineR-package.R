#' monolixPipelineR: multi-start fitting and model comparison for Monolix
#'
#' Monolix itself is reached through `lixoftConnectors`, which is not on CRAN
#' and is only needed when fitting. Call [init_monolix()] once per session
#' before any fitting function.
#'
#' @keywords internal
#' @import dplyr tidyr ggplot2 purrr tibble stringr readr forcats patchwork
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom gridExtra grid.arrange
#' @importFrom tune coord_obs_pred
#' @importFrom ggpmisc use_label stat_correlation
#' @importFrom rlang sym
#' @importFrom stats IQR quantile cor cor.test median sd runif setNames
"_PACKAGE"

# Silence R CMD check notes for column names used in dplyr/ggplot2 code
utils::globalVariables(c(
  ".", "run", "bicc", "estimate", "parameter", "param", "initialValue", "name",
  "value", "iteration", "rses", "cumprob", "w", "t", "best", "phase", "upper",
  "lower", "n_covar", "BICc", "threshold", "treatment", "final_mod", "label",
  "range", "step", "color", "parameter1", "parameter2", "pval", "CV", "y",
  "id", "time", "method", "initial_value", "rse", "profiled_parameter",
  "profiled_value", "initialValue"
))
