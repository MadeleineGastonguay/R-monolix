# Small shared helpers.

# One Monolix seed for every project so SAEM runs are reproducible. Not user-configurable.
MONOLIX_SEED <- 123456

# Folder where Monolix writes results (predictions.txt, tables.txt) for the loaded project.
# Replaces reading `out_folder`/`project_name` from the global environment.
project_results_dir <- function() {
  lixoftConnectors::getProjectSettings()$directory
}

#' Which population parameters are error-model parameters?
#'
#' Monolix names population parameters `<param>_pop`, random-effect SDs
#' `omega_<param>`, covariate effects `beta_<param>_<cov>` and correlations
#' `corr_...`; everything else (e.g. `a`, `b`, `a_y1`) belongs to the error model.
#'
#' @param names Character vector of population parameter names.
#' @return Logical vector, `TRUE` for error-model parameters.
#' @export
is_error_param <- function(names) {
  !grepl("_pop$|^omega_|^omega2_|^beta_|^corr_|^gamma_", names)
}
