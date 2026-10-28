# Connecting to Monolix. Nothing here runs when the package loads.

#' Is lixoftConnectors installed?
#'
#' @return `TRUE` if the lixoftConnectors package can be loaded.
#' @export
monolix_available <- function() {
  requireNamespace("lixoftConnectors", quietly = TRUE)
}

# Stop with a clear message if a Monolix-dependent function is called without it.
require_monolix <- function() {
  if (!monolix_available()) {
    stop("lixoftConnectors is not installed. Install it from your Monolix installation ",
         "and call init_monolix() before fitting.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Initialize the Monolix connection
#'
#' @param path Path to the Monolix installation. Defaults to the
#'   `MONOLIX_PATH` environment variable so machine-specific paths stay out of
#'   the analysis config.
#' @param check If `TRUE`, run [check_installation()] afterwards and print the
#'   demo estimates.
#' @return Invisibly, the demo estimates (or `TRUE` when `check = FALSE`).
#' @export
init_monolix <- function(path = Sys.getenv("MONOLIX_PATH"), check = TRUE) {
  require_monolix()
  if (!nzchar(path)) {
    stop("No Monolix path given. Set MONOLIX_PATH (e.g. in ~/.Renviron) or pass `path`.",
         call. = FALSE)
  }
  lixoftConnectors::initializeLixoftConnectors(software = "monolix", path = path)
  if (!check) return(invisible(TRUE))
  est <- check_installation()
  cat("Monolix works: theophylline demo estimates\n")
  print(est)
  invisible(est)
}

#' Check that Monolix runs, using the theophylline demo project
#'
#' @export
check_installation <- function() {
  require_monolix()
  demo_path <- lixoftConnectors::getDemoPath()
  lixoftConnectors::loadProject(file.path(demo_path, "1.creating_and_using_models",
                                          "1.1.libraries_of_models",
                                          "theophylline_project.mlxtran"))
  lixoftConnectors::runScenario()
  lixoftConnectors::getEstimatedPopulationParameters()
}
