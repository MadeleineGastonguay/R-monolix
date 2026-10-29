# Simulating the fitted individuals with Simulx, for the report's model-variable curves.
# Monolix only writes table variables (tables.txt) at the observation times.

#' Simulate model variables for each fitted individual with Simulx
#'
#' Imports the fitted Monolix project into Simulx and simulates each individual
#' from their estimated parameters (EBEs, i.e. conditional modes; conditional
#' means if there are no EBEs) on a fine time grid over their own observed
#' time range (dose records included). Each individual's doses from the data
#' are kept (the `mlx_` treatment elements made by importing the project).
#' lixoftConnectors is switched to Simulx and back to Monolix afterwards (also
#' on error), so any loaded Monolix project is closed.
#'
#' @param project_file Path to a Monolix project that has been run.
#' @param variables Model variables to simulate. `NULL`: the model's
#'   `table = {...}` variables, or its outputs when there is no table.
#' @param n_points Number of time points per individual.
#' @param monolix_path Monolix install path for re-initializing; `""` uses
#'   the path lixoftConnectors was loaded with.
#' @return The path of the written file
#'   (`<results>/ModelVariables/simulated_model_variables.txt`), invisibly.
#'   Read it with [read_simulated_variables()].
#' @export
simulate_model_variables <- function(project_file, variables = NULL, n_points = 250,
                                     monolix_path = Sys.getenv("MONOLIX_PATH")) {
  require_monolix()
  project_file <- normalizePath(project_file, mustWork = TRUE)
  results_dir <- sub("\\.mlxtran$", "", project_file)

  # Data and model from the Monolix project
  if (isFALSE(lixoftConnectors::loadProject(projectFile = project_file))) {
    stop("Monolix could not load ", project_file, call. = FALSE)
  }
  data_info <- lixoftConnectors::getData()
  if (is.null(variables)) {
    parsed <- tryCatch(parse_mlxtran(lixoftConnectors::getStructuralModel()), error = function(e) NULL)
    variables <- if (length(parsed$tables)) parsed$tables else parsed$outputs
  }
  if (!length(variables)) {
    stop("No variables to simulate: pass `variables`, or list them in the model's `table = {...}`.", call. = FALSE)
  }
  data <- read_monolix_data(data_info$dataFile)
  grid <- simulation_grid(data, id_column = data_info$header[data_info$headerTypes == "id"],
                          time_column = data_info$header[data_info$headerTypes == "time"],
                          n_points = n_points)
  # Subject-specific output times have to come from a file with the id column first
  grid_file <- tempfile("simulx_output_times", fileext = ".csv")
  utils::write.csv(grid, grid_file, row.names = FALSE)

  # Simulx, switching back to Monolix when done
  lixoftConnectors::initializeLixoftConnectors(software = "simulx", path = monolix_path, force = TRUE)
  on.exit(lixoftConnectors::initializeLixoftConnectors(software = "monolix", path = monolix_path, force = TRUE),
          add = TRUE)
  lixoftConnectors::importProject(project_file)
  indiv <- intersect(c("mlx_EBEs", "mlx_CondMean"), names(lixoftConnectors::getIndividualElements()))[1]
  if (is.na(indiv)) {
    stop("No individual estimates (EBEs or conditional means) in ", project_file,
         "; run the conditional mode task first.", call. = FALSE)
  }
  out_elements <- paste0("report_", variables)
  for (i in seq_along(variables)) {
    lixoftConnectors::defineOutputElement(name = out_elements[i],
                                          element = list(data = grid_file, output = variables[i]))
  }
  # doses: importProject() makes one mlx_ treatment element per administration type in the data;
  # set them explicitly so the simulation keeps each individual's doses
  treatments <- grep("^mlx_", names(tryCatch(lixoftConnectors::getTreatmentElements(), error = function(e) list())),
                     value = TRUE)
  group <- lixoftConnectors::getGroups()[[1]]$name
  lixoftConnectors::setGroupElement(group = group, elements = c(indiv, treatments, out_elements))
  lixoftConnectors::setSamplingMethod("keepOrder")
  lixoftConnectors::setGroupSize(group, length(unique(grid$id)))
  lixoftConnectors::runSimulation()
  res <- lixoftConnectors::getSimulationResults()$res

  sims <- do.call(rbind, lapply(variables, function(v) {
    r <- res[[v]]
    if (is.null(r)) stop("No Simulx result for ", v, call. = FALSE)
    # original_id: the id in the data (Simulx renumbers ids when tables have their own)
    id <- if ("original_id" %in% names(r)) r$original_id else r$id
    data.frame(id = id, time = r$time, variable = v, value = r[[v]], parameters = indiv)
  }))
  out_file <- simulated_variables_file(results_dir)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(sims, out_file, row.names = FALSE)
  invisible(out_file)
}
