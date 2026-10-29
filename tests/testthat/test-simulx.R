test_that("simulate_model_variables simulates table variables per individual and switches back to Monolix", {
  skip_if_not_installed("lixoftConnectors")
  proj_dir <- tempfile("proj"); dir.create(proj_dir)
  project <- file.path(proj_dir, "fit.mlxtran"); writeLines("", project)
  calls <- character()
  outputs <- list()
  local_mocked_bindings(
    loadProject = function(...) TRUE,
    getData = function() list(dataFile = test_path("fixtures", "synthetic_vl_data.txt"),
                              header = c("id", "time", "vl"), headerTypes = c("id", "time", "observation")),
    getStructuralModel = function() test_path("fixtures", "2regime_model_delay.txt"),
    initializeLixoftConnectors = function(software, ...) { calls <<- c(calls, software); TRUE },
    importProject = function(...) TRUE,
    getIndividualElements = function() list(mlx_PopIndiv = list(), mlx_EBEs = list()),
    getTreatmentElements = function() list(),  # no doses in this model
    defineOutputElement = function(name, element) outputs[[name]] <<- element,
    getGroups = function() list(list(name = "simulationGroup1")),
    setGroupElement = function(group, elements) calls <<- c(calls, paste(elements, collapse = ",")),
    setSamplingMethod = function(...) TRUE,
    setGroupSize = function(group, size) calls <<- c(calls, paste("size", size)),
    runSimulation = function() TRUE,
    getSimulationResults = function() {
      grid <- read.csv(outputs[[1]]$data)
      res <- lapply(setNames(nm = c("T", "V", "E", "I", "P")), function(v) {
        r <- data.frame(id = match(grid$id, unique(grid$id)), original_id = grid$id, time = grid$time)
        r[[v]] <- seq_len(nrow(grid))
        r
      })
      list(res = res)
    },
    .package = "lixoftConnectors"
  )
  out <- simulate_model_variables(project, n_points = 5)
  # the model's table = {T, V, E, I, P}, one output element each, on each individual's own range
  expect_equal(unname(vapply(outputs, `[[`, "", "output")), c("T", "V", "E", "I", "P"))
  grid <- read.csv(outputs[[1]]$data)
  expect_equal(names(grid), c("id", "time"))
  dat <- read_monolix_data(test_path("fixtures", "synthetic_vl_data.txt"))
  expect_equal(nrow(grid), 5 * length(unique(dat$id)))
  expect_equal(calls, c("simulx", "mlx_EBEs,report_T,report_V,report_E,report_I,report_P",
                        paste("size", length(unique(dat$id))), "monolix"))
  sims <- read_simulated_variables(file.path(proj_dir, "fit"))
  expect_equal(levels(sims$variable), c("T", "V", "E", "I", "P"))
  expect_equal(sort(unique(sims$id)), sort(unique(dat$id)))
  expect_equal(attr(sims, "prediction"), "mode")
})

test_that("simulate_model_variables switches back to Monolix when Simulx fails", {
  skip_if_not_installed("lixoftConnectors")
  proj_dir <- tempfile("proj"); dir.create(proj_dir)
  project <- file.path(proj_dir, "fit.mlxtran"); writeLines("", project)
  calls <- character()
  local_mocked_bindings(
    loadProject = function(...) TRUE,
    getData = function() list(dataFile = test_path("fixtures", "synthetic_vl_data.txt"),
                              header = c("id", "time", "vl"), headerTypes = c("id", "time", "observation")),
    getStructuralModel = function() test_path("fixtures", "2regime_model_delay.txt"),
    initializeLixoftConnectors = function(software, ...) { calls <<- c(calls, software); TRUE },
    importProject = function(...) TRUE,
    getIndividualElements = function() list(mlx_PopIndiv = list()),
    .package = "lixoftConnectors"
  )
  expect_error(simulate_model_variables(project), "No individual estimates")
  expect_equal(calls, c("simulx", "monolix"))
})

test_that("simulate_model_variables keeps the doses imported from the Monolix project", {
  skip_if_not_installed("lixoftConnectors")
  proj_dir <- tempfile("proj"); dir.create(proj_dir)
  project <- file.path(proj_dir, "fit.mlxtran"); writeLines("", project)
  data_file <- file.path(proj_dir, "pk.csv")
  write.csv(data.frame(ID = c(1, 1, 1, 2, 2), TIME = c(0, 2, 12, 0, 6), DOSE = c(25, NA, NA, 20, NA),
                       CONC = c(NA, 17, 9, NA, 12)), data_file, row.names = FALSE, na = ".")
  model <- file.path(proj_dir, "pk.txt")
  writeLines(c("[LONGITUDINAL]", "input = {V, k}", "PK:", "Cc = pkmodel(V, k)", "OUTPUT:", "output = Cc"), model)
  set <- NULL; outputs <- list()
  local_mocked_bindings(
    loadProject = function(...) TRUE,
    getData = function() list(dataFile = data_file, header = c("ID", "TIME", "DOSE", "CONC"),
                              headerTypes = c("id", "time", "amount", "observation")),
    getStructuralModel = function() model,
    initializeLixoftConnectors = function(...) TRUE,
    importProject = function(...) TRUE,
    getIndividualElements = function() list(mlx_EBEs = list()),
    getTreatmentElements = function() list(mlx_Adm1 = list(), my_regimen = list()),
    defineOutputElement = function(name, element) outputs[[name]] <<- element,
    getGroups = function() list(list(name = "simulationGroup1")),
    setGroupElement = function(group, elements) set <<- elements,
    setSamplingMethod = function(...) TRUE,
    setGroupSize = function(...) TRUE,
    runSimulation = function() TRUE,
    getSimulationResults = function() {
      grid <- read.csv(outputs[[1]]$data)
      list(res = list(Cc = data.frame(id = grid$id, time = grid$time, Cc = 1)))
    },
    .package = "lixoftConnectors"
  )
  simulate_model_variables(project, n_points = 4)
  expect_equal(set, c("mlx_EBEs", "mlx_Adm1", "report_Cc"))  # only the imported treatment
  grid <- read.csv(outputs[[1]]$data)
  expect_equal(min(grid$time[grid$id == 1]), 0)  # from the dose, not the first observation
})
