dat <- data.frame(ID = rep(1:4, each = 2), time = rep(0:1, 4),
                  tlr7 = rep(c(0, 1, 0, 1), each = 2), vac = rep(c(0, 0, 1, 1), each = 2))

test_that("report_groups keys groups by stratify_by values and applies labels", {
  g <- report_groups(dat, id_column = "ID", stratify_by = c("tlr7", "vac"),
                     group_labels = list("0.0" = "control", "1.1" = "tlr7 + vac"))
  expect_equal(g$id, 1:4)
  expect_equal(g$group_key, c("0.0", "1.0", "0.1", "1.1"))
  expect_equal(as.character(g$group), c("control", "tlr7=1, vac=0", "tlr7=0, vac=1", "tlr7 + vac"))
  # levels follow the key order
  expect_equal(levels(g$group), c("control", "tlr7=0, vac=1", "tlr7=1, vac=0", "tlr7 + vac"))
})

test_that("report_groups without stratify_by puts everyone in one group", {
  g <- report_groups(dat, id_column = "ID")
  expect_equal(nrow(g), 4)
  expect_equal(as.character(unique(g$group)), "All")
})

test_that("report_groups complains about missing or time-varying columns", {
  expect_error(report_groups(dat, id_column = "ID", stratify_by = "ab"), "not in the data: ab")
  expect_error(report_groups(dat, id_column = "id"), "No id column 'id'")
  dat2 <- dat; dat2$tlr7[2] <- 1
  expect_error(report_groups(dat2, id_column = "ID", stratify_by = "tlr7"), "change within an individual")
})

test_that("transform_estimates rescales estimates and limits and relabels", {
  est <- data.frame(parameter = c("NPs_pop", "p_pop"), estimate = c(2, 3), P2.5 = c(1, 2), P97.5 = c(4, 5))
  out <- transform_estimates(est, list(NPs_pop = list(label = "NP_pop", expr = "x * 10^5")))
  expect_equal(out$reported_parameter, c("NP_pop", "p_pop"))
  expect_equal(out$reported_estimate, c(2e5, 3))
  expect_equal(out$reported_P2.5, c(1e5, 2))
  expect_equal(out$reported_P97.5, c(4e5, 5))
  # without limits or transforms
  out2 <- transform_estimates(est[, 1:2])
  expect_equal(out2$reported_estimate, est$estimate)
  expect_false("reported_P2.5" %in% names(out2))
})

test_that("error-model parameter names are recognised for the report", {
  expect_equal(is_error_param_name(c("a", "b", "a1", "c2", "lambda_vac_1", "a_pop", "beta_a_vac")),
               c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("read_indfits reads Monolix's individual fits chart data", {
  res <- tempfile("res"); dir.create(res)
  dir.create(file.path(res, "ChartsData", "IndividualFits"), recursive = TRUE)
  expect_null(read_indfits(res))
  # columns as Monolix 2024R1 writes them; `pop` uses the median covariates
  write.csv(data.frame(ID = c(1, 1, 2), time = c(0, 0.5, 0), pop = 10:12, popPred = 7:9,
                       indivPredMean = 1:3, indivPredMode = 4:6),
            file.path(res, "ChartsData", "IndividualFits", "vl_fits.txt"), row.names = FALSE)
  f <- read_indfits(res, "vl")
  expect_equal(names(f), c("id", "time", "popPred", "indivPred"))
  expect_equal(f$indivPred, 4:6)
  expect_equal(f$popPred, 7:9)
  expect_equal(attr(f, "indiv_column"), "indivPredMode")
  # several observations: picks the one asked for
  write.csv(data.frame(ID = 1, time = 0, indivPredMode = 0, popPred = 0),
            file.path(res, "ChartsData", "IndividualFits", "y2_fits.txt"), row.names = FALSE)
  expect_equal(read_indfits(res, "vl")$indivPred, 4:6)
  expect_null(read_indfits(res))
})

test_that("read_model_variables reads table variables from tables.txt in long format", {
  res <- tempfile("res"); dir.create(res)
  expect_null(read_model_variables(res))
  write.csv(data.frame(id = c(1, 1, 2), time = c(0, 4, 0),
                       V_mean = 1:3, T_mean = 4:6, V_mode = 7:9, T_mode = 10:12,
                       V_popPred = 0, T_popPred = 0, V_popPred_medianCOV = -1, T_popPred_medianCOV = -1),
            file.path(res, "tables.txt"), row.names = FALSE)
  v <- read_model_variables(res)
  expect_equal(levels(v$variable), c("V", "T"))
  expect_equal(v$value, c(7:9, 10:12))
  expect_equal(v$id, c(1, 1, 2, 1, 1, 2))
  expect_equal(attr(v, "prediction"), "mode")
  expect_equal(unique(read_model_variables(res, "popPred")$value), 0)
  # no _mode columns: falls back to the conditional mean
  write.csv(data.frame(id = 1, time = 0, V_mean = 1), file.path(res, "tables.txt"), row.names = FALSE)
  expect_equal(attr(read_model_variables(res), "prediction"), "mean")
})

test_that("simulation_grid spans each individual's own time range", {
  d <- data.frame(ID = c(1, 1, 1, 2, 2), t = c(0, 5, 10, 2, 4))
  g <- simulation_grid(d, "ID", "t", n_points = 3)
  expect_equal(g$id, c(1, 1, 1, 2, 2, 2))
  expect_equal(g$time, c(0, 5, 10, 2, 3, 4))
})

test_that("read_simulated_variables reads simulate_model_variables() output", {
  res <- tempfile("res"); dir.create(res)
  expect_null(read_simulated_variables(res))
  f <- simulated_variables_file(res); dir.create(dirname(f))
  write.csv(data.frame(id = 1, time = c(0, 1, 0, 1), variable = c("V", "V", "T", "T"), value = 1:4,
                       parameters = "mlx_EBEs"), f, row.names = FALSE)
  v <- read_simulated_variables(res)
  expect_equal(levels(v$variable), c("V", "T"))
  expect_equal(v$value, 1:4)
  expect_equal(attr(v, "prediction"), "mode")
})
