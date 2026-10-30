# The README's example analysis.yaml must stay a valid config. It's checked against a small
# made-up model and data set that use its names (param1-3, y, cov1-3).
test_that("the README's example analysis.yaml passes validate_config() at every step", {
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not found (e.g. during R CMD check)")
  lines <- readLines(readme, warn = FALSE)
  start <- grep("^## Example analysis.yaml", lines)
  expect_length(start, 1)
  open <- start + grep("^``` ?yaml", lines[-seq_len(start)])[1]
  close <- open + grep("^```", lines[-seq_len(open)])[1]
  yaml_txt <- lines[(open + 1):(close - 1)]

  dir <- withr_tempdir()
  writeLines(c("DESCRIPTION: example model", "", "[LONGITUDINAL]", "input = {param1, param2, param3}", "",
               "EQUATION:", "f = param1 * exp(-param3 * t) + param2 / 100", "", "OUTPUT:", "output = f"),
             file.path(dir, "model.txt"))
  d <- expand.grid(time = 0:4, id = 1:20)[, c("id", "time")]
  d$y <- 1 + d$time
  d$cov1 <- d$id %% 2
  d$cov2 <- (d$id %/% 2) %% 2
  d$cov3 <- 60 + d$id
  utils::write.csv(d, file.path(dir, "data.csv"), row.names = FALSE)
  writeLines(yaml_txt, file.path(dir, "analysis.yaml"))

  cfg <- read_analysis_config(file.path(dir, "analysis.yaml"))
  expect_length(attr(cfg, "unknown_keys"), 0)
  for (step in 1:4) expect_invisible(validate_config(cfg, step))

  # the logit-normal parameter and the report's group keys work as the README says
  b <- parameter_bounds(cfg$model$distributions, c("param1", "param2", "param3"))
  expect_equal(b$upper[b$parameter == "param2_pop"], 100)
  g <- report_groups(read_monolix_data(cfg$data$file), "id", cfg$report$stratify_by, cfg$report$group_labels)
  expect_equal(as.character(g$group[g$group_key == "1.0"][1]), "treated")
})
