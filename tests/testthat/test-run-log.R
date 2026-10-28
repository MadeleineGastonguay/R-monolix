test_that("run logs detect changed settings and changed data, but ignore labels", {
  path <- make_config(list("  fixed: ~" = "  fixed: {lambda_pop: 50}"))
  cfg <- read_analysis_config(path)
  out <- file.path(dirname(path), "results", "step2")
  expect_false(check_run_log(cfg, 2, out))
  write_run_log(cfg, 2, out)
  expect_true(file.exists(file.path(out, "run_log.yaml")))
  expect_true(check_run_log(cfg, 2, out))

  # labels, plot options and rationales don't matter
  txt <- readLines(path)
  txt <- sub("  fixed_rationale: ~", "  fixed_rationale: profile likelihood minimum", txt, fixed = TRUE)
  txt <- sub("remove_bicc_outliers: false", "remove_bicc_outliers: true", txt, fixed = TRUE)
  writeLines(txt, path)
  expect_true(check_run_log(read_analysis_config(path), 2, out))

  # a fitting setting does
  writeLines(sub("  wide_nits: 100", "  wide_nits: 50", txt, fixed = TRUE), path)
  expect_error(check_run_log(read_analysis_config(path), 2, out), "step1_igs.wide_nits: 100 -> 50")

  # so does the scenario
  writeLines(sub("  linearization: false", "  linearization: true", txt, fixed = TRUE), path)
  expect_error(check_run_log(read_analysis_config(path), 2, out), "scenario.linearization: FALSE -> TRUE")

  # so does the data file's content
  writeLines(txt, path)
  data_file <- file.path(dirname(path), "synthetic_vl_data.txt")
  write(readLines(data_file)[1:50], data_file)
  expect_error(check_run_log(read_analysis_config(path), 2, out), "files.data")
})

test_that("the run log records seeds, hashes and versions", {
  cfg <- read_analysis_config(make_config())
  out <- tempfile("log")
  log <- yaml::read_yaml(write_run_log(cfg, 1, out))
  expect_equal(log$seeds$monolix, 123456)
  expect_equal(log$seeds$r, 1)
  expect_equal(nchar(log$files$data), 32)
  expect_false("decisions" %in% names(log$settings))
  expect_true(length(log$session_info) > 5)
})
