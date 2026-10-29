fixture <- function(...) testthat::test_path("fixtures", ...)

# Draft a config (VACTLR7 model, synthetic data) in a temp folder, fill the required ranges, and apply
# extra text replacements. Returns the YAML path.
make_config <- function(replace = list(), dir = withr_tempdir()) {
  file.copy(fixture(c("synthetic_vl_data.txt", "2regime_model_delay.txt")), dir)
  path <- file.path(dir, "analysis.yaml")
  suppressMessages(draft_config(file.path(dir, "2regime_model_delay.txt"), file.path(dir, "synthetic_vl_data.txt"),
                                name = "study2", path = path))
  txt <- readLines(path)
  ranges <- c(lambda = "{min: 10, max: 1000, log10: true}", betas = "{min: 0.001, max: 100, log10: true}",
              ta = "{min: 0.01, max: 10, log10: true}", NPs = "{min: 0.001, max: 100, log10: true}",
              p = "{min: 0.01, max: 10, log10: true}", m = "{min: 0.01, max: 100, log10: true}")
  for (k in names(ranges)) {
    txt <- sub(paste0("^    ", k, ": \\{min: ~, max: ~, log10: true\\}.*$"), paste0("    ", k, ": ", ranges[[k]]), txt)
  }
  for (old in names(replace)) {
    hit <- grepl(old, txt, fixed = TRUE)
    if (!any(hit)) stop("no line matches: ", old)
    txt[hit] <- replace[[old]]
  }
  writeLines(txt, path)
  path
}

withr_tempdir <- function() {
  d <- tempfile("cfg")
  dir.create(d)
  d
}
