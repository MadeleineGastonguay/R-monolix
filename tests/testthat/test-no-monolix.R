# Plotting and summary functions must work on saved results without Monolix, so someone
# without a Monolix licence can plot a project's output. This checks every function defined in
# these files, and every package function they call (directly or further down), for
# lixoftConnectors calls.
NO_MONOLIX_FILES <- c("plots_igs.R", "plots_ranges.R", "plots_model_building.R", "plots_fits.R",
                      "plots_model_comparison.R", "plots_profile.R", "diagnostics.R", "summaries.R", "model_building.R",
                      "selection.R", "sampling.R", "schema.R", "report.R")

# Names of the functions defined at the top level of an R file
functions_in <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  nm <- vapply(exprs, function(e) {
    if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) as.character(e[[2]]) else NA_character_
  }, character(1))
  nm[!is.na(nm)]
}

# Package functions reachable from `fns`, following calls inside the package
reachable <- function(fns, ns) {
  pkg_fns <- ls(ns, all.names = TRUE)
  pkg_fns <- pkg_fns[vapply(pkg_fns, function(f) is.function(get(f, ns)), logical(1))]
  seen <- character()
  todo <- fns
  while (length(todo)) {
    f <- todo[1]; todo <- todo[-1]
    if (f %in% seen) next
    seen <- c(seen, f)
    called <- intersect(all.names(body(get(f, ns))), pkg_fns)
    todo <- c(todo, setdiff(called, seen))
  }
  seen
}

test_that("plotting and summary functions don't call lixoftConnectors", {
  r_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "package R/ sources not available")
  ns <- asNamespace("monolixPipelineR")
  fns <- unlist(lapply(file.path(r_dir, NO_MONOLIX_FILES), functions_in))
  expect_gt(length(fns), 30)  # the file list is being read
  all_fns <- reachable(fns, ns)
  uses <- all_fns[vapply(all_fns, function(f) {
    any(grepl("lixoftConnectors::|require_monolix", deparse(get(f, ns))))
  }, logical(1))]
  expect(length(uses) == 0,
         paste0("These functions call Monolix (directly or through another package function), ",
                "so they can't be used on saved results without it: ", paste(uses, collapse = ", ")))
})

test_that("the check catches a Monolix call further down", {
  ns <- new.env()
  ns$plot_a <- function(x) helper_b(x)
  ns$helper_b <- function(x) lixoftConnectors::getPopulationParameterInformation()
  ns$plot_c <- function(x) x
  expect_setequal(reachable("plot_a", ns), c("plot_a", "helper_b"))
  expect_equal(reachable("plot_c", ns), "plot_c")
})
