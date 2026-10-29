# Column names used in every table the package produces or accepts.
#
#   parameter, run, estimate, rse, ofv, aic, bic, bicc   (estimates)
#   parameter, run, initial_value, method                (initial values)
#   model                                                (tables spanning several candidate models)
#   profiled_parameter, profiled_value                   (profile likelihood)
#
# Tables passed to Monolix itself keep Monolix's own names (`name`, `initialValue`).

# Old (SIV) column name -> package column name
OLD_COLUMN_NAMES <- c(rses = "rse", initialValue = "initial_value", am = "model")

#' Rename old column names to the package's names
#'
#' Converts tables that use older column names: `rses` becomes `rse`, `initialValue` becomes `initial_value` and `am`
#' becomes `model`. Old profile likelihood tables stored the profiled value in
#' a column named after the parameter (e.g. `param1_pop`); pass that name as
#' `profiled_parameter` to convert it to `profiled_parameter`/`profiled_value`.
#'
#' @param df A data frame.
#' @param profiled_parameter For old profile likelihood tables, the name of the
#'   profiled parameter.
#' @return `df` with package column names. Tables already in the package
#'   format are returned unchanged.
#' @export
standardize_columns <- function(df, profiled_parameter = NULL) {
  if (!is.data.frame(df)) return(df)
  old <- intersect(names(OLD_COLUMN_NAMES), names(df))
  clash <- OLD_COLUMN_NAMES[old][OLD_COLUMN_NAMES[old] %in% names(df)]
  if (length(clash) > 0) {
    stop("Table has both old and new names for: ", paste(clash, collapse = ", "), call. = FALSE)
  }
  names(df)[match(old, names(df))] <- unname(OLD_COLUMN_NAMES[old])

  if (!is.null(profiled_parameter)) {
    if (!profiled_parameter %in% names(df)) {
      stop("No column named '", profiled_parameter, "' to convert.", call. = FALSE)
    }
    if ("profiled_value" %in% names(df)) {
      stop("Table already has a profiled_value column.", call. = FALSE)
    }
    names(df)[names(df) == profiled_parameter] <- "profiled_value"
    df$profiled_parameter <- profiled_parameter
    df <- df[, c("profiled_parameter", "profiled_value", setdiff(names(df), c("profiled_parameter", "profiled_value")))]
  }
  df
}

#' Load saved initial guess sampling results
#'
#' Reads an `.RData` file of initial guess sampling results (e.g. `IGS.RData`,
#' `log10_IGS.RData`) and returns its objects as a list, with
#' column names converted by [standardize_columns()]. The file itself is not
#' changed.
#'
#' @param path Path to the `.RData` file.
#' @return Named list (e.g. `tabestimates`, `tabiters`, `phase`, `initialvalues`).
#' @export
load_igs_results <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  out <- mget(ls(env, sorted = TRUE), envir = env)
  lapply(out, standardize_columns)
}
