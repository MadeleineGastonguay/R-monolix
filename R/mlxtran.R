# Reading the structural model and data files without Monolix, so configs can be
# drafted and checked on any machine.

#' Read the parts of an Mlxtran structural model that the pipeline needs
#'
#' Pure text parsing; Monolix is not needed. Comments (`;`) are ignored.
#'
#' @param model_file Path to the Mlxtran model file.
#' @return A list with `parameters` (individual parameters from `input = {...}`,
#'   excluding inputs declared as regressors), `regressors`, `outputs`
#'   (`output = {...}`), `tables` (`table = {...}`) and `rescalings` (inputs
#'   only multiplied or divided by a constant in the EQUATION block, see
#'   [guess_report_transforms()]).
#' @export
parse_mlxtran <- function(model_file) {
  if (!file.exists(model_file)) stop("Model file not found: ", model_file, call. = FALSE)
  lines <- readLines(model_file, warn = FALSE)
  lines <- sub(";.*$", "", lines)               # drop comments
  text <- paste(lines, collapse = "\n")

  braced <- function(key) {
    m <- regmatches(text, regexec(paste0("(?m)^\\s*", key, "\\s*=\\s*\\{([^}]*)\\}"), text, perl = TRUE))[[1]]
    if (length(m) == 0) {
      # single value without braces, e.g. `output = Cc`
      m <- regmatches(text, regexec(paste0("(?m)^\\s*", key, "\\s*=\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*$"), text, perl = TRUE))[[1]]
      if (length(m) == 0) return(character())
    }
    items <- trimws(strsplit(m[2], ",")[[1]])
    items[nzchar(items)]
  }

  inputs <- braced("input")
  if (length(inputs) == 0) stop("No `input = {...}` found in ", model_file, call. = FALSE)
  regressors <- inputs[vapply(inputs, function(p) {
    grepl(paste0("(?m)^\\s*", p, "\\s*=\\s*\\{\\s*use\\s*=\\s*regressor\\s*\\}"), text, perl = TRUE)
  }, logical(1))]

  outputs <- braced("output")
  if (length(outputs) == 0) stop("No `output = ...` found in ", model_file, call. = FALSE)

  parameters <- setdiff(inputs, regressors)
  list(
    parameters = parameters,
    regressors = regressors,
    outputs = outputs,
    tables = braced("table"),
    rescalings = find_rescalings(lines, parameters)
  )
}

# Equation lines of the form `new = param * c`, `new = c * param` or `new = param / c`,
# where c is a numeric constant such as 10^5 or 1e-7. First definition per parameter wins.
find_rescalings <- function(lines, parameters) {
  out <- list()
  const_ok <- function(s) grepl("^[0-9.eE+^() -]+$", s)
  eval_const <- function(s) {
    v <- tryCatch(eval(parse(text = s), envir = baseenv()), error = function(e) NA)
    if (is.numeric(v) && length(v) == 1 && is.finite(v) && v != 0) v else NA
  }
  for (ln in lines) {
    m <- regmatches(ln, regexec("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+?)\\s*$", ln))[[1]]
    if (length(m) == 0) next
    lhs <- m[2]; rhs <- m[3]
    hit <- NULL
    if (grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*\\*", rhs)) {
      p <- sub("\\s*\\*.*$", "", rhs); cst <- sub("^[^*]*\\*\\s*", "", rhs)
      if (p %in% parameters && const_ok(cst)) hit <- list(p, cst, "*")
    }
    if (is.null(hit) && grepl("\\*\\s*[A-Za-z_][A-Za-z0-9_]*$", rhs)) {
      p <- sub("^.*\\*\\s*", "", rhs); cst <- sub("\\s*\\*[^*]*$", "", rhs)
      if (p %in% parameters && const_ok(cst)) hit <- list(p, cst, "*")
    }
    if (is.null(hit) && grepl("^[A-Za-z_][A-Za-z0-9_]*\\s*/", rhs)) {
      p <- sub("\\s*/.*$", "", rhs); cst <- sub("^[^/]*/\\s*", "", rhs)
      if (p %in% parameters && const_ok(cst)) hit <- list(p, cst, "/")
    }
    if (is.null(hit) || hit[[1]] %in% names(out) || lhs == hit[[1]]) next
    if (is.na(eval_const(hit[[2]]))) next
    out[[hit[[1]]]] <- list(derived = lhs, expr = paste("x", hit[[3]], trimws(hit[[2]])))
  }
  out
}

#' Read a Monolix data file
#'
#' Guesses the delimiter (tab, comma, semicolon or whitespace) from the header
#' line. Monolix's `.` for missing values is read as `NA`.
#'
#' @param data_file Path to the data file.
#' @return A data frame.
#' @export
read_monolix_data <- function(data_file) {
  if (!file.exists(data_file)) stop("Data file not found: ", data_file, call. = FALSE)
  header <- readLines(data_file, n = 1, warn = FALSE)
  counts <- c("\t" = lengths(regmatches(header, gregexpr("\t", header))),
              "," = lengths(regmatches(header, gregexpr(",", header))),
              ";" = lengths(regmatches(header, gregexpr(";", header))))
  if (max(counts) > 0) {
    as.data.frame(readr::read_delim(data_file, delim = names(which.max(counts)), na = c("", "NA", "."),
                                    show_col_types = FALSE, progress = FALSE))
  } else {
    utils::read.table(data_file, header = TRUE, na.strings = c("NA", "."), stringsAsFactors = FALSE)
  }
}

# Monolix column types the pipeline understands
MONOLIX_HEADER_TYPES <- c("id", "time", "observation", "obsid", "cens", "limit", "catcov", "contcov",
                          "regressor", "ignore", "amount", "admid", "rate", "tinf", "mdv", "evid",
                          "addl", "ii", "ss", "occ")

#' Guess the Monolix column type of each data column
#'
#' Rules, in order: `id` and `time` by name; `cens` for 0/1 columns named like
#' cens/censor/blq/bql; `observation` by name (dv, y, obs, ...), otherwise the
#' numeric column that varies most within individuals; columns constant within
#' every individual become `catcov` (10 or fewer levels, or non-numeric) or
#' `contcov`; anything else is `ignore` and flagged for you to decide.
#'
#' @param data Data frame as read by [read_monolix_data()].
#' @return Data frame with `column`, `type` and `reason`.
#' @export
guess_header_types <- function(data) {
  cols <- names(data)
  type <- setNames(rep(NA_character_, length(cols)), cols)
  reason <- setNames(rep("", length(cols)), cols)
  lc <- tolower(cols)

  id <- cols[lc %in% c("id", "subject", "subj", "usubjid", "patient", "animal")][1]
  if (is.na(id)) { id <- cols[1]; reason[id] <- "no column named id; using the first column" } else reason[id] <- "named like an id"
  type[id] <- "id"

  tm <- cols[lc %in% c("time", "t", "day", "days", "tad", "hours", "weeks")]
  tm <- setdiff(tm, id)[1]
  if (!is.na(tm)) { type[tm] <- "time"; reason[tm] <- "named like time" }

  ids <- data[[id]]
  varies_within_id <- function(x) {
    any(tapply(x, ids, function(v) length(unique(v[!is.na(v)])) > 1), na.rm = TRUE)
  }
  is_binary <- function(x) all(x[!is.na(x)] %in% c(0, 1, -1))

  for (cc in cols[is.na(type)]) {
    if (grepl("cens|blq|bql", tolower(cc)) && is_binary(data[[cc]])) {
      type[cc] <- "cens"; reason[cc] <- "0/1 column named like censoring"
    }
  }

  obs <- cols[is.na(type) & lc %in% c("dv", "y", "obs", "observation", "conc", "vl", "lvl", "response")]
  if (length(obs) > 0) {
    type[obs] <- "observation"; reason[obs] <- "named like an observation"
  } else {
    cand <- cols[is.na(type)]
    cand <- cand[vapply(cand, function(cc) is.numeric(data[[cc]]) && !is_binary(data[[cc]]) &&
                          varies_within_id(data[[cc]]), logical(1))]
    if (length(cand) > 0) {
      n_unique <- vapply(cand, function(cc) length(unique(data[[cc]])), numeric(1))
      pick <- cand[which.max(n_unique)]
      type[pick] <- "observation"; reason[pick] <- "numeric column varying most within individuals"
    }
  }

  for (cc in cols[is.na(type)]) {
    x <- data[[cc]]
    if (!varies_within_id(x)) {
      n_lev <- length(unique(x[!is.na(x)]))
      if (!is.numeric(x) || n_lev <= 10) {
        type[cc] <- "catcov"; reason[cc] <- paste0("constant within each id, ", n_lev, " levels")
      } else {
        type[cc] <- "contcov"; reason[cc] <- paste0("constant within each id, ", n_lev, " distinct values")
      }
    } else {
      type[cc] <- "ignore"; reason[cc] <- "varies within individuals but isn't the observation (regressor?) - decide"
    }
  }
  data.frame(column = cols, type = unname(type), reason = unname(reason), stringsAsFactors = FALSE)
}

#' Guess report transforms from a parsed structural model
#'
#' An input that the model only rescales by a constant (e.g. `NP = NPs*10^5`)
#' is reported on the rescaled scale: `NPs_pop` is shown as `NP_pop = x * 10^5`.
#'
#' @param parsed Output of [parse_mlxtran()].
#' @return Named list, one entry per `<param>_pop` with `label` and `expr`.
#' @export
guess_report_transforms <- function(parsed) {
  out <- lapply(parsed$rescalings, function(r) list(label = paste0(r$derived, "_pop"), expr = r$expr))
  if (length(out)) names(out) <- paste0(names(parsed$rescalings), "_pop")
  out
}

#' Names Monolix gives error-model parameters
#'
#' With one observation the error parameters are `a`, `b` (and `c`); with
#' several they're numbered by observation in the order given: `a1`, `b1`,
#' `a2`, ... (checked on Monolix 2024R1).
#'
#' @param error_model Named list or vector, observation column -> error model
#'   (`model.error_model` in the config).
#' @return Character vector of parameter names.
#' @export
error_param_names <- function(error_model) {
  em <- unlist(error_model)
  if (length(em) == 1) return(error_params_for(em[[1]]))
  unlist(lapply(seq_along(em), function(i) {
    p <- error_params_for(em[[i]])
    if (length(p)) paste0(p, i) else character()
  }))
}

# Names Monolix reserves for error parameters with n observations; a structural parameter
# with one of these names makes newProject() fail ("Name 'a' is already used as error
# parameter"). Checked for `a`; b and c are reserved the same way to be safe.
reserved_error_names <- function(n_observations) {
  if (n_observations <= 1) c("a", "b", "c") else paste0(rep(c("a", "b", "c"), each = n_observations), seq_len(n_observations))
}

# Error-model parameters Monolix creates for one observation
error_params_for <- function(error_model) {
  switch(error_model,
         constant = "a",
         proportional = "b",
         combined1 = c("a", "b"),
         combined2 = c("a", "b"),
         character())
}

#' Initial-value range spanning a fold change around a rough value
#'
#' @param x Rough value of the parameter.
#' @param fold Fold change on each side.
#' @return Named numeric vector `c(min = x / fold, max = x * fold)`.
#' @export
range_around <- function(x, fold = 100) {
  stopifnot(is.numeric(x), x > 0, fold > 1)
  c(min = x / fold, max = x * fold)
}
