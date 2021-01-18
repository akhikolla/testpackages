#' Read a hierarchical fixed width data file
#'
#' Analogous to [readr::read_fwf()] but allowing for
#' hierarchical fixed width data files (where the data file has rows of
#' different record types, each with their own variables and column
#' specifications). `hipread_long()` reads hierarchical data into "long"
#' format, meaning that there is one row per observation, and variables
#' that don't apply to the current observation receive missing values.
#' Alternatively, `hipread_list()` reads hierarchical data into "list"
#' format, which returns a list that has one data.frame per record type.
#'
#' @param file A filename
#' @param var_info Variable information, specified by either [`hip_fwf_positions()`]
#'   or `hip_fwf_widths()`. For hierarchical data files, there should be a named list,
#'   where the name is the value indicated by the record type variable and there is
#'   one variable information per record type.
#' @param rt_info A record type information object, created by [`hip_rt()`], which
#'   contains information about the location of the record type variable that
#'   defines the record type for each observation. The default contains width
#'   0, which indicates that there the data is rectangular and does not have
#'   a record type variable.
#' @param skip Number of lines to skip at the start of the data (defaults to 0).
#' @param n_max Maximum number of lines to read. Negative numbers (the default)
#'   reads all lines.
#' @param encoding (Defaults to UTF-8) A string indicating what encoding to use
#'   when reading the data, but like readr, the data will always be converted to
#'   UTF-8 once it is imported. Note that UTF-16 and UTF-32 are not supported for
#'   non-character columns.
#' @param compression If `NULL`, guesses the compression from the
#'   file extension (if extension is "gz" uses gzip, otherwise
#'   treats as plain text), can specify it with a string ("txt"
#'   indicates plain text and "gz" for gzip).
#' @param progress A logical indicating whether progress should be
#'   displayed on the screen, defaults to showing progress unless
#'   the current context is non-interactive or in a knitr document or
#'   if the user has turned off readr's progress by default using
#'   the option `options("readr.show_progress")`.
#'
#' @return A `tbl_df` data frame
#' @export
#'
#' @examples
#' # Read an example hierarchical data.frame into long format
#' data <- hipread_long(
#'   hipread_example("test-basic.dat"),
#'   list(
#'     H = hip_fwf_positions(
#'       c(1, 2, 5, 8),
#'       c(1, 4, 7, 10),
#'       c("rt", "hhnum", "hh_char", "hh_dbl"),
#'       c("c", "i", "c", "d")
#'     ),
#'     P = hip_fwf_widths(
#'       c(1, 3, 1, 3, 1),
#'       c("rt", "hhnum",  "pernum", "per_dbl", "per_mix"),
#'       c("c", "i", "i", "d", "c")
#'     )
#'   ),
#'   hip_rt(1, 1)
#' )
#'
#' # Read an example hierarchical data.frame into list format
#' data <- hipread_list(
#'   hipread_example("test-basic.dat"),
#'   list(
#'     H = hip_fwf_positions(
#'       c(1, 2, 5, 8),
#'       c(1, 4, 7, 10),
#'       c("rt", "hhnum", "hh_char", "hh_dbl"),
#'       c("c", "i", "c", "d")
#'     ),
#'     P = hip_fwf_widths(
#'       c(1, 3, 1, 3, 1),
#'       c("rt", "hhnum",  "pernum", "per_dbl", "per_mix"),
#'       c("c", "i", "i", "d", "c")
#'     )
#'   ),
#'   hip_rt(1, 1)
#' )
#'
#' # Read a rectangular data.frame
#' data_rect <- hipread_long(
#'   hipread_example("test-basic.dat"),
#'   hip_fwf_positions(
#'     c(1, 2),
#'     c(1, 4),
#'     c("rt", "hhnum"),
#'     c("c", "i")
#'   )
#' )
hipread_long <- function(
  file, var_info, rt_info = hip_rt(1, 0), compression = NULL,
  skip = 0, n_max = -1, encoding = "UTF-8", progress = show_progress()
) {
  file <- check_file(file)
  isgzipped <- is_gzip_compression(compression, file)
  var_info <- add_level_to_rect(var_info)
  var_names <- get_var_names(var_info)
  var_pos_info <- get_var_pos(var_info, var_names)
  var_types <- get_var_types(var_info, var_names)
  var_opts <- get_var_opts(var_info, var_names)
  skip <- check_skip(skip)
  n_max <- check_n_max(n_max)
  if (is.null(encoding)) encoding <- "UTF-8"

  read_long(
    file, var_names, var_types, rt_info, var_pos_info,
    var_opts, skip, n_max, isgzipped, encoding, progress
  )
}


#' @rdname hipread_long
#' @export
hipread_list <- function(
  file, var_info, rt_info = hip_rt(1, 0), compression = NULL,
  skip = 0, n_max = -1, encoding = "UTF-8", progress = show_progress()
) {
  file <- check_file(file)
  isgzipped <- is_gzip_compression(compression, file)
  var_info <- add_level_to_rect(var_info)
  var_names <- get_vinfo_col_as_list(var_info, "col_names")
  var_pos_info <- get_var_pos(var_info)
  var_types <- get_vinfo_col_as_list(var_info, "col_types")
  var_opts <- get_var_opts_list(var_info)
  skip <- check_skip(skip)
  n_max <- check_n_max(n_max)
  if (is.null(encoding)) encoding <- "UTF-8"

  read_list(
    file, var_names, var_types, rt_info, var_pos_info,
    var_opts, skip, n_max, isgzipped, encoding, progress
  )
}
