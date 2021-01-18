#' Specify column-specific options for hipread
#'
#' Specify column specifications analogous to `readr::fwf_positions()`.
#' However, unlike in readr, the column type information is specified
#' alongside the column positions and there are two extra options that
#' can be specified (`trim_ws` gives control over trimming whitespace
#' in character columns, and `imp_dec` allows for implicit decimals in
#' double columns).
#'
#' @param start,end A vector integers describing the start and end positions
#'   of each field
#' @param col_names A character vector of variable names
#' @param col_types A vector of column types (specified as either
#'   "c" or "character" for character, "d" or "double" for double and
#'   "i" or "integer" for integer).
#' @param trim_ws A logical vector, indicating whether to trim whitespace
#'   on both sides of character columns (Defaults to `TRUE`, ignored on
#'   non-character columns).
#' @param imp_dec An integer vector, indicating the number of implicit decimals
#'   on a double variable (Defaults to 0, ignored on non-double columns).
#'
#' @return A data.frame containing the column specifications
#' @export
#'
#' @examples
#' # 3 Columns, specified by position
#' hip_fwf_positions(
#'   c(1, 3, 7),
#'   c(2, 6, 10),
#'   c("Var1", "Var2", "Var3"),
#'   c("c", "i", "d")
#' )
#'
#' # The same 3 columns, specified by width
#' hip_fwf_widths(
#'   c(2, 4, 4),
#'   c("Var1", "Var2", "Var3"),
#'   c("c", "i", "d")
#' )
#'
hip_fwf_positions <- function(
  start, end, col_names, col_types, trim_ws = TRUE, imp_dec = 0
) {
  if (!is_integerish(start)) stop("fwf column start positions in must be integers.")
  if (!all(start > 0)) stop("fwf column start positions must be greater or equal to 1.")
  if (!is_integerish(end)) stop("fwf column end positions in must be integers.")
  if (!all(end >= start)) stop("fwf column end positions must be at least as large as the start positions")
  if (!is.logical(trim_ws)) stop("trim_ws must be a logical vector (TRUE/FALSE).")
  if (!is_integerish(imp_dec)) stop("imp_dec must be an integer.")

  out <- tibble::tibble(
    start = start - 1,
    end = end,
    col_names = col_names,
    col_types = standardize_col_types(col_types),
    trim_ws = trim_ws,
    imp_dec = imp_dec
  )

  class(out) <- c("hip_pos", class(out))
  out
}

#' @rdname hip_fwf_positions
#' @param widths A vector of integer widths for each field (assumes that
#'   columns are consecutive - that there is no overlap or gap between fields)
#' @export
hip_fwf_widths <- function(
  widths, col_names, col_types, trim_ws = TRUE, imp_dec = 0
) {
  if (!is_integerish(widths)) stop("fwf column widths must be integers.")
  if (!all(widths > 0)) stop("fwf column widths must be greater than 0.")

  pos <- cumsum(c(1L, widths))

  hip_fwf_positions(
    pos[-length(pos)], pos[-1] - 1L, col_names, col_types, trim_ws, imp_dec
  )
}



#' Create a record type information object
#'
#' Create a record type information object for hipread to use
#' when reading hierarchical files. A width of 0 indicates that
#' the file is rectangular (eg a standard fixed width file).
#'
#' @param start Start position of the record type variable
#' @param width The width of the record type variable
#' @param warn_on_missing Whether to warn when encountering a
#'   record type that is not specified
#'
#' @return A list, really only intended to be used internally by hipread
#' @export
hip_rt <- function(start, width, warn_on_missing = TRUE) {
  if (length(start) > 1) stop(paste0(
    "rectype start must be a single integer, but is of size ", length(start)
  ))

  if (length(width) > 1) stop(paste0(
    "rectype width must be a single integer, but is of size ", length(width)
  ))

  if (!is_integerish(start)) stop(paste0(
    "rectype start must be an integer but is ", start
  ))

  if (!is_integerish(width)) stop(paste0(
    "rectype width must be an integer but is ", width
  ))

  if (start < 1) stop(paste0(
    "rectype start must be greater than 1 but is ", start
  ))

  list(
    start = start - 1, width = width,
    verbose_warning = as.logical(warn_on_missing)
  )
}
