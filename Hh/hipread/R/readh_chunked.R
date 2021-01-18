#' Read a hierarchical fixed width data file, in chunks
#'
#' Analogous to [readr::read_fwf()], but with chunks, and allowing for
#' hierarchical fixed width data files (where the data file has rows of
#' different record types, each with their own variables and column
#' specifications). `hipread_long_chunked()` reads hierarchical data into "long"
#' format, meaning that there is one row per observation, and variables
#' that don't apply to the current observation receive missing values.
#' Alternatively, `hipread_list_chunked()` reads hierarchical data into "list"
#' format, which returns a list that has one data.frame per record type.
#'
#' @inheritParams hipread_long
#' @param callback A [`callback`] function, allowing you to perform a
#'   function on each chunk.
#' @param chunk_size The size of the chunks that will be read as a
#'   single unit (defaults to 10000)
#'
#' @return Depends on the type of [`callback`] function you use
#' @export
#'
#' @examples
#' # Read in a data, filtering out hhnum == "002"
#' data <- hipread_long_chunked(
#'   hipread_example("test-basic.dat"),
#'   HipDataFrameCallback$new(function(x, pos) x[x$hhnum != 2, ]),
#'   4,
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
hipread_long_chunked <- function(
  file, callback, chunk_size, var_info, rt_info = hip_rt(1, 0),
  compression = NULL, skip = 0, encoding = "UTF-8",
  progress = show_progress()
) {
  file <- check_file(file)
  isgzipped <- is_gzip_compression(compression, file)
  var_info <- add_level_to_rect(var_info)
  var_names <- get_var_names(var_info)
  var_pos_info <- get_var_pos(var_info, var_names)
  var_types <- get_var_types(var_info, var_names)
  var_opts <- get_var_opts(var_info, var_names)
  skip <- check_skip(skip)
  if (is.null(encoding)) encoding <- "UTF-8"

  callback <- as_chunk_callback(callback)
  on.exit(callback$finally(), add = TRUE)

  read_chunked_long(
    file, callback, chunk_size, var_names, var_types, rt_info,
    var_pos_info, var_opts, skip, isgzipped, encoding, progress
  )

  return(callback$result())
}


#' @rdname hipread_long_chunked
#' @export
hipread_list_chunked <- function(
  file, callback, chunk_size, var_info, rt_info = hip_rt(1, 0),
  compression = NULL, skip = 0, encoding = "UTF-8",
  progress = show_progress()
) {
  file <- check_file(file)
  isgzipped <- is_gzip_compression(compression, file)
  var_info <- add_level_to_rect(var_info)
  var_names <- get_vinfo_col_as_list(var_info, "col_names")
  var_pos_info <- get_var_pos(var_info)
  var_types <- get_vinfo_col_as_list(var_info, "col_types")
  var_opts <- get_var_opts_list(var_info)
  skip <- check_skip(skip)
  if (is.null(encoding)) encoding <- "UTF-8"

  callback <- as_chunk_callback(callback)
  on.exit(callback$finally(), add = TRUE)

  read_chunked_list(
    file, callback, chunk_size, var_names, var_types, rt_info,
    var_pos_info, var_opts, skip, isgzipped, encoding, progress
  )

  return(callback$result())
}

# Copied from readr
as_chunk_callback <- function(x) UseMethod("as_chunk_callback")
as_chunk_callback.function <- function(x) {
  HipSideEffectChunkCallback$new(x)
}
as_chunk_callback.R6ClassGenerator <- function(x) {
  as_chunk_callback(x$new())
}
as_chunk_callback.ChunkCallback <- function(x) {
  x
}
as_chunk_callback.HipChunkCallback <- function(x) {
  x
}
