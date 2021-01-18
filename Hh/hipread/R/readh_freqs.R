#' Calculate frequencies from fixed width file without loading into memory
#'
#' Calculate the frequency of values in all variables in a fixed width file.
#' Does so without holding the whole data in memory or creating a full
#' R data.frame and calling R code on interim pieces. (Probably only
#' useful inside IPUMS HQ).
#'
#' @inheritParams hipread_long
#'
#' @return A list of frequencies
#' @export
#' @keywords internal
hipread_freqs <- function(
  file, var_info, rt_info = hip_rt(1, 0),
  compression = NULL, progress = show_progress()
) {
  file <- check_file(file)
  isgzipped <- is_gzip_compression(compression, file)
  var_info <- add_level_to_rect(var_info)
  var_names <- get_var_names(var_info)
  var_pos_info <- get_var_pos(var_info, var_names)
  var_types <- get_var_types(var_info, var_names)

  read_freqs(
    file, var_names, rt_info, var_pos_info, isgzipped, progress
  )
}
