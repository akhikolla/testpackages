#' Read a hierarchical fixed width data file, in yields
#'
#' Enhances `hipread_long()` or `hipread_list()` to allow you to read
#' hierarchical data in pieces (called 'yields') and allow your code to
#' have full control between reading pieces, allowing for more freedom
#' than the 'callback' method introduced in the chunk functions (like
#' [`hipread_long_chunked()`]).
#'
#' These functions return a HipYield R6 object which have the following
#' methods:
#' - `yield(n = 10000)` A function to read the next 'yield' from the data,
#'   returns a `tbl_df` (or list of `tbl_df` for `hipread_list_yield()`)
#'   with up to n rows (it will return NULL if no rows are left, or all
#'   available ones if less than n are available).
#' - `reset()` A function to reset the data so that the next yield will
#'   read data from the start.
#' - `is_done()` A function that returns whether the file has been completely
#'   read yet or not.
#' - `cur_pos` A property that contains the next row number that will be
#'    read (1-indexed).
#'
#' @inheritParams hipread_long
#' @return A HipYield R6 object (See 'Details' for more information)
#' @export
#' @examples
#' library(hipread)
#' data <- hipread_long_yield(
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
#' # Read the first 4 rows
#' data$yield(4)
#'
#' # Read the next 2 rows
#' data$yield(2)
#'
#' # Reset and then read the first 4 rows again
#' data$reset()
#' data$yield(4)
hipread_long_yield <- function(
  file, var_info, rt_info = hip_rt(1, 0), compression = NULL,
  skip = 0, encoding = "UTF-8"
) {
  HipLongYield$new(
    file,
    var_info,
    rt_info,
    compression,
    skip,
    encoding
  )
}

#' @export
#' @rdname hipread_long_yield
hipread_list_yield <- function(
  file, var_info, rt_info = hip_rt(1, 0), compression = NULL,
  skip = 0, encoding = "UTF-8"
) {
  HipListYield$new(
    file,
    var_info,
    rt_info,
    compression,
    skip,
    encoding
  )
}
