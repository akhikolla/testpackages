#' @title Bind rows of multiple data frames with zero fill
#'
#' @description A version of bind_rows out of dplyr that
#'     fills non-common columns with zeroes instead of NA.
#'     Gives an error if any of the input data contains NA already.
#'
#' @param ... Data frames to combine, passed into bind_rows (see dplyr documentation)
#'
#' @template start-example
#' @examples
#' library(tibble)
#'
#' dose_info_A <- tibble(
#'    group_id = "hist_A",
#'    drug_A = 1
#'  )
#'
#'  dose_info_B <- tibble(
#'    group_id = "hist_B",
#'    drug_B = 100 * (1:2)
#'  )
#'
#' bind_rows_0(dose_info_A, dose_info_B)
#'
#' @template stop-example
#' @export
bind_rows_0 <- function(...) {
  dot.args <- list(...)
  full_args <- try(lapply(dot.args, na.fail))
  if (inherits(full_args, "try-error")) {
    stop("One or more data frame(s) passed into bind_rows_0 contain NA.")
    full_args <- dot.args
  }
  temp <- do.call("bind_rows", full_args)

  temp[is.na(temp)] <- 0
  temp
}
