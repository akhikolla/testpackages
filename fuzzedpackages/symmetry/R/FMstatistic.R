#' @template test_stat
#' @templateVar name The Cabilio–Masaro test statistic
#' @templateVar cite Cabilio1996
#' @templateVar formula to be added
#' @export
FM <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  FM_Cpp(X)
}
