#' @template test_stat
#' @templateVar name The \eqn{\sqrt{b_1}}  test statistic
#' @templateVar cite Cabilio1996
#' @templateVar formula to be added
#' @export
B1 <- function(X) {
  if (!is.numeric((X)) && !is.logical((X))) {
    warning("Argument is not numeric or logical: returning NA")
    return(NA)
  }
  B1_Cpp(X)
}
