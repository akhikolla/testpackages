# =============================== all_max_rcpp ============================== #

#' Sliding and disjoint block maxima
#'
#' Calculates the (sliding) maxima of all blocks of \code{b} contiguous values
#' and all sets of the maxima of disjoint blocks of \code{b} contiguous values
#' in the vector \code{x}.  This provides the first step of computations in
#' \code{\link{spm}}.
#'
#' @param x A numeric vector of raw observations.
#' @param b A numeric scalar.  The block size.
#' @param which_dj A character scalar.  Determines Which sets of disjoint
#'   maxima are calculated: \code{"all"}, all sets; \code{"first"}, only the
#'   set whose first block starts on the first observation in \code{x};
#'   \code{"last"}, only the set whose last block end on the last observation
#'   in \code{x}.
#' @param ... Further arguments to be passed to
#'   \code{\link[RcppRoll:RcppRoll-exports]{roll_max}}.
#' @details \strong{Sliding maxima.} The function
#'   \code{\link[RcppRoll:RcppRoll-exports]{roll_max}} in the \code{RcppRoll}
#'   package is used.
#'
#'   \strong{Disjoint maxima.}  If \code{n = length(x)} is an integer
#'   multiple of \code{b}, or if \code{which_dj = "first"} or
#'   \code{which_dj = "last"} then only one set of \code{n / b} disjoint
#'   block maxima are returned.
#'   Otherwise, \code{n - floor(n / b) * b + 1} sets of \code{floor(n / b)}
#'   disjoint block maxima are returned.  Set \code{i} are the disjoint maxima
#'   of \code{x[i:(i + floor(n / b) * b - 1)]}.  That is, all possible sets
#'   of contiguous disjoint maxima achieving the maxima length of
#'   \code{floor(n / b)} are calculated.
#'
#'   In both instances \code{na.rm = TRUE} is passed to \code{\link{max}} so
#'   that blocks containing missing values produce a non-missing result.
#'
#'   Also returned are the values in \code{x} that contribute to each set
#'   of block maxima.
#' @return A list containing
#'     \item{\code{ys} }{a numeric vector containing one set of sliding
#'       block maxima.}
#'     \item{\code{xs}}{a numeric vector containing the values that
#'       contribute to \code{ys}, that is, the whole input vector \code{x}.}
#'     \item{\code{yd} }{if \code{which_dj = "all"} a \code{floor(n / b)}
#'       by \code{n - floor(n / b) * b + 1} numeric matrix.  Each column
#'       contains a set of disjoint maxima.  Otherwise, a \code{floor(n / b)}
#'       by 1 numeric matrix containing one set of block maxima.}
#'     \item{\code{xd} }{if \code{which_dj = "all"} a
#'       \code{floor(n / b) * b} by \code{n - floor(n / b) * b + 1} numeric
#'       matrix.  Each column contains the values in \code{x} that contribute
#'       to the corresponding column in \code{yd}.  Otherwise, a
#'       \code{floor(n / b)} by 1 numeric matrix containing one the one set of
#'       the values in \code{x} that contribute to \code{yd}.}
#' @seealso \code{\link{spm}} for semiparametric estimation of the
#'   extremal index based on block maxima.
#' @examples
#' x <- 1:11
#' all_max_rcpp(x, 3)
#' all_max_rcpp(x, 3, which_dj = "first")
#' all_max_rcpp(x, 3, which_dj = "last")
#' @export
all_max_rcpp <- function(x, b = 1, which_dj = c("all", "first", "last"), ...){
  which_dj <- match.arg(which_dj)
  # First calculate the sliding block maxima.  All the disjoint maxima that
  # we need are contained in s_max, and we need the sliding maxima anyway
  ys <- RcppRoll::roll_max(x = x, n = b, ...)
  # The number of maxima of blocks of length b
  n <- length(x)
  n_max <- floor(n / b)
  # Set the possible first indices
  first_value <- switch(which_dj,
                        all = 1:(n - n_max * b + 1),
                        first = 1,
                        last = n - n_max * b + 1)
  # A function to return block maxima and contributing values starting from
  # the first value first
  get_maxima <- function(first) {
    s_ind <- seq.int(from = first, by = b, length.out = n_max)
    return(c(ys[s_ind], x[first:(first + n_max * b - 1)]))
  }
  temp <- vapply(first_value, FUN = get_maxima, numeric(n_max * (b + 1)))
  yd <- temp[1:n_max, , drop = FALSE]
  xd <- temp[-(1:n_max), , drop = FALSE]
  return(list(ys = ys, xs = x, yd = yd, xd = xd))
}
