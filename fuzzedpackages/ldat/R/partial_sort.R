
#' Partial sort an lvec
#'
#' @param x an object of type \code{\link{lvec}}
#' @param pivots a numeric vector with indices at which the vector will be 
#'   sorted. See details for more information.
#' @param clone clone the vector first before sorting; or sort (and therefore
#'   modify) the input vector directly.
#'
#' @details
#' After partial sorting the vector values at the pivots are the same as the
#' vector values of a completely sorted vector. Furthermore, for each pivot
#' \code{i} all elements \code{x[j]; j < i} are smaller or equal to than 
#' \code{x[i]} and all elements \code{x[j]; j > i} are larger than or equal to
#' \code{x[i]}.
#' 
#' The speed of this operation should be O(n, k) with n the size of the 
#' \code{lvec}  and k the number of pivots. 
#'
#' @examples
#' x <- as_lvec(rnorm(100))
#' y <- partial_sort(x, c(10, 50, 90))
#' x_sorted <- sort(x)
#' stopifnot(all(y[c(10, 50, 90)] == x_sorted[c(10, 50, 90)]))
#' stopifnot(max(y[1:9]) <= min(y[11:100]))
#' stopifnot(max(y[1:49]) <= min(y[51:100]))
#' stopifnot(max(y[1:89]) <= min(y[91:100]))
#' 
#' @rdname partial_sort
#' @useDynLib ldat
#' @importFrom Rcpp sourceCpp
#' @export
partial_sort <- function(x, pivots, clone = TRUE) {
  if (clone) x = clone(x)
  pivots <- round(as.numeric(pivots))
  if (min(pivots) < 1) stop("Pivots smaller than 0 found.")
  if (max(pivots) > length(x)) stop("Pivots larger than vector length found.")
  pivots <- sort(unique(pivots))
  requireNamespace("lvec")
  #.Call("partial_sort", x, pivots, PACKAGE = "ldat")
  partial_sort_cpp(x, pivots)
  x
}

#' @rdname partial_sort
#' @useDynLib ldat
#' @export
partial_order <- function(x, pivots) {
  pivots <- round(as.numeric(pivots))
  if (min(pivots) < 1) stop("Pivots smaller than 0 found.")
  if (max(pivots) > length(x)) stop("Pivots larger than vector length found.")
  pivots <- sort(unique(pivots))
  requireNamespace("lvec")
  #structure(.Call("partial_order", x, pivots, PACKAGE = "ldat"), 
  structure(partial_order_cpp(x, pivots), class = "lvec")
}


#' Calculate the quantiles of an lvec
#'
#' @param x an object of type \code{\link{lvec}}.
#' @param probs a numeric vector with probabilities ([0,1]).
#' @param names add names to the result vector.
#' @param na.rm remove missing values before calculating the quantiles
#' @param true_probs add an attribute with the probabilities at the chosen 
#'   pivots.
#' @param ... ignored.
#'
#' @details
#' This function uses a more simple method than that used by the regular 
#' \code{\link{quantile}} method. It sorts the vector (using 
#' \code{\link{partial_sort}} for speed) and selects elements from \code{x} 
#' that correspond to the given probabilities. For example, when \code{x} has 
#' length of 11 and \code{prob} equal to 0.5, it selects the 6th element from 
#' the (partially) sorted \code{x}. For large enough vectors this is a 
#' reasonable approach. 
#'
#' @importFrom stats quantile
#' @export
quantile.lvec <- function(x, probs = seq(0, 1, 0.25), names = TRUE, 
    na.rm = TRUE, true_probs = FALSE, ...) {
  # Based on function with same name from stats
  format_perc <- function(x, digits = max(2L, getOption("digits"))) {
    if (length(x)) {
      x <- x * 100
      paste0(format(x, trim = TRUE, digits = digits), "%")
    } else character(0)
  }

  pivots <- round((length(x) - 1) * probs + 1)
  y <- if (na.rm) lget(x, !is.na(x)) else clone(x)
  y <- partial_sort(y, pivots = pivots, clone = FALSE)
  result <- as_rvec(lget(y, pivots))
  if (true_probs) {
    attr(result, "probs") <- (pivots - 1) / (length(x) - 1)
  }
  if (names) {
    names(result) <- format_perc(probs)
  }
  result

}

#' Calculate the median of an lvec
#' 
#' @param x an object of type \code{\link{lvec}}.
#' @param na.rm remove missing values before calculating the quantiles
#' @param ... ignored.
#'
#' @seealso
#' For more details see \code{\link{quantile.lvec}}.
#'
#' @importFrom stats median
#' @export
median.lvec <- function(x, na.rm = TRUE, ...) {
  res <- quantile(x, probs = 0.5, na.rm = na.rm, names = FALSE)
  as.numeric(res)
}

