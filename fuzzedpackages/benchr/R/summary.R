#' @title Summary method for the `benchmark` timings.
#'
#' @description
#' This method computes aggregated statistics (quantiles, means and sums)
#' for each expression.
#'
#' @param object An object of class `benchmark`.
#' @param relative Character. The name or index of the column whose values
#' will be used to compute relative timings.
#' @param \dots Not currently used.
#'
#' @return The method returns a `data.frame` with additional attributes,
#' which contains these columns:
#' \item{expr}{The deparsed expression as passed to
#' `benchmark` or the name of the argument if the expression was
#' passed as a named argument.}
#' \item{n.eval}{Number of successful measurements.}
#' \item{min}{Minimal timing measurement for this expression.}
#' \item{lw.qu}{First quartile of measurements for this expression.}
#' \item{mean}{Sample mean of measurements for this expression.}
#' \item{median}{Sample median of measurements for this expression.}
#' \item{up.qu}{Third quartile of measurements for this expression.}
#' \item{max}{Maximal timing measurement for this expression.}
#' \item{total}{Total (summed) measured time for this expression.}
#' \item{relative}{Relative difference across expressions compared to a minimal
#'  value in the column, specified by the argument `relative`.}
#'
#' Additional attributes:
#' \item{units}{Units for time intervals.}
#'
#' @importFrom stats quantile
#' @include utils.R
#' @method summary benchmark
#' @export
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@gmail.com}
#'
#' @seealso
#' [mean.benchmark()]
#'
#' @examples
#' timings <- benchmark(
#'   rchisq(100, 0), rchisq(100, 1), rchisq(100, 2),
#'   rchisq(100, 3), rchisq(100, 5),
#'   times = 1000L
#' )
#' summary(timings)
summary.benchmark <- function(object, relative = "median", ...) {
  cols <- c("min", "lw.qu", "median", "mean", "up.qu", "max", "total")
  fun <- function(x) {
    if (anyNA(x)) x <- x[!is.na(x)]
    qq <- quantile(x)
    s <- sum(x)
    n <- length(x)
    c(n, qq[1L:3L], s / n, qq[4L:5L], s)
  }
  summarise(object, cols, fun, relative)
}
