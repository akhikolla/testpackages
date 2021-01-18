#' @title Mean method for the `benchmark` timings.
#'
#' @description
#' This method computes aggregated statistics (sample mean and confidence
#' intervals) for each expression.
#'
#' @param x An object of class `benchmark`.
#' @param trim Numeric. The fraction (0 to 0.5) of observations to be trimmed
#' before the mean is computed.
#' @param conf.level Numeric. Confidence level of the interval.
#' @param relative Character. The name or index of the column whose values
#' will be used to compute relative timings.
#' @param \dots Not currently used.
#'
#' @return The method returns a `data.frame` with additional attributes,
#' which contains these columns:
#' \item{expr}{The deparsed expression as passed to
#' `benchmark` or the name of the argument if the expression was
#' passed as a named argument.}
#' \item{mean}{Sample mean for timing results.}
#' \item{trimmed}{Trimmed sample mean for timing results (a fraction of
#' observations to be trimmed is defined by the argument `trim`).}
#' \item{lw.ci}{Lower boundary for the confidence level (confidence level is
#' specified by the argument `conf.level`).}
#' \item{up.ci}{Upper boundary for the confidence level (confidence level is
#' specified by the argument `conf.level`).}
#' \item{relative}{Relative difference across expressions compared to a minimal
#' value in the column, specified by the argument `relative`.}
#'
#' Additional attributes:
#' \item{units}{Units for time intervals.}
#' \item{conf.level}{Confidence level.}
#' \item{trim}{Fraction of observations that was trimmed before the trimmed
#' mean was computed.}
#'
#' @importFrom stats wilcox.test
#' @include utils.R
#' @method mean benchmark
#' @export
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@gmail.com}
#'
#' @seealso
#' [summary.benchmark()]
#'
#' @examples
#' timings <- benchmark(
#'   rchisq(100, 0), rchisq(100, 1), rchisq(100, 2),
#'   rchisq(100, 3), rchisq(100, 5),
#'   times = 1000L
#' )
#' mean(timings)
mean.benchmark <- function(x, trim = 0.05, conf.level = 0.95, relative = "mean", ...) {
  if (attr(x, "times") >= 5L) {
    cols <- c("mean", "trimmed", "lw.ci", "up.ci")
    fun <- function(x) {
      if (anyNA(x)) x <- x[!is.na(x)]
      s <- sum(x)
      n <- length(x)
      c(
        n, s / n, mean.default(x, trim = trim),
        wilcox.test(x, conf.int = TRUE, conf.level = conf.level)$conf.int
      )
    }
  } else {
    cols <- c("mean", "trimmed")
    fun <- function(x) {
      if (anyNA(x)) x <- x[!is.na(x)]
      s <- sum(x)
      n <- length(x)
      c(n, s / n, mean.default(x, trim = trim))
    }
  }
  summarise(x, cols, fun, relative)
}
