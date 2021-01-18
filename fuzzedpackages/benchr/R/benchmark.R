#' @title High Precise Measurement of R Expressions Execution Time
#'
#' @description
#' `benchmark` serves as a more accurate replacement of the often
#' seen `system.time(replicate(1000, expr))` expression. It tries hard to
#' accurately measure only the time it takes to evaluate `expr`.
#' To achieve this, the sub-millisecond (supposedly nanosecond) accurate
#' timing functions most modern operating systems provide are used.
#' Additionally all evaluations of the expressions are done in C++ code to
#' minimize any measurment error.
#'
#' @param ... Captures any number of unevaluated expressions passed to
#' benchmark as named or unnamed arguments.
#' @param times Integer. Number of times to evaluate each expression.
#' @param order Character. The order in which the expressions are evaluated.
#' @param envir The environment in which the expressions will be evaluated.
#' @param gcFirst Logical. Should a garbage collection be performed immediately
#' before the timing?
#' @param gcDuring Logical. Should a garbage collection be performed immediately
#' before each iteration, as produced by `times`? (very slow)
#' @param progress Logical. Show progress bar during expressions evaluation.
#'
#' @details
#' Before evaluating each expression `times` times, the overhead of
#' calling the timing functions and the C++ function call overhead are
#' estimated. This estimated overhead is subtracted from each measured
#' evaluation time. Should the resulting timing be negative, a warning
#' is thrown and the respective value is replaced by `0`. If the timing
#' is zero, a warning is also raised. Should all evaluations result in one of
#' the two error conditions described above, an error is raised.
#'
#' @section The order in which the expressions are evaluated:
#' \describe{
#'  \item{\dQuote{random}}{(the default) randomizes the execution order}
#'  \item{\dQuote{inorder}}{executes each expression in order}
#'  \item{\dQuote{block}}{executes all repetitions of each expression
#'     as one block.}
#' }
#'
#' @return Object of class `benchmark`, which is a `data.frame` with
#' a number of additional attributes and contains the following columns:
#' \item{expr}{The deparsed expression as passed to
#' `benchmark` or the name of the argument if the expression was
#' passed as a named argument.}
#' \item{time}{The measured execution time of the expression in seconds.
#' The order of the observations in the data frame is the order in which they
#' were executed.}
#' An object of class benchmark also contains the following attributes:
#' \item{precision}{Timer precision in seconds.}
#' \item{error}{Timer error (overhead) in seconds.}
#' \item{units}{Units for time intervals (by default, "s" -- seconds).}
#' \item{times}{Number of repeats for each measurement.}
#' \item{order}{Execution regime.}
#' \item{gc}{Whether garbage collection took place before each execution.}
#'
#' @include RcppExports.R
#' @include utils.R
#' @include order.R
#' @export
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@gmail.com}
#'
#' @seealso
#' [summary.benchmark()],
#' [mean.benchmark()],
#' [print.benchmark()],
#' [plot.benchmark()],
#' [boxplot.benchmark()]
#'
#' @examples
#' ## Measure the time it takes to dispatch a simple function call
#' ## compared to simply evaluating the constant NULL
#' f <- function() NULL
#' res <- benchmark(NULL, f(), times = 1000L)
#'
#' ## Print results:
#' print(res)
#'
#' ## Plot results
#' boxplot(res)
benchmark <- function(..., times = 100L, order = c("random", "inorder", "block"),
                      envir = parent.frame(), progress = TRUE,
                      gcFirst = TRUE, gcDuring = FALSE) {
  if (length(list(...)) == 0) {
    stop("No expressions to benchmark.")
  }
  if (times <= 0L) {
    stop("Argument 'times' should be positive integer.")
  }
  order <- match.arg(order)
  exprs <- named_dots(...)
  exprs_order <- make_order(exprs, times, order)
  warmup <- getOption("benchr.warmup", 2e5L)
  if (gcFirst) gc(FALSE)
  error <- timer_error(warmup)
  timings <- do_benchmark(exprs, parent.frame(), exprs_order, gcDuring, progress)
  timings <- timings - error
  timings[timings < 0] <- NA_real_
  if (anyNA(timings)) {
    nas <- sum(is.na(timings))
    if (nas == length(timings)) {
      stop("All timed evaluations were either smaller than the estimated \
                 timer error or zero. The most likely cause is a low resolution clock.")
    } else {
      warning(sprintf("Could not measure a positive execution time for \
                            %i evaluations.", nas))
    }
  }
  structure(
    list(expr = exprs_order, time = timings),
    class = c("benchmark", "data.frame"),
    row.names = c(NA_integer_, -length(timings)),
    units = "s",
    error = error,
    precision = timer_precision(),
    order = order,
    gc = gcDuring,
    times = times
  )
}
