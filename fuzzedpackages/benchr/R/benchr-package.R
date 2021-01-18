#' @title High Precise Measurement of R Expressions Execution Time
#'
#' @description
#' Package \pkg{benchr} provides an infrastructure (framework) for precise
#' measurement of R expressions execution time.
#'
#' @details
#' To measure execution time, \pkg{benchr} provides function
#' [benchmark()], as well as a number of additional methods for
#' analysis and representation of results.
#'
#' For precise time measurement we use a cross-platform monotone clock, provided
#' by C++11 standard in header file `chrono`. The timer accuracy depends on the
#' implementation by the compiler in use, the OS and the hardware. Most
#' commonly, the precision is one micro- or nanosecond. We provide the
#' opportunity to get the timer accuracy (time interval between two consecutive
#' timer ticks) via function [timer_precision()]. This accuracy is
#' also listed in the output of implicit or explicit `print` call.
#'
#' We estimate the timer overhead before the actual measurement by running
#' multiple (`2*10^5` by default) calls to an empty function. By doing that, we not
#' only estimate the overhead, but also produce a warm-up effect on the
#' processor, taking it out from idle state. After the actual measurement
#' results are adjusted by the timer overhead.
#'
#' Time intervals are measured in seconds and stored as `long double`,
#' which lets us capture a wide range of possible values from
#' `.Machine$double.xmin` to `.Machine$double.xmax`. This is quite
#' enough to operate both within very small (nanoseconds) and very big time
#' frames (e.g. the maximum time interval possible is
#' `.Machine$double.xmax / 3600` hours).
#'
#' It should be noted that the R session is not an isolated container with
#' strictly bounded resources, therefore the execution time can be influenced by
#' various factors, which may lead to outliers. In order to increase measurement
#' reliability, we repeat executions multiple times (100 repetitions for each
#' expression by default). This approach allows to collect enough data for
#' statistical analysis in time difference.
#'
#' We have also implemented several execution regimes in order to minimize the
#' probability of systematic errors in measurements. By default, a random order
#' of execution is being used. There is also a block order of execution, when
#' the first expression is repeated a fixed number of times, then the second and
#' so on. In such regime one can decrease the influence of allocators and
#' garbage collection, since the memory is allocated only at the beginning of
#' each block. The third option is to execute expressions in the order, provided
#' by the user.
#'
#' Note that we do not make any checks regarding return objects, i.e. one can
#' compare not only algorithms with the same result, but also the same algorithm
#' with different input parameters (e.g. input data sets of different size).
#'
#' We also do not check whether the expressions are `language` objects (see
#' [is.language()]) and do not coerce to that type.
#'
#' @section Package options:
# " Some of the available functionality is hidden from the user and not
#' accessible through function arguments. We allow to modify these parameters
#' via package options. We have tried to set optimal default values, which you
#' may consider changing in some cases. Here's a complete list of package options:
#' \describe{
#'   \item{\option{benchr.warmup}}{Number of iterations for timer overhead
#'   estimation (`2*10^5` by default).}
#'   \item{\option{benchr.print_details}}{Whether additional information
#'   on the measurement parameters should be displayed (`FALSE` by default).}
#'   \item{\option{benchr.use_ggplot}}{Whether \pkg{ggplot2} package should be
#'   used to produce plots, if the package is installed (`TRUE` by default).}
#' }
#'
#' @examples
#' # Benchmark expressions
#' res <- benchmark(
#'   rep(1:10, each = 10),
#'   rep.int(1:10, rep.int(10, 10))
#' )
#' # Aggregated statistics
#' mean(res)
#' summary(res)
#' # Plot results
#' boxplot(res)
#' @name benchr
#' @docType package
#'
#' @useDynLib benchr, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppProgress
#'
"_PACKAGE"
