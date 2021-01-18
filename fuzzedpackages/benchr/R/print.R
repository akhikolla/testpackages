#' @title Print method for the `benchmark` timings.
#'
#' @description
#' This is universal method of measurement results representation, which can be
#' called either implicitly or explicitly. The method uses `summary` method
#' to compute aggregated statistics for benchmarking results. `print` also
#' provides additional information about the timer precision and overhead, the
#' execution regime and the number of repeats.
#'
#' @param x An object of class `benchmark`, `summaryBenchmark` or
#' `meanBenchmark`.
#' @param units Character. The units to be used in printing the timings.
#' The available units are nanoseconds (`"ns"`), microseconds
#' (`"us"`), milliseconds (`"ms"`), seconds (`"s"`).
#' @param order Character. Order results according to this column of the output.
#' @param relative Character. The name or index of the column whose values will
#' be used to compute relative timings.
#' @param details Logical. Show additional detauls about measurement process.
#' @param ... Arguments passed on to [print.data.frame()].
#'
#' @return Apart from the table output produced by `summary`, the method
#' also prints additional information about the benchmarking process (with
#' `details = TRUE`):
#' \item{Timer precision}{Timer precision in seconds.}
#' \item{Timer error}{Timer error (overhead) in seconds.}
#' \item{Replications}{Number of repeats for each expression.}
#' \item{Expressions order}{Execution regime.}
#' \item{Garbage collection}{Whether garbage collection took place before each
#' execution.}
#'
#' @include summary.R
#' @method print benchmark
#' @export
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@gmail.com}
#'
#' @seealso
#' [summary.benchmark()],
#' [mean.benchmark()]
#'
#' @examples
#' a1 <- a2 <- a3 <- a4 <- numeric(0)
#' res <- benchmark(
#'   a1 <- c(a1, 1),
#'   a2 <- append(a2, 1),
#'   a3[length(a3) + 1] <- 1,
#'   a4[[length(a4) + 1]] <- 1
#' )
#' print(res)
print.benchmark <- function(x, units = "auto", order = "none",
                            relative = "median", details = FALSE, ...) {
  if (getOption("benchr.print_details", details)) {
    cat("Benchmark details:\n")
    print_details(x)
  }
  cat("Benchmark summary:\n")
  print(summary(x, relative), units, order, ...)
  invisible(x)
}

# Print a benchmark details info
#' @include units.R
print_details <- function(x) {
  fields <- c(
    "Timer precision",
    "Timer error",
    "Replications",
    "Expressions order",
    "Garbage collection"
  )
  values <- c(
    format_units(attr(x, "precision")),
    format_units(attr(x, "error")),
    attr(x, "times"),
    attr(x, "order"),
    attr(x, "gc")
  )
  cat(paste(format(fields, justify = "right"), ":", values), sep = "\n")
  invisible(x)
}

#' @method print summaryBenchmark
#' @export
#' @rdname print.benchmark
#' @include units.R
print.summaryBenchmark <- function(x, units = "auto", order = "none", ...) {
  nm <- names(x)
  cols <- nm[!nm %in% c("expr", "n.eval", "relative")]
  xx <- convert_units(x, units)
  if (order != "none") {
    xx <- xx[do.call(base::order, xx[order]), ]
  }
  xx[cols] <- signif(xx[cols], 3L)
  cat("Time units", ":", pretty_unit(attr(xx, "units")), "\n")
  print.data.frame(xx, ..., row.names = FALSE)
  invisible(x)
}
