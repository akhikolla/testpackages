# Print time unit with the full name
pretty_unit <- function(x) {
  stopifnot(is.character(x))
  stopifnot(length(x) == 1L)
  switch(x,
    ns = "nanoseconds",
    us = "microseconds",
    ms = "milliseconds",
    s = "seconds",
    m = "minutes",
    h = "hours",
    default = stop("Unknown unit ', x, '.")
  )
}

# Find the appropriate time unit for the min value
find_unit <- function(x) {
  stopifnot(is.numeric(x))
  if (anyNA(x)) {
    x <- x[!is.na(x)]
  }
  xx <- min(x[x > 0])
  if (xx < 1E-6) {
    units <- "ns"
  } else if (xx < 1E-3) {
    units <- "us"
  } else if (xx < 1E-0) {
    units <- "ms"
  } else if (xx < 60) {
    units <- "s"
  } else if (xx < 3600) {
    units <- "m"
  } else {
    units <- "h"
  }
  units
}

# Time units intervals
time_intervals <- c(
  ns = 1E-9,
  us = 1E-6,
  ms = 1E-3,
  s = 1E0,
  m = 60,
  h = 3600
)

# Convert time units
convert_units <- function(x, units = "auto") {
  stopifnot(is.data.frame(x))
  from <- attr(x, "units")
  if (is.null(from)) from <- "s"
  units <- match.arg(units, c("auto", names(time_intervals)))
  nm <- names(x)
  cols <- nm[!nm %in% c("expr", "n.eval", "relative")]
  if (units == "auto") {
    units <- find_unit(unlist(x[cols]) * time_intervals[from])
  }
  num <- time_intervals[from]
  denum <- time_intervals[units]
  x[cols] <- x[cols] * num / denum
  attr(x, "units") <- units
  x
}

format_units <- function(x, units = "auto", ...) {
  stopifnot(is.numeric(x))
  from <- attr(x, "units")
  if (is.null(from)) from <- "s"
  units <- match.arg(units, c("auto", names(time_intervals)))
  if (units == "auto") {
    units <- find_unit(x * time_intervals[from])
  }
  num <- time_intervals[from]
  denum <- time_intervals[units]
  x <- x * num / denum
  paste(x, pretty_unit(units))
}
