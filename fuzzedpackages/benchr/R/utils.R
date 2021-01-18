# Capture unevaluated dots
dots <- function(...) {
  eval(substitute(alist(...)))
}

# Capture unevaluated dots with names
named_dots <- function(...) {
  args <- dots(...)
  nm <- names(args)
  deparse2 <- function(e) {
    res <- deparse(e, 500L)
    res <- gsub("^\\s+", "", res)
    res <- gsub("([^{}])$", "\\1;", res)
    res <- paste(res, collapse = " ")
    res <- gsub(";( [{}]|$)", "\\1", res)
    res
  }
  deparsed <- vapply(args, deparse2, character(1L), USE.NAMES = FALSE)
  ifelse(is.null(nm), nm <- deparsed, nm[nm == ""] <- deparsed[nm == ""])
  names(args) <- nm
  args
}

# Check if a package is installed
is_installed <- function(pkg) {
  nzchar(system.file(package = pkg))
}

# Grouped summary stats for the benchmark object
summarise <- function(x, cols, fun, relative) {
  res <- lapply(split(x$time, x$expr), fun)
  res <- structure(
    cbind(names(res), as.data.frame(do.call(rbind, res))),
    class = c("summaryBenchmark", "data.frame"),
    names = c("expr", "n.eval", cols),
    row.names = c(NA, -length(res)),
    units = attr(x, "units")
  )
  if (nrow(res) > 1L && !is.null(relative)) {
    nm <- names(res)
    if (is.numeric(relative)) {
      relative <- nm[relative]
    }
    relative <- match.arg(relative, cols)
    res$relative <- signif(res[[relative]] / min(res[[relative]]), 3)
  }
  res$expr <- as.factor(res$expr)
  res
}
