
#' Print an object of class \code{srank}
#' @param x An object of class \code{srank}
#' @param digits the number of significant digits to use when printing
#' @param ... further arguments passed to or from other methods
#' @export
print.srank <- function(x, digits = 3, ...) {

    cat("\n")
    if (is.null(digits))
        digits <- getOption("digits")
    print.data.frame(x$summary, digits = digits)

    invisible(x)
}

#' Print an object of class \code{strat}
#' @param x An object of class \code{strat}
#' @param digits the number of significant digits to use when printing
#' @param ... further arguments passed to or from other methods
#' @export
print.strat <- function(x, digits = 3, ...) {

    cat("\n")
    if (is.null(digits))
        digits <- getOption("digits")
    cat("overall stratification:\n\n")
    print.default(x$overall, digits = digits)

    if(!is.null(x$decomposition)){
      cat("\n")
      cat(paste("decomposition by ", colnames(x$within_group)[1],
          ":\n\n", sep = ""))
      print.default(x$decomposition, digits = digits)
    }

    invisible(x)
}
