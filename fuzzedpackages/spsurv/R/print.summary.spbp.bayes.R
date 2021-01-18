#' Bernstein Polynomial Based Regression Object Summary Bayes
#'
#' @export
#' @param x a summary.spbp.bayes object
#' @param digits number of digits to display.
#' @param ... further arguments passed to or from other methods
#' @method print summary.spbp.bayes
#' @return none

print.summary.spbp.bayes <- ## summary printings
  function(x, digits = max(getOption('digits')-4, 3), ...){
    if (!is.null(x$call)) {
      cat("\n")
      cat("Call:\n")
      dput(x$call)
      cat("\n")
    }
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    cat("  n=", x$n)
    if (!is.null(x$nevent)) cat(", number of events=", x$nevent, "\n")
    else cat("\n")

    cat("\n")
    print(x$summary_chain)

    cat("---\n")
    print(x$summary_exp)

    cat("---\n")
    print(x$waic)
    print(x$loo)

    invisible()
  }
