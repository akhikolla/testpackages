#' @export
summary.migpd <-
function(object, ...){
  d <- dim(object$data)[2]
  conv <- sapply(object$models, function(z) z$convergence)

  if (sum(conv) != 0)
    conv <- (1:d)[conv != 0]
  else conv <- NULL

  co <- coef(object)

  res <- list(d=d,conv=conv,penalty=object$penalty,co=co)
  oldClass(res) <- "summary.migpd"
  res
}

#` @export 
print.summary.migpd <- function(object, ...){
    cat("\nA collection of", object$d, "generalized Pareto models.\n")
        
    if (is.null(object$conv)) cat("All models converged.\n")
    else cat("The following model(s) did not converge:", paste(object$conv, collapse=","), "\n")
        
    cat("Penalty to the likelihood:", object$penalty)
    cat("\n\nSummary of models:\n")
    print(object$co, ...)
    cat("\n")
    invisible(object)
}