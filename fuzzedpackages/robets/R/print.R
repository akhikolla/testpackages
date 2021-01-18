#' Print robets model
#'
#' @param x An object of class \code{robets}.
#' @param ... Other undocumented arguments.
#' 
#' @examples
#' model <- robets(nottem)
#' print(model)
#' @export
print.robets <- function(x, ...)
{
  cat(paste(x$method, "\n\n"))
  cat(paste("Call:\n", deparse(x$call), "\n\n"))
  ncoef <- length(x$initstate)
  if(!is.null(x$lambda))
    cat("  Box-Cox transformation: lambda=",round(x$lambda,4), "\n\n")
  
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", round(x$par["alpha"], 4), "\n"))
  if(x$components[2]!="N")
    cat(paste("    beta  =", round(x$par["beta"], 4), "\n"))
  if(x$components[3]!="N")
    cat(paste("    gamma =", round(x$par["gamma"], 4), "\n"))
  if(x$components[4]!="FALSE")
    cat(paste("    phi   =", round(x$par["phi"], 4), "\n"))
  
  cat("\n  Initial states:\n")
  cat(paste("    sigma =", round(x$initstate[1], 4), "\n"))
  cat(paste("    l =", round(x$initstate[2], 4), "\n"))
  if (x$components[2]!="N")
    cat(paste("    b =", round(x$initstate[3], 4), "\n"))
  else
  {
    x$initstate <- c(x$initstate[1:2], NA, x$initstate[3:ncoef])
    ncoef <- ncoef+1
  }
  if (x$components[3]!="N")
  {
    cat("    s=")
    if (ncoef <= 9)
      cat(round(x$initstate[4:ncoef], 4))
    else
    {
      cat(round(x$initstate[4:9], 4))
      cat("\n           ")
      cat(round(x$initstate[10:ncoef], 4))
    }
    cat("\n")
  }
  
  cat("\n  sigma:  ")
  cat(round(sqrt(x$sigma2),4))
  if(!is.null(x$aic))
  {
    stats <- c(x$robaic,x$robaicc,x$robbic)
    names(stats) <- c("robAIC","robAICc","robBIC")
    cat("\n\n")
    print(stats)
  }
  #    cat("\n  AIC:    ")
  #    cat(round(x$aic,4))
  #    cat("\n  AICc:   ")
  #    cat(round(x$aicc,4))
  #    cat("\n  BIC:    ")
  #    cat(round(x$bic,4))
}