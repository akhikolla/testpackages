
#---------------------------------------------

#' Print the summary.bellreg output
#'
#' @export
#' @param x an object of the class summary.bellreg.
#' @param ... further arguments passed to or from other methods.
#' @return a summary of the fitted model.
#'
print.summary.bellreg <- function(x, ...){
  if(x$approach == "mle"){
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Coefficients:\n")
    stats::printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    cat("\n")
    # cat("----------------------- \n")
    # cat("\n")
    cat("logLik =", x$logLik, " ", "AIC =", x$AIC,"\n")
  }else{
    cat("\n")
    print(x$call)
    cat("\n")
    print(x$coefficients)
    cat("\n")
    cat("Inference for Stan model: ", x$model_name, '.\n', sep = '')
    cat(x$chains, " chains, each with iter=", x$iter,
        "; warmup=", x$warmup, "; thin=", x$thin, "; \n",
        "post-warmup draws per chain=", x$n_kept[1], ", ",
        "total post-warmup draws=", sum(x$n_kept), ".\n\n", sep = '')
  }

}


#---------------------------------------------

#' Summary for the bellreg model
#'
#' @aliases summary.bellreg
#' @export
#' @param object an objecto of the class 'bellreg'.
#' @param ... further arguments passed to or from other methods.
#'
#'
summary.bellreg <- function(object, ...){

  if(object$approach == "mle"){
    p <- object$p
    labels <- object$labels
    coefficients <- object$fit$par
    V <- vcov(object)

    se <- sqrt(diag(V))
    zval <- coefficients / se
    TAB <- cbind(Estimate = coefficients,
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*stats::pnorm(-abs(zval)))

    if(p==1)
    {
      TAB <- t(as.matrix(TAB[1:p,]))
    }else{
      TAB <- TAB[1:p,]
    }

    rownames(TAB) <- labels
    res <- list(call=object$call,
                coefficients=TAB,
                logLik = object$logLik,
                AIC = object$AIC)
  }else{
    labels <- object$labels
    s <- rstan::summary(object$fit, pars=c("beta"))
    TAB <- round(s$summary, digits = 3)
    rownames(TAB) <- labels
    n_kept <- object$fit@sim$n_save - object$fit@sim$warmup2

    res <- list(call=object$call, coefficients=TAB,
                n_kept=n_kept, model_name=object$fit@model_name,
                chains=object$fit@sim$chains, warmup=object$fit@sim$warmup,
                thin=object$fit@sim$thin, iter=object$fit@sim$iter, approach=object$approach)

  }
  res$approach <- object$approach
  class(res) <- "summary.bellreg"
  return(res)
}

