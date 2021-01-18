#' Print the summary.ypbp output
#'
#' @export
#' @param x an object of the class summary.ypbp.
#' @param ... further arguments passed to or from other methods.
#' @return a summary of the fitted model.
print.summary.ypbp <- function(x, ...){
  if(x$approach=="mle"){
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Short-term coefficients:\n")
    stats::printCoefmat(x$coefficients1, P.value=TRUE, has.Pvalue=TRUE)
    cat("\n")
    # cat("----------------------- \n")
    # cat("\n")
    cat("Long-term coefficients:\n")
    stats::printCoefmat(x$coefficients2, P.value=TRUE, has.Pvalue=TRUE)
    if(is.null(x$coefficients3)){
      cat("\n")
      cat("--- \n")
      cat("loglik =", x$loglik, " ", "AIC =", x$AIC,"\n")
    }else{
      cat("\n")
      cat("Proportional hazards coefficients:\n")
      stats::printCoefmat(x$coefficients3, P.value=TRUE, has.Pvalue=TRUE)
      cat("\n")
      cat("\n")
      cat("--- \n")
      cat("loglik =", x$loglik, " ", "AIC =", x$AIC,"\n")
    }

  }else{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Short-term coefficients:\n")
    print(x$coefficients1)
    cat("\n")
    cat("Long-term coefficients:\n")
    print(x$coefficients2)
    cat("\n")
    if(!is.null(x$coefficients3)){
      cat("Proportional hazards coefficients:\n")
      print(x$coefficients3)
    }
    cat("--- \n")
    cat("Inference for Stan model: ", x$model_name, '.\n', sep = '')
    cat(x$chains, " chains, each with iter=", x$iter,
        "; warmup=", x$warmup, "; thin=", x$thin, "; \n",
        "post-warmup draws per chain=", x$n_kept[1], ", ",
        "total post-warmup draws=", sum(x$n_kept), ".\n\n", sep = '')

  }

}


#---------------------------------------------

#' Summary for the ypbp model
#'
#' @aliases summary.ypbp
#' @export
#' @param object an objecto of the class 'ypbp'.
#' @param ... further arguments passed to or from other methods.
#'
summary.ypbp <- function(object, ...){
  p <- object$p
  q <- object$q
  if(object$approach=="mle"){
    n <- object$n
    k <- p+2*q+object$degree

    loglik <- object$fit$value
    AIC <- -2*loglik + 2*k
    BIC <- -2*loglik + k*log(n)

    labels <- object$labels
    coefficients <- object$fit$par[1:(2*q+p)]
    vcov <- MASS::ginv(-object$fit$hessian)[1:(2*q+p),1:(2*q+p)]

    se <- sqrt(diag(vcov))
    zval <- coefficients / se
    TAB <- cbind(Estimate = coefficients,
                 StdErr = se,
                 z.value = zval,
                 p.value = 2*stats::pnorm(-abs(zval)))

    if(q==1){
      TAB1 <- t(as.matrix(TAB[1:q,]))
      TAB2 <- t(as.matrix(TAB[(q+1):(2*q),]))
    }else{
      TAB1 <- TAB[1:q,]
      TAB2 <- TAB[(q+1):(2*q),]
    }
    rownames(TAB1) <- labels
    rownames(TAB2) <- labels

    if(p==1){
      TAB3 <- t(as.matrix(TAB[(2*q+1):(2*q+p),]))
      rownames(TAB3) <- object$labels.ph
      res <- list(call=object$call,
                  coefficients1=TAB1, coefficients2=TAB2, coefficients3=TAB3,
                  loglik=loglik, AIC=AIC, approach=object$approach)
    }
    if(p>1){
      TAB3 <- TAB[(2*q+1):(2*q+p),]
      rownames(TAB3) <- object$labels.ph
      res <- list(call=object$call,
                  coefficients1=TAB1, coefficients2=TAB2, coefficients3=TAB3,
                  loglik=loglik, AIC=AIC, approach=object$approach)
    }
    if(p==0){
      res <- list(call=object$call,
                  coefficients1=TAB1, coefficients2=TAB2,
                  loglik=loglik, AIC=AIC, approach=object$approach)
    }

   # Bayesiam output:
  }else{
    labels <- object$labels
    s1 <- rstan::summary(object$fit, pars=c("psi"))
    s2 <- rstan::summary(object$fit, pars=c("phi"))
    TAB1 <- round(s1$summary, digits = 3)
    TAB2 <- round(s2$summary, digits = 3)
    rownames(TAB1) <- labels
    rownames(TAB2) <- labels


    n_kept <- object$fit@sim$n_save - object$fit@sim$warmup2

    if(p==0){
      res <- list(call=object$call,
                  coefficients1=TAB1, coefficients2=TAB2,
                  n_kept=n_kept, model_name=object$fit@model_name,
                  chains=object$fit@sim$chains, warmup=object$fit@sim$warmup,
                  thin=object$fit@sim$thin, iter=object$fit@sim$iter, approach=object$approach)
    }else{
      if(p>0){
        s3 <- rstan::summary(object$fit, pars=c("beta"))
        TAB3 <- round(s3$summary, digits = 3)
        rownames(TAB3) <- object$labels.ph
        res <- list(call=object$call,
                    coefficients1=TAB1, coefficients2=TAB2, coefficients3=TAB3,
                    n_kept=n_kept, model_name=object$fit@model_name,
                    chains=object$fit@sim$chains, warmup=object$fit@sim$warmup,
                    thin=object$fit@sim$thin, iter=object$fit@sim$iter, approach=object$approach)
      }
    }

  }


  class(res) <- "summary.ypbp"

  return(res)
}


