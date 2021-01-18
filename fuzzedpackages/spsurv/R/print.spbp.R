#' Bernstein Polynomial Based Regression Object Print
#'
#' @export
#' @param x an object of class spbp
#' @param digits number of digits to display
#' @param signif.stars see \code{\link{getOption}}
#' @param ... further arguments passed to or from other methods
#' @method print spbp
#' @return none


print.spbp <-
  function(x, digits = max(getOption('digits')-4, 3),
           signif.stars = getOption("show.signif.stars"), ...){

  savedig <- options(digits = digits)
  on.exit(options(savedig))

  if (!is.null(x$call)) {
    cat("\n---")
    cat("\n")
    cat("Call:\n")
    dput(x$call)
    cat("\n")
  }

  if(x$call$approach == "mle"){

    coef <- as.array(x$coefficients[1:x$q])
    var <- as.array(x$var[1:x$q, 1:x$q])

    ### Error handling ###
    # Null model
    if (is.null(x$coefficients)) return(x)

    coef2 <- coef[!(is.na(coef))] #non-missing coefs
    if(is.null(coef) | is.null(var )) stop("Input is not valid")

    se <- as.array(suppressWarnings(sqrt(diag(x$var)[1:x$q])))

    Coefmat  <- cbind(coef, exp(coef), se, coef/se,
                                  pchisq((coef/ se)^2, 1, lower.tail=FALSE))
    dimnames(Coefmat) <- list(names(coef), c("coef", "exp(coef)",
                                                         "se(coef)", "z", "Pr(>|z|)"))

    if(!is.null(x$coefficients)) {
      cat("\n")
      printCoefmat(Coefmat, digits = digits,
                   signif.stars = signif.stars, ...)
    }
    if(!is.null(x$loglik)) {
      cat("\n Loglik(model)= ", x$loglik[2])
      cat("      Loglik(no predictors)= ", x$loglik[1], "\n")
    }
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    cat("      Chisq= ", logtest," on ", x$q, " degrees of freedom ",
        pchisq(logtest, x$q, lower.tail=FALSE), "\n")

    if(!is.null(x$n)) {
        cat("n= ", x$n)
    }
  }
  else{
    summarise <- rstan::summary(x$stanfit, pars = "beta")$summary
    design <- as.matrix(model.matrix(x))
    p <- ncol(design)
    Coef <- cbind(matrix(summarise[, 1], nrow = p),
                  matrix(coda::HPDinterval(coda::mcmc(rstan::extract(x$stanfit, "beta")$beta)), nrow = p),
                  matrix(summarise[, -c(1, 5, 7, 9, 10)], nrow = p))

    rownames(Coef) <-  colnames(design)
    colnames(Coef) <- c("mean", "lowerHPD", "upperHPD", colnames(summarise)[-c(1, 5, 7, 9, 10)])
    print(Coef, digits = digits)

    cat("---\n")
    cat("\n WAIC Estimate= ", sprintf('%.3f', x$waic$estimates[3,1]))
    cat("      WAIC SE= ", sprintf('%.3f', x$waic$estimates[3,2], "\n"))
    cat("\n LOOIC Estimate= ", sprintf('%.3f', x$loo$estimates[3,1]))
    cat("      LOOIC SE= ", sprintf('%.3f', x$loo$estimates[3,2], "\n"))
    cat("\n---\n")
    cat("\n")
  }
}

