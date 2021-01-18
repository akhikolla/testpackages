

# Computes the Bernstein polynomial's bases. (note: for computation stability, b is not divided by tau here)
bp <- function(time, degree, tau) {
  n <- length(time)
  y <- time/tau
  b <- matrix(nrow=n, ncol=degree)
  B <- matrix(nrow=n, ncol=degree)
  for(k in 1:degree)
  {
    b[,k] <- stats::dbeta(y, k, degree - k + 1)
    B[,k] <- stats::pbeta(y, k, degree - k + 1)
  }
  return(list(b=b, B=B))
}

#---------------------------------------------
#' Variance-covariance matrix for a ypbp model
#'
#' @aliases vcov.ypbp
#' @description This function extracts and returns the variance-covariance matrix associated with the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class ypbp.
#' @param ... further arguments passed to or from other methods.
#' @return  the variance-covariance matrix associated with the regression coefficients.
#'
#'
vcov.ypbp <- function(object, ...){
  p <- object$p
  q <- object$q
  V <- MASS::ginv(-object$fit$hessian)[1:(2*q+p), 1:(2*q+p)]
  colnames(V) <- names(object$fit$par)[1:(2*q+p)]
  rownames(V) <- names(object$fit$par)[1:(2*q+p)]
  return(V)
}


#---------------------------------------------
#' Estimated regression coefficients
#'
#' @aliases coef.ypbp
#' @description This function returns the estimated regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class ypbp.
#' @param ... further arguments passed to or from other methods.
#' @return  the estimated regression coefficients.
#' @examples
#' \donttest{
#' fit <- ypbp(Surv(time, status)~arm, data=ipass)
#' coef(fit)
#' }
#'
coef.ypbp <- function(object, ...){
  p <- object$p
  q <- object$q
  coeffs <- object$fit$par[1:(2*q+p)]
  return(coeffs)
}


#---------------------------------------------
#' Generic S3 method confint
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the confidence intervals for the regression coefficients
#'
confint <- function(object, ...) UseMethod("confint")

#---------------------------------------------
#' Confidence intervals for the regression coefficients
#'
#' @aliases confint.ypbp
#' @description This function returns the estimated confidence intervals for the regression coefficients when the maximum likelihood estimation approach is used in the model fitting.
#' @export
#' @param object an object of the class ypbp.
#' @param level the confidence level required.
#' @param ... further arguments passed to or from other methods.
#' @return  A matrix (or vector) with columns giving lower and upper confidence limits for the regression coefficients. These will be labeled as (1-level)/2 and 1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
#' @examples
#' \donttest{
#' fit <- ypbp(Surv(time, status)~arm, data=ipass)
#' confint(fit)
#'}
#'
confint.ypbp <- function(object, level=0.95, ...){
  p <- object$p
  q <- object$q
  V <- vcov.ypbp(object)
  par.hat <- object$fit$par[1:(2*q+p)]
  alpha <- 1-level
  d <- stats::qnorm(1 - alpha/2)*sqrt(diag(V))
  lower <- par.hat - d
  upper <- par.hat + d
  CI <- cbind(lower, upper)
  labels <- round(100*(c(alpha/2, 1-alpha/2)),1)
  colnames(CI) <- paste0(labels, "%")
  return(CI)
}


#---------------------------------------------
#' Model.matrix method for ypbp models
#'
#' @aliases model.matrix.ypbp
#' @description Reconstruct the model matrix (or matrices if the alternative formulation of the YP model is used) for a ypbp model.
#' @export
#' @param object an object of the class ypbp.
#' @param ... further arguments passed to or from other methods.
#' @return  The model matrix (or matrices) for the fit.
#' @examples
#' \donttest{
#' fit <- ypbp(Surv(time, status)~arm, data=ipass)
#' model.matrix(fit)
#'}
#'
model.matrix.ypbp <- function(object, ...){
  formula <- Formula::Formula(object$formula)
  mf <- object$mf
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- suppressWarnings(try( stats::model.matrix(formula, data = mf, rhs = 2), TRUE))
  attrZ <- attributes(Z)
  attrZ$dim[2] <- ncol(Z) - 1
  attrZ$assign <- attrZ$assign[-1]
  attrZ$dimnames[[2]] <- attrZ$dimnames[[2]][-1]
  Z <- matrix(Z[,-1], ncol=ncol(Z) - 1)
  attributes(Z) <- attrZ
  if(ncol(X)>0){
    attrX <- attributes(X)
    attrX$dim[2] <- ncol(X) - 1
    attrX$assign <- attrX$assign[-1]
    attrX$dimnames[[2]] <- attrX$dimnames[[2]][-1]
    X <- matrix(X[,-1], ncol=ncol(X) - 1)
    attributes(X) <- attrX
  }
  if(ncol(X)>0){
    out <- list(Z=Z, X=X)
  }else{
    out <- Z
  }

  return(out)
}

