#' make predictions from a BVCfit object
#'
#' make predictions from a BVCfit object
#'
#' @param object BVCfit object.
#' @param X.new a matrix of new values for X at which predictions are to be made.
#' @param Z.new a vector of new values for Z at which predictions are to be made.
#' @param E.new a vector of new values for E at which predictions are to be made.
#' @param clin.new a vector or matrix of new values for clin at which predictions are to be made.
#' @param Y.new a vector of the response of new observations. If provided, the prediction mean squared error (PMSE) will be computed based on Y.new.
#' @param ... other predict arguments
#'
#' @details X.new (clin.new) must have the same number of columns as X (clin) used for fitting the model. If E and clin are provided when fit the model, E.new and clin.new
#' must not be NULL, and vice versa. The predictions are made based on the posterior estimates of coefficients in the BVCfit object.
#' Note that the main effects of environmental exposures Z and E are not subject to selection.
#'
#' @return  an object of class "BVCfit.pred" is returned, which is a list with components:
#' \item{pmse}{predictions mean squared error. pmse is NULL is Y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.BVCfit
#' @method predict BVCfit
#' @seealso \code{\link{BVCfit}}
#'
#' @examples
#' data(gExp)
#' spbayes=BVCfit(X, Y, Z, E, clin)
#' spbayes
#'
#' data(gExp.new)
#' pred = predict(spbayes, X.new, Z.new, E.new, clin.new, Y.new)
#' pred$pmse
#' # pred$y.pred
#'
#' @export
predict.BVCfit=function(object, X.new, Z.new, E.new=NULL, clin.new=NULL, Y.new = NULL, ...){
  # print(class(object))
  NextMethod()
}


#' @rdname predict.BVCfit
#' @method predict VarLin
#' @export
predict.VarLin=function(object, X.new, Z.new, E.new, clin.new=NULL, Y.new=NULL, ...){
  x = as.matrix(X.new)
  CLC=cbind(E.new, clin.new)
  EX = as.matrix(as.numeric(E.new) * x)
  x = cbind(1, x) # add intercept
  design = Design.matrix(Z.new, x, object$basis$kn, object$basis$degree)
  xx = design$Xns

  coeff = c(object$coefficient$ZX)
  coeff.clc = c(object$coefficient$E, object$coefficient$clin)
  coeff.zeta = object$coefficient$EX

  if(length(coeff)!=ncol(xx)){
    stop(paste("number of columns of X.new dose not match the length of the estimates."))
  }

  if(length(coeff.clc)!=ncol(CLC)){
    stop(paste("incorrect number of clinical covariates (", ncol(CLC)-1, "), supposed to be ", length(coeff.clc)-1, sep = ""))
  }

  y.pred = xx %*% coeff + CLC %*% coeff.clc + EX %*% coeff.zeta
  if(is.null(Y.new)){
    pmse = NULL
  }else{
    pmse = sum((Y.new - y.pred)^2)/length(Y.new)
  }

  pred = list(pmse=pmse, y.pred=y.pred)
  class(pred) = "BVCfit.pred"
  pred
}


#' @rdname predict.BVCfit
#' @method predict VarOnly
#' @export
predict.VarOnly=function(object, X.new, Z.new, clin.new=NULL, Y.new=NULL, ...)
{
  x = as.matrix(X.new)
  coeff = c(object$coefficient$ZX)
  coeff.clc = object$coefficient$clin

  if(is.null(clin.new)){
    if(!is.null(coeff.clc)) stop(paste("please provide the data of ", length(coeff.clc), " clinical covariates", sep = ""))

  }else{
    CLC=as.matrix(clin.new)
    if(is.null(coeff.clc) || length(coeff.clc)!=ncol(CLC)){
      stop(paste("incorrect number of clinical covariates (", ncol(CLC)-1, "), supposed to be ", length(coeff.clc), sep = ""))
    }
  }
  x = cbind(1, x) # add intercept
  design = Design.matrix(Z.new, x, object$basis$kn, object$basis$degree)
  xx = design$Xns

  if(length(coeff)!=ncol(xx)){
    stop(paste("number of columns of X.new dose not match the length of the estimates."))
  }

  y.pred = xx %*% coeff
  if(!is.null(coeff.clc)) y.pred = y.pred + CLC %*% coeff.clc

  if(is.null(Y.new)){
    pmse = NULL
  }else{
    pmse = sum((Y.new - y.pred)^2)/length(Y.new)
  }

  pred = list(pmse=pmse, y.pred=y.pred)
  class(pred) = "BVCfit.pred"
  pred
}


#' @rdname predict.BVCfit
#' @method predict LinOnly
#' @export
predict.LinOnly=function(object, X.new, Z.new, E.new=NULL, clin.new=NULL, Y.new=NULL, ...){
  x = as.matrix(X.new)
  CLC = cbind(E.new, clin.new, as.matrix(Z.new, ncol=1))
  ZX = as.matrix(as.numeric(Z.new) * x);
  xx = cbind(1, x) # add intercept

  coeff = c(object$coefficient$ZX$intercept, object$coefficient$ZX$Main)
  coeff.clc = c(object$coefficient$E, object$coefficient$clin, object$coefficient$ZX$Z)
  coeff.EX = object$coefficient$EX
  coeff.ZX = object$coefficient$ZX$Interaction

  if(length(coeff)!=ncol(xx)){
    stop(paste("number of columns of X.new dose not match the length of the estimates."))
  }

  if(is.null(E.new)){
    if(!is.null(coeff.EX)) stop(paste("please provide the data of environmental factor E", sep = ""))
    EX = NULL
  }else{
    if(is.null(coeff.EX)) stop(paste("environmental factor E should be NULL", sep = ""))
    EX = as.matrix(as.numeric(E.new) * x)
  }

  if(length(coeff.clc)!=ncol(CLC)){
    stop(paste("incorrect number of clinical covariates (", ncol(CLC)-(!is.null(E.new))-1, "), supposed to be ", length(coeff.clc)-(!is.null(object$coefficient$E))-1, sep = ""))
  }

  y.pred = xx %*% coeff + CLC %*% coeff.clc + ZX %*% coeff.ZX

  if(!is.null(EX) && !is.null(coeff.EX)) y.pred = y.pred + EX %*% coeff.EX

  if(is.null(Y.new)){
    pmse = NULL
  }else{
    pmse = sum((Y.new - y.pred)^2)/length(Y.new)
  }

  pred = list(pmse=pmse, y.pred=y.pred)
  class(pred) = "BVCfit.pred"
  pred
}
