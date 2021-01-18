#' make predictions from a roben object
#'
#' make predictions from a roben object
#'
#' @param object roben object.
#' @param X.new a matrix of new values for X at which predictions are to be made.
#' @param E.new a vector of new values for E at which predictions are to be made.
#' @param clin.new a vector or matrix of new values for clin at which predictions are to be made.
#' @param Y.new a vector of the response of new observations. If provided, the prediction error will be computed based on Y.new.
#' @param ... other predict arguments
#'
#' @details X.new (E.new) must have the same number of columns as X (E) used for fitting the model. If clin was provided when fit the model, clin.new
#' must not be NULL, and vice versa. The predictions are made based on the posterior estimates of coefficients in the roben object.
#' Note that the main effects of environmental exposures E are not subject to selection.
#'
#' If Y.new is provided, the prediction error will be computed. For robust methods, the prediction mean absolute deviations (PMAD) will be computed.
#' For non-robust methods, the prediction mean squared error (PMSE) will be computed.
#'
#' @return  an object of class `roben.pred' is returned, which is a list with components:
#' \item{error}{prediction error. error is NULL is Y.new=NULL.}
#' \item{y.pred}{predicted values of the new observations.}
#'
#' @rdname predict.roben
#' @method predict roben
#' @seealso \code{\link{roben}}
#'
#' @examples
#' data(GxE_small)
#' test=sample((1:nrow(X)), floor(nrow(X)/5))
#' fit=roben(X[-test,], Y[-test,], E[-test,], clin[-test,], iterations=5000)
#' predict(fit, X[test,], E[test,], clin[test,], Y[test,])
#'
#' @export
predict.roben=function(object, X.new, E.new, clin.new=NULL, Y.new=NULL, ...){

  intercept = TRUE
  dat = Data.matrix(X.new, Y.new, E.new, clin.new, intercept)
  xx = dat$xx
  y.new = dat$y
  CLC = dat$CLC

  coeff = c(object$coefficient$GE)
  coeff.clc = c(object$coefficient$Int, object$coefficient$clin, object$coefficient$E)

  if(length(coeff)!=ncol(xx)){
    stop(paste("number of columns of X.new dose not match the length of the estimates."))
  }

  if(length(coeff.clc)!=ncol(CLC)){
    stop(paste("incorrect number of clinical covariates (", ncol(CLC)-1, "), supposed to be ", length(coeff.clc)-1, sep = ""))
  }

  y.pred = xx %*% coeff + CLC %*% coeff.clc
  error = NULL

  if(inherits(object, "RBVS")){
    error = sum(abs(y.new - y.pred))/length(y.new)
    # error.type = "PMAD"
    names(error) = "PMAD"
  }else{
    error = sum((y.new - y.pred)^2)/length(y.new)
    # error.type = "PMSE"
    names(error) = "PMSE"
  }

  pred = list(error=error, y.pred=y.pred)
  class(pred) = "roben.pred"
  pred
}


