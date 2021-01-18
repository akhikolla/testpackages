##' Joint modelling for longitutal and censored data with competing risks
##' @title Coefficients of longitudinal/survival sub-model
##' @param object  The JMcmprsk object returned by either jmo or jmc function.
##' @param coeff The coefficients returned by the JMcmprsk object.
##' @param ... further arguments passed to or from other methods.
##' @seealso \code{\link{jmc}}
##' @return Return estimates fixed effects with variable names
##' @export
coef.JMcmprsk <-
  function (object, coeff=c("beta", "alpha", "gamma"), ...) {
    if (!inherits(object, "JMcmprsk"))
      stop("Use only with 'JMcmprsk' objects.\n")
    if (object$type == "jmo") {
      if (coeff == "beta") {
        betas <- object$betas
        out=betas
        return(out)
      } else if (coeff == "alpha") {
        alphas <- object$alphamatrix
        out=alphas
        return(out)
      } else if (coeff == "gamma") {
        gammas <- object$gamma_matrix
        out=gammas
        return(out)
      } else if (coeff == "theta") {
        thetas <- object$thetas
        out=thetas
        return(out)
      } else {
        stop("Unexpected arguments! Must choose one of the following options: beta, alpha, gamma, theta.")
      }
    } else if (object$type == "jmc")  {
      if (coeff == "beta") {
        betas <- object$betas
        out=betas
        return(out)
      } else if (coeff == "gamma") {
        gammas <- object$gamma_matrix
        out=gammas
        return(out)
      } else {
        stop("Unexpected arguments! Must choose one of the following options: beta, gamma.")
      }

    }
}
